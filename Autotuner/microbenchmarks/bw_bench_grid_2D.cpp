///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief A transfer microbenchmark aiming to capture how link overlap effects their bandwidth
///

#include <unistd.h>
#include <cassert>

#include <numa.h>

#include "chl_smart_wrappers.hpp"
#include "chl_grid_amalgamation.hpp"
#include "microbenchmarks.hpp"

int main(const int argc, const char *argv[]) {

	int ctr = 1, host_loc = CHL_MEMLOCS -1, elemSize = 8, 
		worker_case_id, h_case_id, h_rev_case_id = worker_case_id = h_case_id = CHL_WORKERS, log_results = 0;

	switch (argc) {
	case (6):
		host_loc = atoi(argv[ctr++]);
		worker_case_id = atoi(argv[ctr++]);
		h_case_id = atoi(argv[ctr++]);
		h_rev_case_id = atoi(argv[ctr++]);
		log_results = atoi(argv[ctr++]);
		break;
	default:
		error("Incorrect input arguments. Usage: ./correct_run host_loc worker_case_id h_case_id h_rev_case_id log_results:\n"
		"host_loc: the host memory location that will be used for the benchmark\n"
		"worker_case_id: The list of active workers (destinations)\n"
		"h_case_id: Copy queues for host -> dev \n"
		"h_rev_case_id: Copy queues for dev -> host \n"
		"log_results: If !=0 log results to file \n");
  	}

	if(host_loc < CHL_WORKERS || host_loc >= CHL_MEMLOCS) 
		error("bw_bench_bidirectional_grid_2D: Unsupported host_loc = %d\n", host_loc);
	if(!is_subset(h_case_id, worker_case_id)) error("bw_bench_bidirectional_grid_2D: "
		"h_case_id = %d is not a subset of worker_case_id = %d\n", h_case_id, worker_case_id);
	if(!is_subset(h_rev_case_id, worker_case_id)) error("bw_bench_bidirectional_grid_2D: "
		"h_rev_case_id = %d is not a subset of worker_case_id = %d\n", h_rev_case_id, worker_case_id);

	char *filename;
	FILE* fp;
	if(log_results){
		filename = (char *) malloc(1024 * sizeof(char));
    	sprintf(filename, "%s/Database/microbenchmarks/bw_bench_bidirectional_grid_2D_%d_%d_%d_%d.log", 
			DEPLOYDB, host_loc, worker_case_id, h_case_id, h_rev_case_id);
		fp = fopen(filename, "r");
    	if(fp){
			warning("bw_bench_bidirectional_grid_2D: filename %s exists, quiting\n", filename);
			return 1;
		}
		fp = fopen(filename, "w+");
    	if(!fp) error("bw_bench_bidirectional_grid_2D: File path %s is wrong or write permission missing\n", filename);
	}
	int active_unit_num = 0, active_unit_id_list[CHL_WORKERS];
	translate_binary_to_unit_list(worker_case_id, &active_unit_num, active_unit_id_list);
	int active_memloc_num = active_unit_num + 1, active_memloc_id_list[active_memloc_num];
	for(int dev_id_idx = 0 ; dev_id_idx < active_unit_num; dev_id_idx++)
		active_memloc_id_list[dev_id_idx] = active_unit_id_list[dev_id_idx];
	active_memloc_id_list[active_memloc_num-1] = host_loc;

	int active_h2d_queues = 0, active_h2d_queue_ids[CHL_WORKERS];
	translate_binary_to_unit_list(h_case_id, &active_h2d_queues, active_h2d_queue_ids);
	int active_d2h_queues = 0, active_d2h_queue_ids[CHL_WORKERS];
	translate_binary_to_unit_list(h_rev_case_id, &active_d2h_queues, active_d2h_queue_ids);

	int load_mult_h2d = active_unit_num/active_h2d_queues, load_mult_d2h = active_unit_num/active_d2h_queues;

	int maxDim = CHLGetMaxDimSqAsset2D(2*active_unit_num, elemSize, TILE_MAX, host_loc);
	for(int dev_id_idx = 0 ; dev_id_idx < active_unit_num; dev_id_idx++) 
		CHLEnableLinks(active_unit_id_list[dev_id_idx], active_unit_num);
	for(int dev_id_idx = 0 ; dev_id_idx < active_unit_num; dev_id_idx++){
		maxDim = std::min(maxDim, (int) CHLGetMaxDimSqAsset2D(2*active_unit_num, elemSize, TILE_MAX, active_unit_id_list[dev_id_idx]));
	}
	long long ldim = maxDim;
	fprintf(stderr,"\nbw_bench_bidirectional_grid_2D: \nSystem = %s\nmaxDim = %d, ldim = %lld\n", 
		TESTBED, maxDim, ldim);
	fprintf(stderr,"-------------------------------------------------------------------------------"
		"-----------------------------------------------------------------------\n");

	double timer = csecond();
	void* loc_buffs[active_memloc_num][active_memloc_num][2];
	for(int dev_id_idx = 0 ; dev_id_idx < active_memloc_num; dev_id_idx++){
		for(int dev_id_idy = 0 ; dev_id_idy < active_memloc_num; dev_id_idy++){
			int loc_src = active_memloc_id_list[dev_id_idy];
			int loc_dest = active_memloc_id_list[dev_id_idx];
			if(loc_src == loc_dest){
				loc_buffs[dev_id_idx][dev_id_idy][0] = loc_buffs[dev_id_idx][dev_id_idy][1] = NULL;
				continue;
			}
			else if((loc_src == host_loc && is_in_list(loc_dest, active_h2d_queue_ids, active_h2d_queues)) ||
				(loc_dest == host_loc && is_in_list(loc_src, active_d2h_queue_ids, active_d2h_queues)) ||
				(loc_src < CHL_WORKERS && loc_dest < CHL_WORKERS)){
				loc_buffs[dev_id_idx][dev_id_idy][0] = CHLMalloc(ldim*ldim*elemSize, loc_dest, 1);
				if (loc_dest == CHL_WORKERS) CHLTouche((double*) loc_buffs[dev_id_idx][dev_id_idy][0], ldim*ldim, elemSize);
				loc_buffs[dev_id_idx][dev_id_idy][1] = CHLMalloc(ldim*ldim*elemSize, loc_src, 1);
				if (loc_src == CHL_WORKERS) CHLTouche((double*) loc_buffs[dev_id_idx][dev_id_idy][1], ldim*ldim, elemSize);
			}
			else loc_buffs[dev_id_idx][dev_id_idy][0] = loc_buffs[dev_id_idx][dev_id_idy][1] = NULL;
		}
	}
	timer = csecond() - timer;
	fprintf(stderr, "Allocation buffer size = (%lld x %lld) x %d complete:\t alloc_timer=%lf ms\n", ldim, ldim, elemSize, timer  * 1000);
	fprintf(stderr, "-------------------------------------------------------------------------------"
		"-----------------------------------------------------------------------\n");
	CQueue_p queue_list[active_memloc_num][active_memloc_num];
	Event_timer_p device_timer[active_memloc_num][active_memloc_num];
	for(int dev_id_idx = 0 ; dev_id_idx < active_memloc_num; dev_id_idx++){
		for(int dev_id_idy = 0 ; dev_id_idy < active_memloc_num; dev_id_idy++){
		//printf("dev_id = %d, dev_id_idx = %d, dev_id_idy = %d, CHL_WORKERS = %d\n", dev_id, dev_id_idx, dev_id_idy, CHL_WORKERS);
			int loc_src = active_memloc_id_list[dev_id_idy];
			int loc_dest = active_memloc_id_list[dev_id_idx]; 
			int queue_id = (loc_src >= CHL_WORKERS || loc_src < 0)? loc_dest : loc_src;
			queue_list[dev_id_idx][dev_id_idy] = new CommandQueue(queue_id, COMMUNICATION);
			device_timer[dev_id_idx][dev_id_idy] = new Event_timer(queue_id);
		}
	}

	fprintf(stderr, "Warming up");
	/// Warmup.
	for (int it = 0; it < 10; it++){
		fprintf(stderr, ".");
		for(int dev_id_idx = 0 ; dev_id_idx < active_memloc_num; dev_id_idx++){
			for(int dev_id_idy = 0 ; dev_id_idy < active_memloc_num; dev_id_idy++){
				int loc_src = active_memloc_id_list[dev_id_idy];
				int loc_dest = active_memloc_id_list[dev_id_idx]; 
				if(loc_buffs[dev_id_idx][dev_id_idy][0]) queue_list[dev_id_idx][dev_id_idy]->memcpy2DAsync(
					loc_buffs[dev_id_idx][dev_id_idy][0], ldim, loc_buffs[dev_id_idx][dev_id_idy][1], 
					ldim, TILE_MAX, TILE_MAX, elemSize, loc_dest, loc_src, 1);
				queue_list[dev_id_idx][dev_id_idy]->sync_barrier();
			}
		}
	}
	fprintf(stderr, " complete.\n-------------------------------------------------------------------------------"
		"-----------------------------------------------------------------------\n");
	CHLSyncCheckErr();
	int chunk_dim_num = maxDim/TILE_MAX;
	for (int dim = TILE_MAX; dim <= TILE_MAX; dim*=2){
		int sample_sz;
		double dev_t[active_memloc_num][active_memloc_num], transfer_t_vals[active_memloc_num][active_memloc_num][MICRO_MAX_ITER], 
			transfer_t_sum[active_memloc_num][active_memloc_num], transfer_t_mean[active_memloc_num][active_memloc_num], 
			error_margin[active_memloc_num][active_memloc_num], bench_t = 0;
		for(int dev_id_idx = 0; dev_id_idx < active_memloc_num; dev_id_idx++)
				for(int dev_id_idy = 0 ; dev_id_idy < active_memloc_num; dev_id_idy++){
					dev_t[dev_id_idx][dev_id_idy] = transfer_t_sum[dev_id_idx][dev_id_idy] = 
					transfer_t_mean[dev_id_idx][dev_id_idy] = error_margin[dev_id_idx][dev_id_idy] = 0;
					for(int idx = 0 ; idx < MICRO_MAX_ITER; idx++) transfer_t_vals[dev_id_idx][dev_id_idy][idx] = 0;
				}
		for (sample_sz = 1; sample_sz < MICRO_MAX_ITER + 1; sample_sz++){
			timer = csecond();
			int runctr = 0, case_ran[active_memloc_num][active_memloc_num];
			for(int dev_id_idx = 0; dev_id_idx < active_memloc_num; dev_id_idx++)
				for(int dev_id_idy = 0 ; dev_id_idy < active_memloc_num; dev_id_idy++) case_ran[dev_id_idx][dev_id_idy] = 0;
			while (runctr < active_memloc_num*active_memloc_num){
				int dev_id_idx = ((int) rand()) % active_memloc_num;
				int dev_id_idy = ((int) rand()) % active_memloc_num;
				if (case_ran[dev_id_idx][dev_id_idy]) continue;
				else{
					if(loc_buffs[dev_id_idx][dev_id_idy][0]){
						int loc_src = active_memloc_id_list[dev_id_idy];
						int loc_dest = active_memloc_id_list[dev_id_idx]; 
						int queue_id = (loc_src >= CHL_WORKERS || loc_src < 0)? loc_dest : loc_src;
						CHLSelectDevice(queue_id);
						int load_itter = 1; 
						//if(loc_src == host_loc) load_itter = load_mult_h2d;
						//else if (loc_dest == host_loc) load_itter = load_mult_d2h;
						device_timer[dev_id_idx][dev_id_idy]->start_point(queue_list[dev_id_idx][dev_id_idy]);
						for(long d1 = 0; d1< chunk_dim_num; d1++)
							for(long d2 = 0; d2< chunk_dim_num; d2++){
								long addroffset = d1*dim + d2*dim*ldim;
								queue_list[dev_id_idx][dev_id_idy]->memcpy2DAsync(loc_buffs[dev_id_idx][dev_id_idy][0] + addroffset, ldim,
								loc_buffs[dev_id_idx][dev_id_idy][1] + addroffset, ldim, dim, dim, elemSize, loc_dest, loc_src, 1);
							}
						device_timer[dev_id_idx][dev_id_idy]->stop_point(queue_list[dev_id_idx][dev_id_idy]);
					}
					case_ran[dev_id_idx][dev_id_idy] = 1;
					runctr++;
				}
			}
			for(int dev_id_idx = 0; dev_id_idx < active_memloc_num; dev_id_idx++)
				for(int dev_id_idy = 0 ; dev_id_idy < active_memloc_num; dev_id_idy++)
					if(loc_buffs[dev_id_idx][dev_id_idy][0]) queue_list[dev_id_idx][dev_id_idy]->sync_barrier();
			for(int dev_id_idx = 0; dev_id_idx < active_memloc_num; dev_id_idx++){
				for(int dev_id_idy = 0 ; dev_id_idy < active_memloc_num; dev_id_idy++){
					if(loc_buffs[dev_id_idx][dev_id_idy][0]){
						int loc_src = active_memloc_id_list[dev_id_idy];
						int loc_dest = active_memloc_id_list[dev_id_idx]; 
						int queue_id = (loc_src >= CHL_WORKERS || loc_src < 0)? loc_dest : loc_src;
						dev_t[dev_id_idx][dev_id_idy] = device_timer[dev_id_idx][dev_id_idy]->sync_get_time()/1000/(chunk_dim_num*chunk_dim_num);
					}
				}
			}
			CHLSyncCheckErr();
			int complete_flag = 1;
			for(int dev_id_idx = 0; dev_id_idx < active_memloc_num; dev_id_idx++)
				for(int dev_id_idy = 0 ; dev_id_idy < active_memloc_num; dev_id_idy++)
					if(loc_buffs[dev_id_idx][dev_id_idy][0])
						if(!confidence_interval_5_percent(sample_sz, dev_t[dev_id_idx][dev_id_idy], transfer_t_vals[dev_id_idx][dev_id_idy], &transfer_t_sum[dev_id_idx][dev_id_idy], 
							&transfer_t_mean[dev_id_idx][dev_id_idy], &error_margin[dev_id_idx][dev_id_idy])) complete_flag = 0;
			if (sample_sz > MICRO_MIN_ITER && complete_flag) break;
			timer = csecond() - timer;
			bench_t += timer;
			if(bench_t > 60){
				fprintf(stderr, "Microbench itter ran %lf sec ( > 1 min): Stopping sampling at %d/%d\n", 
				bench_t, sample_sz, MICRO_MAX_ITER);
				break;
			}	
		}
		double grid_bw[active_memloc_num][active_memloc_num], sum_bw = 0;
		for(int dev_id_idx = 0; dev_id_idx < active_memloc_num; dev_id_idx++)
			for(int dev_id_idy = 0 ; dev_id_idy < active_memloc_num; dev_id_idy++)
				if(loc_buffs[dev_id_idx][dev_id_idy][0]){
					int loc_src = active_memloc_id_list[dev_id_idy];
					int loc_dest = active_memloc_id_list[dev_id_idx]; 
					int load_itter = 1; 
					//if(loc_src == host_loc) load_itter = load_mult_h2d;
					//else if (loc_dest == host_loc) load_itter = load_mult_d2h;
					sum_bw += grid_bw[dev_id_idx][dev_id_idy] = Gval_per_s(load_itter*dim*dim*elemSize, transfer_t_mean[dev_id_idx][dev_id_idy]);

				}
				else grid_bw[dev_id_idx][dev_id_idy] = -1; 
//#ifdef PDEBUG
		fprintf(stderr, "Ran %d itterations for convergence (bench_t = %.3lf s)\n"
			"-> dim = %d, active_memloc_id_list = %s, active_h2d_queue_ids = %s, active_d2h_queue_ids = %s:"
			"\n\tBW(sum) = %.2lf Gb/s\n", sample_sz, bench_t, dim, printlist(active_memloc_id_list, active_memloc_num),
			printlist(active_h2d_queue_ids, active_h2d_queues), printlist(active_d2h_queue_ids, active_d2h_queues), sum_bw);
		for(int dev_id_idx = 0; dev_id_idx < active_memloc_num; dev_id_idx++)
			fprintf(stderr, "\tdest(%d) : %s\n",  active_memloc_id_list[dev_id_idx], printlist(grid_bw[dev_id_idx], active_memloc_num));
		fprintf(stderr,"-------------------------------------------------------------------------------"
			"-----------------------------------------------------------------------\n");
//#endif
		if(log_results) for(int dev_id_idx = 0; dev_id_idx < active_memloc_num; dev_id_idx++)
			for(int dev_id_idy = 0 ; dev_id_idy < active_memloc_num; dev_id_idy++)
			fprintf(fp, "%d,%d, %d,%d,%d,%lf\n", active_memloc_id_list[dev_id_idx], active_memloc_id_list[dev_id_idy],  dim, dim, elemSize, grid_bw[dev_id_idx][dev_id_idy]);
	}
	timer = csecond();
	for(int dev_id_idx = 0; dev_id_idx < active_memloc_num; dev_id_idx++)
		for(int dev_id_idy = 0 ; dev_id_idy < active_memloc_num; dev_id_idy++)
			if(loc_buffs[dev_id_idx][dev_id_idy][0]){
				int loc_src = active_memloc_id_list[dev_id_idy];
				int loc_dest = active_memloc_id_list[dev_id_idx]; 
				CHLFree(loc_buffs[dev_id_idx][dev_id_idy][0],ldim*ldim*elemSize, loc_dest);
				CHLFree(loc_buffs[dev_id_idx][dev_id_idy][1],ldim*ldim*elemSize, loc_src);
			}
	timer = csecond() - timer;
	fprintf(stderr, "Free buffers = (%lld x %lld) x %d complete:\t alloc_timer=%lf ms\n", ldim, ldim, elemSize, timer  * 1000);
	fprintf(stderr,"-------------------------------------------------------------------------------"
		"-----------------------------------------------------------------------\n");
	if(log_results) fclose(fp);
  	return 0;
}
