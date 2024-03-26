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

	int ctr = 1, loc = CHL_MEMLOCS -1, elemSize = 8, source_locs_binary = CHL_WORKERS, log_results = 0;

	switch (argc) {
	case (4):
		loc = atoi(argv[ctr++]);
		source_locs_binary = atoi(argv[ctr++]);
		log_results = atoi(argv[ctr++]);
		break;
	default:
		error("Incorrect input arguments. Usage: ./correct_run gather_loc source_locs_binary log_results:\n"
		"loc: the dest memory location for the gather\n"
		"source_locs_binary: the gather sources\n"
		"log_results: If !=0 log results to file \n");
  	}

	if(loc < 0 || loc >= CHL_MEMLOCS) 
		error("bw_bench_gather_2D: Unsupported src loc = %d\n", loc);
	char *filename;
	FILE* fp;
	if(log_results){
		filename = (char *) malloc(1024 * sizeof(char));
    	sprintf(filename, "%s/Database/microbenchmarks/bw_bench_gather_2D_%d_%d.log", DEPLOYDB, loc, source_locs_binary);
		fp = fopen(filename, "r");
    	if(fp){
			warning("bw_bench_gather_2D: filename %s exists, quiting\n", filename);
			return 1;
		}
		fp = fopen(filename, "w+");
    	if(!fp) error("bw_bench_gather_2D: File path %s is wrong or write permission missing\n", filename);
	}
	int active_unit_num = 0, active_unit_id_list[CHL_WORKERS];
	translate_binary_to_unit_list(source_locs_binary, &active_unit_num, active_unit_id_list);
	int maxDim = CHLGetMaxDimSqAsset2D(active_unit_num, elemSize, TILE_MAX, loc);
	CHLEnableLinks(loc, CHL_WORKERS);
	for(int dev_id_idx = 0 ; dev_id_idx < active_unit_num; dev_id_idx++){
		int loc_dest = active_unit_id_list[dev_id_idx]; 
		maxDim = std::min(maxDim, (int) CHLGetMaxDimSqAsset2D(4, elemSize, TILE_MAX, loc_dest));
	}
	long long ldim = maxDim;
	fprintf(stderr,"\nbw_bench_gather_2D: \nSystem = %s\nmaxDim = %d, ldim = %lld src loc = %d, dest = %s\n", 
		TESTBED, maxDim, ldim, loc, printlist(active_unit_id_list, active_unit_num));
	fprintf(stderr,"-------------------------------------------------------------------------------"
		"-----------------------------------------------------------------------\n");
	double timer = csecond();
	void* loc_buffs[active_unit_num], *worker_buffs[active_unit_num];
	for(int dev_id_idx = 0 ; dev_id_idx < active_unit_num; dev_id_idx++){
		int loc_dest = active_unit_id_list[dev_id_idx]; 
		loc_buffs[dev_id_idx] = CHLMalloc(ldim*ldim*elemSize, loc, 1);
		worker_buffs[dev_id_idx] = CHLMalloc(ldim*ldim*elemSize, loc_dest, 1);
	}
	timer = csecond() - timer;
	fprintf(stderr, "Allocation buffer size = (%lld x %lld) x %d complete:\t alloc_timer=%lf ms\n", ldim, ldim, elemSize, timer  * 1000);
	fprintf(stderr, "-------------------------------------------------------------------------------"
		"-----------------------------------------------------------------------\n");
	CQueue_p queue_list[active_unit_num];
	Event_timer_p device_timer[active_unit_num];
	for(int dev_id_idx = 0 ; dev_id_idx < active_unit_num; dev_id_idx++){
		//printf("dev_id = %d, dev_id_idx = %d, dev_id_idy = %d, CHL_WORKERS = %d\n", dev_id, dev_id_idx, dev_id_idy, CHL_WORKERS);
		int loc_dest = active_unit_id_list[dev_id_idx]; 
		int queue_id = (loc >= CHL_WORKERS || loc < 0)? loc_dest : loc;
		queue_list[dev_id_idx] = new CommandQueue(queue_id, COMMUNICATION);
		device_timer[dev_id_idx] = new Event_timer(queue_id);
	}

	fprintf(stderr, "Warming up");
	/// Warmup.
	for (int it = 0; it < 10; it++){
		fprintf(stderr, ".");
		for(int dev_id_idx = 0 ; dev_id_idx < active_unit_num; dev_id_idx++){
			int loc_dest = active_unit_id_list[dev_id_idx]; 
			if(loc!=loc_dest) queue_list[dev_id_idx]->memcpy2DAsync(worker_buffs[dev_id_idx], ldim,
				loc_buffs[dev_id_idx], ldim, TILE_MAX, TILE_MAX, elemSize, loc_dest, loc, 1);
			queue_list[dev_id_idx]->sync_barrier();
		}
	}
	fprintf(stderr, " complete.\n-------------------------------------------------------------------------------"
		"-----------------------------------------------------------------------\n");
	CHLSyncCheckErr();
	int chunk_dim_num = maxDim/TILE_MAX;
	for (int dim = TILE_MAX; dim <= TILE_MAX; dim*=2){
		int sample_sz;
		double dev_t[active_unit_num], transfer_t_vals[active_unit_num][MICRO_MAX_ITER] = {0}, 
			transfer_t_sum[active_unit_num] = {0}, transfer_t_mean[active_unit_num] = {0}, 
			error_margin[active_unit_num] = {0}, bench_t = 0;
		for (sample_sz = 1; sample_sz < MICRO_MAX_ITER + 1; sample_sz++){
			timer = csecond();
			int runctr = 0, case_ran[active_unit_num] = {0};
			while (runctr < active_unit_num){
				int dev_id_idx = ((int) rand()) % active_unit_num;
				if (case_ran[dev_id_idx]) continue;
				else{
					if(loc != active_unit_id_list[dev_id_idx]){
						int loc_dest = active_unit_id_list[dev_id_idx]; 
						int queue_id = (loc >= CHL_WORKERS || loc < 0)? loc_dest : loc;
						CHLSelectDevice(queue_id);
						device_timer[dev_id_idx]->start_point(queue_list[dev_id_idx]);
						for(long d1 = 0; d1< chunk_dim_num; d1++)
							for(long d2 = 0; d2< chunk_dim_num; d2++){
								long addroffset = d1*dim + d2*dim*ldim;
								queue_list[dev_id_idx]->memcpy2DAsync(worker_buffs[dev_id_idx] + addroffset, ldim,
								loc_buffs[dev_id_idx] + addroffset, ldim, dim, dim, elemSize, loc_dest, loc, 1);
							}
						device_timer[dev_id_idx]->stop_point(queue_list[dev_id_idx]);
					}
					case_ran[dev_id_idx] = 1;
					runctr++;
				}
			}
			for(int dev_id_idx = 0; dev_id_idx < active_unit_num; dev_id_idx++)
			
			for(int dev_id_idx = 0; dev_id_idx < active_unit_num; dev_id_idx++) 
			if(loc != active_unit_id_list[dev_id_idx])
				queue_list[dev_id_idx]->sync_barrier();
			for(int dev_id_idx = 0; dev_id_idx < active_unit_num; dev_id_idx++) 
			if(loc != active_unit_id_list[dev_id_idx]){
				int loc_dest = active_unit_id_list[dev_id_idx]; 
				int queue_id = (loc >= CHL_WORKERS || loc < 0)? loc_dest : loc;
				CHLSelectDevice(queue_id);
				dev_t[dev_id_idx] = device_timer[dev_id_idx]->sync_get_time()/1000/(chunk_dim_num*chunk_dim_num);
			}
			CHLSyncCheckErr();
			int complete_flag = 1;
			for(int dev_id_idx = 0; dev_id_idx < active_unit_num; dev_id_idx++) 
			if(loc != active_unit_id_list[dev_id_idx]){
				if(!confidence_interval_5_percent(sample_sz, dev_t[dev_id_idx], transfer_t_vals[dev_id_idx], &transfer_t_sum[dev_id_idx], 
					&transfer_t_mean[dev_id_idx], &error_margin[dev_id_idx])) complete_flag = 0;
			}
			if (sample_sz > MICRO_MIN_ITER && complete_flag) break;
			timer = csecond() - timer;
			bench_t += timer;
			if(bench_t > 60){
				fprintf(stderr, "Microbench itter ran %lf sec ( > 1 min): Stopping sampling at %d/%d\n", 
				bench_t, sample_sz, MICRO_MAX_ITER);
				break;
			}		
		}
		double gather_bw[active_unit_num] = {0}, sum_bw = 0;
		for(int dev_id_idx = 0; dev_id_idx < active_unit_num; dev_id_idx++) 
			if(loc != active_unit_id_list[dev_id_idx]) sum_bw += gather_bw[dev_id_idx] = Gval_per_s(dim*dim*elemSize, transfer_t_mean[dev_id_idx]);
			else gather_bw[dev_id_idx] = -1; 
//#ifdef PDEBUG
		fprintf(stderr, "Ran %d itterations for convergence (bench_t = %.3lf s)\n"
			"-> dim = %d (chunk_dim_num = %d), active_unit_id_list = %s :\n\t BW(link) = %s Gb/s\n\t BW(sum) = %.2lf Gb/s\n",
			sample_sz, bench_t, dim, chunk_dim_num, printlist(active_unit_id_list, active_unit_num), printlist(gather_bw, active_unit_num), sum_bw);
		fprintf(stderr,"-------------------------------------------------------------------------------"
			"-----------------------------------------------------------------------\n");
//#endif
		if(log_results) for(int dev_id_idx = 0; dev_id_idx < active_unit_num; dev_id_idx++)
			fprintf(fp, "%d,%d,%d,%d,%lf\n", active_unit_id_list[dev_id_idx], dim, dim, elemSize, gather_bw[dev_id_idx]);
	}
	timer = csecond();
	for(int dev_id_idx = 0; dev_id_idx < active_unit_num; dev_id_idx++){
			int loc_dest = active_unit_id_list[dev_id_idx];
  			CHLFree(loc_buffs[dev_id_idx], ldim*ldim*elemSize, loc);
  			CHLFree(worker_buffs[dev_id_idx], ldim*ldim*elemSize, loc_dest);
	}
	timer = csecond() - timer;
	fprintf(stderr, "Free buffers = (%lld x %lld) x %d complete:\t alloc_timer=%lf ms\n", ldim, ldim, elemSize, timer  * 1000);
	fprintf(stderr,"-------------------------------------------------------------------------------"
		"-----------------------------------------------------------------------\n");
	if(log_results) fclose(fp);
  	return 0;
}
