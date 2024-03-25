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

	int ctr = 1, loc = CHL_MEMLOCS -1, elem_size = 8, case_id = CHL_WORKERS, log_results = 0;

	switch (argc) {
	case (4):
		loc = atoi(argv[ctr++]);
		case_id = atoi(argv[ctr++]);
		log_results = atoi(argv[ctr++]);
		break;
	default:
		error("Incorrect input arguments. Usage: ./correct_run loc case_id log_results:\n"
		"loc: the src memory location for the broadcast\n"
		"case_id: the broadcast destinations\n"
		"log_results: If !=0 log results to file \n");
  	}

	if(loc < 0 || loc >= CHL_MEMLOCS) error("bw_bench_broadcast_2D: Unsupported src loc = %d\n", loc);
	char *filename;
	FILE* fp;
	if(log_results){
		filename = (char *) malloc(1024 * sizeof(char));
    	sprintf(filename, "%s/Database/microbenchmarks/bw_bench_broadcast_2D_%d_%d.log", DEPLOYDB, loc, case_id);
		fp = fopen(filename, "r");
    	if(fp){
			warning("bw_bench_broadcast_2D: filename %s exists, quiting\n", filename);
			return 1;
		}
		fp = fopen(filename, "w+");
    	if(!fp) error("bw_bench_broadcast_2D: File path %s is wrong or write permission missing\n", filename);
	}
	int active_unit_num = 0, active_unit_id_list[CHL_WORKERS];
	translate_binary_to_unit_list(case_id, &active_unit_num, active_unit_id_list);
	int maxDim = std::min(MAX_DIM_TRANS, (int) CHLGetMaxDimSqAsset2D(active_unit_num, elem_size, STEP_TRANS, loc));
	CHLEnableLinks(loc, active_unit_num);
	for(int dev_id_idx = 0 ; dev_id_idx < active_unit_num; dev_id_idx++){
		maxDim = std::min(maxDim, (int) CHLGetMaxDimSqAsset2D(1, elem_size, STEP_TRANS, (dev_id_idx)));
	}
	long long ldim = maxDim;
	fprintf(stderr,"\nbw_bench_broadcast_2D: \nSystem = %s\nmaxDim = %d, ldim = %lld\n", 
		TESTBED, maxDim, ldim);
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
	for (int it = 0; it < 3; it++){
		fprintf(stderr, ".");
		for(int dev_id_idx = 0 ; dev_id_idx < active_unit_num; dev_id_idx++){
			int loc_dest = active_unit_id_list[dev_id_idx]; 
			if(loc!=loc_dest) queue_list[dev_id_idx]->memcpy2DAsync(worker_buffs[dev_id_idx], ldim,
				loc_buffs[dev_id_idx], ldim, maxDim, maxDim, elemSize, loc_dest, loc, 1);
			queue_list[dev_id_idx]->sync_barrier();
		}
	}
	fprintf(stderr, " complete.\n-------------------------------------------------------------------------------"
		"-----------------------------------------------------------------------\n");
	CHLSyncCheckErr();
	
	for (int dim = 256; dim <= maxDim; dim*=2){
		int sample_sz;
		double dev_t[active_unit_num], transfer_t_vals[active_unit_num][MICRO_MAX_ITER] = {0}, 
			transfer_t_sum[active_unit_num] = {0}, transfer_t_mean[active_unit_num] = {0}, 
			error_margin[active_unit_num] = {0}, bench_t = csecond();
		for (sample_sz = 1; sample_sz < MICRO_MAX_ITER + 1; sample_sz++){
			for(int dev_id_idx = 0; dev_id_idx < active_unit_num; dev_id_idx++)
			if(loc != active_unit_id_list[dev_id_idx]){
				int loc_dest = active_unit_id_list[dev_id_idx]; 
				int queue_id = (loc >= CHL_WORKERS || loc < 0)? loc_dest : loc;
				CHLSelectDevice(queue_id);
				device_timer[dev_id_idx]->start_point(queue_list[dev_id_idx]);
				queue_list[dev_id_idx]->memcpy2DAsync(worker_buffs[dev_id_idx], ldim,
					loc_buffs[dev_id_idx], ldim, dim, dim, elemSize, loc_dest, loc, 1);
				device_timer[dev_id_idx]->stop_point(queue_list[dev_id_idx]);
			}
			for(int dev_id_idx = 0; dev_id_idx < active_unit_num; dev_id_idx++) 
			if(loc != active_unit_id_list[dev_id_idx])
				queue_list[dev_id_idx]->sync_barrier();
			for(int dev_id_idx = 0; dev_id_idx < active_unit_num; dev_id_idx++) 
			if(loc != active_unit_id_list[dev_id_idx]){
				int loc_dest = active_unit_id_list[dev_id_idx]; 
				int queue_id = (loc >= CHL_WORKERS || loc < 0)? loc_dest : loc;
				CHLSelectDevice(queue_id);
				dev_t[dev_id_idx] = device_timer[dev_id_idx]->sync_get_time()/1000;
			}
			CHLSyncCheckErr();
			int complete_flag = 1;
			for(int dev_id_idx = 0; dev_id_idx < active_unit_num; dev_id_idx++) 
			if(loc != active_unit_id_list[dev_id_idx]){
				if(!confidence_interval_5_percent(sample_sz, dev_t[dev_id_idx], transfer_t_vals[dev_id_idx], &transfer_t_sum[dev_id_idx], 
					&transfer_t_mean[dev_id_idx], &error_margin[dev_id_idx])) complete_flag = 0;
			}
			if (sample_sz > MICRO_MIN_ITER && complete_flag) break;
		}
		double broadcast_bw[active_unit_num] = {0}, sum_bw = 0;
		for(int dev_id_idx = 0; dev_id_idx < active_unit_num; dev_id_idx++) 
			if(loc != active_unit_id_list[dev_id_idx]) sum_bw += broadcast_bw[dev_id_idx] = Gval_per_s(dim*dim*elemSize, transfer_t_mean[dev_id_idx]);
			else broadcast_bw[dev_id_idx] = -1; 
//#ifdef PDEBUG
		fprintf(stderr, "Ran %d itterations for convergence.\n"
			"-> dim = %d, active_unit_id_list = %s :\n\t BW(link) = %s Gb/s\n\t BW(sum) = %.2lf Gb/s\n",
			sample_sz, dim, printlist(active_unit_id_list, active_unit_num), printlist(broadcast_bw, active_unit_num), sum_bw);
		fprintf(stderr,"-------------------------------------------------------------------------------"
			"-----------------------------------------------------------------------\n");
//#endif
		if(log_results) for(int dev_id_idx = 0; dev_id_idx < active_unit_num; dev_id_idx++)
			fprintf(fp, "%d,%d,%d,%d,%lf\n", active_unit_id_list[dev_id_idx], dim, dim, elemSize, broadcast_bw[dev_id_idx]);
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
