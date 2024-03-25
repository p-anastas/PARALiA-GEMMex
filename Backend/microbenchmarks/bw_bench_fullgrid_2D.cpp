///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief A transfer microbenchmark aiming to capture how link overlap effects their bandwidth
///

#include <unistd.h>
#include <cassert>

#include <numa.h>

#include "smart_wrappers.hpp"
#include "grid_amalgamation.hpp"
#include "microbenchmarks.hpp"

int main(const int argc, const char *argv[]) {

	int ctr = 1, loc = CHL_MEMLOCS -1, elem_size = 8;

	switch (argc) {
	case (2):
		loc = atoi(argv[ctr++]);
		break;
	default:
		error("Incorrect input arguments. Usage: ./correct_run host_loc h_case_id h_rev_case_id:\n"
		"host_loc: the host memory location that will be used  for the benchmark\n",
		"h_case_id: Copy queues for host -> dev \n",
		"h_rev_case_id: Copy queues for dev -> host \n");
  	}

	if(loc < 0 || loc >= CHL_MEMLOCS) error("bw_bench_broadcast_2D: Unsupported src loc = %d\n", loc);

	CHLEnableLinks(loc, CHL_WORKERS);

	int maxDim = std::min(MAX_DIM_TRANS, (int) CHLGetMaxDimSqAsset2D(CHL_WORKERS, elem_size, STEP_TRANS, loc));
	for(int dev_id_idx = 0 ; dev_id_idx < CHL_WORKERS; dev_id_idx++){
		maxDim = std::min(maxDim, (int) CHLGetMaxDimSqAsset2D(1, elem_size, STEP_TRANS, (dev_id_idx)));
	}
	long long ldim = maxDim;
	fprintf(stderr,"\nbw_bench_broadcast_2D: \nSystem = %s\nmaxDim = %d, ldim = %lld\n", 
		TESTBED, maxDim, ldim);
	fprintf(stderr,"-------------------------------------------------------------------------------"
		"-----------------------------------------------------------------------\n");
	int active_unit_num = CHL_WORKERS, active_unit_id_list[CHL_WORKERS];
	double temp_bw, temp_dummy;
	for(int dev_id_idx = 0 ; dev_id_idx < CHL_WORKERS; dev_id_idx++) 
		active_unit_id_list[dev_id_idx] = dev_id_idx;

	double timer = csecond();
	void* loc_buffs[CHL_WORKERS], *worker_buffs[CHL_WORKERS];
	for(int dev_id_idx = 0 ; dev_id_idx < CHL_WORKERS; dev_id_idx++){
		loc_buffs[dev_id_idx] = CHLMalloc(ldim*ldim*elemSize, loc, 1);
		worker_buffs[dev_id_idx] = CHLMalloc(ldim*ldim*elemSize, dev_id_idx, 1);
	}
	timer = csecond() - timer;
	fprintf(stderr, "Allocation buffer size = (%lld x %lld) x %d complete:\t alloc_timer=%lf ms\n", ldim, ldim, elemSize, timer  * 1000);
	fprintf(stderr, "-------------------------------------------------------------------------------"
		"-----------------------------------------------------------------------\n");
	CQueue_p queue_list[CHL_WORKERS];
	Event_timer_p device_timer[CHL_WORKERS];
	for(int dev_id_idx = 0 ; dev_id_idx < CHL_WORKERS; dev_id_idx++){
		//printf("dev_id = %d, dev_id_idx = %d, dev_id_idy = %d, CHL_WORKERS = %d\n", dev_id, dev_id_idx, dev_id_idy, CHL_WORKERS);
		int queue_id = (loc >= CHL_WORKERS || loc < 0)? dev_id_idx : loc;
		queue_list[dev_id_idx] = new CommandQueue(queue_id, COMMUNICATION);
		device_timer[dev_id_idx] = new Event_timer(queue_id);
	}

	fprintf(stderr, "Warming up");
	/// Warmup.
	for (int it = 0; it < 3; it++){
		fprintf(stderr, ".");
		for(int dev_id_idx = 0 ; dev_id_idx < CHL_WORKERS; dev_id_idx++){
			if(loc!=dev_id_idx) queue_list[dev_id_idx]->memcpy2DAsync(worker_buffs[dev_id_idx], ldim,
				loc_buffs[dev_id_idx], ldim, maxDim, maxDim, elemSize, dev_id_idx, loc, 1);
			queue_list[dev_id_idx]->sync_barrier();
		}
	}
	fprintf(stderr, " complete.\n-------------------------------------------------------------------------------"
		"-----------------------------------------------------------------------\n");
	CHLSyncCheckErr();
	
	for (int dim = 256; dim <= maxDim; dim*=2){
		int sample_sz;
		double dev_t[CHL_WORKERS], transfer_t_vals[CHL_WORKERS][MICRO_MAX_ITER] = {0}, 
			transfer_t_sum[CHL_WORKERS] = {0}, transfer_t_mean[CHL_WORKERS] = {0}, 
			error_margin[CHL_WORKERS] = {0}, bench_t = csecond();
		for (sample_sz = 1; sample_sz < MICRO_MAX_ITER + 1; sample_sz++){
			for(int active_unit_idx = 0; active_unit_idx < CHL_WORKERS; active_unit_idx++) if(loc != active_unit_idx){ 
				int queue_id = (loc >= CHL_WORKERS || loc < 0)? active_unit_idx : loc;
				CHLSelectDevice(queue_id);
				device_timer[active_unit_idx]->start_point(queue_list[active_unit_idx]);
				queue_list[active_unit_idx]->memcpy2DAsync(worker_buffs[active_unit_idx], ldim,
					loc_buffs[active_unit_idx], ldim, dim, dim, elemSize, active_unit_idx, loc, 1);
				device_timer[active_unit_idx]->stop_point(queue_list[active_unit_idx]);
			}
			for(int active_unit_idx = 0; active_unit_idx < CHL_WORKERS; active_unit_idx++) if(loc != active_unit_idx)
				queue_list[active_unit_idx]->sync_barrier();
			for(int active_unit_idx = 0; active_unit_idx < CHL_WORKERS; active_unit_idx++) if(loc != active_unit_idx){
				int queue_id = (loc >= CHL_WORKERS || loc < 0)? active_unit_idx : loc;
				CHLSelectDevice(queue_id);
				dev_t[active_unit_idx] = device_timer[active_unit_idx]->sync_get_time()/1000;
			}
			CHLSyncCheckErr();
			int complete_flag = 1;
			for(int active_unit_idx = 0; active_unit_idx < CHL_WORKERS; active_unit_idx++) if(loc != active_unit_idx){
				if(!confidence_interval_5_percent(sample_sz, dev_t[active_unit_idx], transfer_t_vals[active_unit_idx], &transfer_t_sum[active_unit_idx], 
					&transfer_t_mean[active_unit_idx], &error_margin[active_unit_idx])) complete_flag = 0;
			}
			if (sample_sz > MICRO_MIN_ITER && complete_flag) break;
		}
		double broadcast_bw[CHL_WORKERS] = {0}, sum_bw = 0;
		for(int active_unit_idx = 0; active_unit_idx < CHL_WORKERS; active_unit_idx++) 
			if(loc != active_unit_idx) sum_bw += broadcast_bw[active_unit_idx] = Gval_per_s(dim*dim*elemSize, transfer_t_mean[active_unit_idx]);
			else broadcast_bw[active_unit_idx] = -1; 
//#ifdef PDEBUG
		fprintf(stderr, "Ran %d itterations for convergence.\n"
			"-> dim = %d, active_unit_id_list = %s :\n\t h2d_bw(link) = %s Gb/s\n\t h2d_bw(sum) = %.2lf Gb/s\n",
			sample_sz, dim, printlist(active_unit_id_list, CHL_WORKERS), printlist(broadcast_bw, CHL_WORKERS), sum_bw);
		fprintf(stderr,"-------------------------------------------------------------------------------"
			"-----------------------------------------------------------------------\n");
//#endif
	}
  	return 0;
}
