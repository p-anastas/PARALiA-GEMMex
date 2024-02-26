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

void log_results(char* filename, int numDev, int* case_id_list, double case_id_bw_list[][3], int &case_id_lists_len){
	FILE* fp = fopen(filename,"w");
	if (!fp) error("changelink_select_grids: log_results: LogFile failed to open");
	for (int idx = 0; idx < case_id_lists_len; idx++)
		fprintf(fp,"%d, %lf,%lf,%lf\n", case_id_list[idx], case_id_bw_list[idx][0], case_id_bw_list[idx][1], case_id_bw_list[idx][2]);
    fclose(fp);
}

void log_results_bid(char* filename, int numDev, int* case_id_list, int* rev_case_id_list, double case_id_bw_list_bid[][3],  double case_id_bw_list_h2d[][3],  double case_id_bw_list_d2h[][3], int &case_id_lists_len){
	FILE* fp = fopen(filename,"w");
	if (!fp) error("changelink_select_grids: log_results: LogFile failed to open");
	for (int idx = 0; idx < case_id_lists_len; idx++)
		fprintf(fp,"%d,%d, %lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", case_id_list[idx], rev_case_id_list[idx], case_id_bw_list_bid[idx][0], case_id_bw_list_bid[idx][1], case_id_bw_list_bid[idx][2]
		, case_id_bw_list_h2d[idx][0], case_id_bw_list_h2d[idx][1], case_id_bw_list_h2d[idx][2]
		, case_id_bw_list_d2h[idx][0], case_id_bw_list_d2h[idx][1], case_id_bw_list_d2h[idx][2]);
    fclose(fp);
}

void remove_case_id(int case_id, int* case_id_list, double case_id_bw_list[][3], int case_id_lists_len){
	int remove_idx = -1; 
	for(int idx = 0; idx < case_id_lists_len; idx++){
		if (case_id_list[idx] == case_id){
			remove_idx = idx;
			break;
		}
	}
	if(-1 != remove_idx){
		fprintf(stderr, "remove_case_id: removing case_id = %d from index = %d\n", case_id, remove_idx);
		for(int idx = remove_idx; idx < case_id_lists_len -1; idx++){
			case_id_list[idx] = case_id_list[idx+1];
			case_id_bw_list[idx][0] = case_id_bw_list[idx+1][0];
			case_id_bw_list[idx][1] = case_id_bw_list[idx+1][1];
			case_id_bw_list[idx][2] = case_id_bw_list[idx+1][2];
		}
	}
}

int main(const int argc, const char *argv[]) {

	int ctr = 1, samples, dev_id, dev_count, log_result_flag = 0;
	int minDim = MIN_DIM_TRANS, maxDim = 0, step = STEP_TRANS, mem_loc = -1;
	if (minDim < 1) error("Transfer Microbench: Bytes must be > 0");

	switch (argc) {
	case (3):
		mem_loc = atoi(argv[ctr++]);
		log_result_flag = atoi(argv[ctr++]);
		break;
	default:
		error("Incorrect input arguments. Usage: ./correct_run mem_loc log_result_flag(=0 or 1)\n");
  	}

	int benchmarks_exist = 0;
	char *filename = NULL;
	if (log_result_flag) {
		filename = (char *) malloc(1024 * sizeof(char));
		for (int dev_idx = 0; dev_idx < CHL_WORKERS; dev_idx++){
			sprintf(filename, "%s/Benchmark-Results/H2D_chl_microbench_devNum-%d_ver-%s.log", DEPLOYDB, dev_idx+1, VERSION);
			benchmarks_exist += check_benchmark(filename);
			sprintf(filename, "%s/Benchmark-Results/D2H_chl_microbench_devNum-%d_ver-%s.log", DEPLOYDB, dev_idx+1, VERSION);
			benchmarks_exist += check_benchmark(filename);
			sprintf(filename, "%s/Benchmark-Results/HbidD_chl_microbench_devNum-%d_ver-%s.log", DEPLOYDB, dev_idx+1, VERSION);
			benchmarks_exist += check_benchmark(filename);
		}
		if (benchmarks_exist == 3*(CHL_WORKERS)){
			fprintf(stderr, "changelink_select_grids: All microbenchmarks are already complete\n");
			exit(1); 
		}
	}

	for(short dev_id_idx = 0 ; dev_id_idx < CHL_WORKERS; dev_id_idx++){
		//printf("dev_id = %d, dev_id_idx = %d, dev_id_idy = %d, CHL_WORKERS = %d\n", dev_id, dev_id_idx, dev_id_idy, CHL_WORKERS);
		short queue_id = (dev_id_idx);
		h2d_queue_list[dev_id_idx] = new CommandQueue(queue_id, COMMUNICATION);
		d2h_queue_list[dev_id_idx] = new CommandQueue(queue_id, COMMUNICATION);
	}

	// Define the max size of a benchmark kernel to run on this machine.
	maxDim = std::min(MAX_DIM_TRANS, (int) CHLGetMaxDimSqAsset2D(2*(CHL_WORKERS), sizeof(double), STEP_TRANS, -1));

	for(short dev_id_idx = 0 ; dev_id_idx < CHL_WORKERS; dev_id_idx++){
		maxDim = std::min(maxDim, (int)CHLGetMaxDimSqAsset2D(2, sizeof(double), STEP_TRANS, (dev_id_idx)));
	}
	maxDim/=8;
  	long long ldhost = std::min((int) (2*maxDim), (int) CHLGetMaxDimSqAsset2D(2*2*(CHL_WORKERS), sizeof(double), STEP_TRANS, -1)),
		lddev = maxDim;
	fprintf(stderr,"\nchangelink_select_grids: \nSystem = %s\nminDim = %d, maxDim = %d, step = %d(adaptive)\nldhost = %d, lddev = %d\n", 
		TESTBED, minDim, maxDim, step, ldhost, lddev);
				fprintf(stderr,"------------------------------------------------------------------------------------------------------------------------------------------------------\n");

	void* host_buffs[CHL_WORKERS*2];
	void* dev_buffs[CHL_WORKERS*2];
	double timer = csecond();
	for(int dev_id_idx = 0 ; dev_id_idx < CHL_WORKERS; dev_id_idx++)
	{
		host_buffs[2*dev_id_idx + 0] = CHLMalloc(ldhost*ldhost*elemSize, mem_loc, 1);
		host_buffs[2*dev_id_idx + 1] = CHLMalloc(ldhost*ldhost*elemSize, mem_loc, 1);
	}
	timer = csecond() - timer;
	fprintf(stderr, "Allocation host_buffs size = (%d x %d) x %d complete:\t alloc_timer=%lf ms\n", ldhost, ldhost, elemSize, timer  * 1000);

	timer = csecond();
	for(short dev_id_idx = 0 ; dev_id_idx < CHL_WORKERS; dev_id_idx++){
		dev_buffs[2*dev_id_idx + 0] = CHLMalloc(lddev*lddev*elemSize, (dev_id_idx), 1);
		dev_buffs[2*dev_id_idx + 1] = CHLMalloc(lddev*lddev*elemSize, (dev_id_idx), 1);
	}
	timer = csecond() - timer;
	fprintf(stderr, "Allocation dev_buffs size = (%d x %d) x %d complete:\t alloc_timer=%lf ms\n", lddev, lddev, elemSize, timer  * 1000);
				fprintf(stderr,"------------------------------------------------------------------------------------------------------------------------------------------------------\n");

	if(mem_loc == -1) mem_loc = CHL_MEMLOCS - 1;
	fprintf(stderr, "Warming up");
	/// Warmup.
	for (int it = 0; it < 3; it++){
		fprintf(stderr, ".");
		for(short dev_id_idx = 0 ; dev_id_idx < CHL_WORKERS; dev_id_idx++){
		h2d_queue_list[dev_id_idx]->memcpy2DAsync(dev_buffs[2*dev_id_idx], lddev,
			host_buffs[2*dev_id_idx], ldhost,
			maxDim, maxDim, elemSize,
			(dev_id_idx), CHL_WORKERS, 1);
			h2d_queue_list[dev_id_idx]->sync_barrier();
		}
		for(short dev_id_idx = 0 ; dev_id_idx < CHL_WORKERS; dev_id_idx++){
			d2h_queue_list[dev_id_idx]->memcpy2DAsync(host_buffs[2*dev_id_idx + 1], ldhost,
				dev_buffs[2*dev_id_idx + 1], lddev,
				maxDim, maxDim, elemSize,
				CHL_WORKERS, (dev_id_idx), 1);
			d2h_queue_list[dev_id_idx]->sync_barrier();
		}
	}
	fprintf(stderr, " complete.\n");
	CHLSyncCheckErr();
	int dim;
	double temp_bw, temp_dummy;

	int active_unit_num, active_unit_id_list[CHL_WORKERS] = {-42};
	double bw_per_devnum_devwise[CHL_WORKERS][MAX_WORKER_CONFIG][CHL_WORKERS];
	double bw_per_devnum_total[CHL_WORKERS][MAX_WORKER_CONFIG];
	for(int idx = 0; idx < CHL_WORKERS; idx++) 
	for(int idx1 = 0; idx1 < MAX_WORKER_CONFIG; idx1++){
		for(int idx2 = 0; idx2 < CHL_WORKERS; idx2++) 
			bw_per_devnum_devwise[idx][idx1][idx2] = -1;
		bw_per_devnum_total[idx][idx1] = -1;
	}

	for(int wk_idx = 0; wk_idx < CHL_WORKERS; wk_idx++)
	for (int case_idx = 0; case_idx < MAX_WORKER_CONFIG; case_idx++){
			int case_id = CHL_INPUT_QUEUES_CASE_IDS[wk_idx][case_idx];
			translate_binary_to_unit_list(case_id, &active_unit_num, active_unit_id_list);
			if(active_unit_num != wk_idx+1) error("clh_2D_microbench_grids: Something was input from CHL_INPUT_QUEUES_CASE_IDS incorrectly\n");
			dim = (int) maxDim/sqrt(active_unit_num);
			perform_microbenchmark_2D(host_buffs, dev_buffs, mem_loc, (((long long) maxDim*maxDim)*elemSize*active_unit_num), lddev, ldhost, 
				dim, active_unit_id_list, active_unit_num, bw_per_devnum_devwise[wk_idx][case_idx], -1, NULL, 0, NULL, &temp_bw, &temp_dummy, &temp_dummy);
			bw_per_devnum_total[wk_idx][case_idx] = temp_bw;
#ifdef CLDEBUG
			fprintf(stderr,"------------------------------------------------------------------------------------------------------------------------------------------------------\n");
#endif
	}
	
	for(int wk_idx = 0; wk_idx < CHL_WORKERS; wk_idx++){
		fprintf(stderr,"------------------------------------------------------------------------------------------------------------------------------------------------------\n\n");
		for (int case_idx = 0; case_idx < MAX_WORKER_CONFIG; case_idx++){
			translate_binary_to_unit_list(CHL_INPUT_QUEUES_CASE_IDS[wk_idx][case_idx], &active_unit_num, active_unit_id_list);
			fprintf(stderr, "H2D Candidate %d for CHL_WORKERS = %d : %s -> MEMLOC %d (name=%s) BW = %.2lf Gb/s\n", 
			case_idx, wk_idx + 1, printlist<int>(active_unit_id_list, active_unit_num), mem_loc, mem_name(mem_loc), bw_per_devnum_total[wk_idx][case_idx]);
		}
		//if(log_result_flag){
		//	sprintf(filename, "%s/Benchmark-Results/H2D_chl_microbench_devNum-%d_ver-%s.log", DEPLOYDB, wk_idx+1, VERSION);
		//	log_results(filename, wk_idx + 1, CHL_INPUT_QUEUES_CASE_IDS[wk_idx], bw_per_devnum_total[wk_idx], MAX_WORKER_CONFIG);
		//}
		fprintf(stderr,"\n------------------------------------------------------------------------------------------------------------------------------------------------------\n");

	}
	return 1;
	/*
	int rev_active_unit_num, rev_active_unit_id_list[CHL_WORKERS] = {-42}, rev_dim;
	int rev_id_per_devnum[CHL_WORKERS][explored_cases] = {-42}, rev_case_id_ctr[CHL_WORKERS] = {0};
	double rev_bw_per_devnum[CHL_WORKERS][explored_cases][3] = {0};
	int rev_peak_bw_idx[CHL_WORKERS];
	for(int idx = 0; idx < CHL_WORKERS; idx++){
		for(int idx1 = 0; idx1 < explored_cases; idx1++){
			rev_id_per_devnum[idx][idx1] = -1;
			for(int idx2 = 0; idx2 < 3; idx2++) rev_bw_per_devnum[idx][idx1][idx2] = -1;
		}
		rev_peak_bw_idx[idx] = 0;
	}
	for (int buff_idx = 0; buff_idx < host_memloc_num; buff_idx++){
		for (int idx = 0; idx < CHL_WORKERS; idx++) rev_case_id_ctr[idx] = 0;
		for (int case_id = 1; case_id < explored_cases; case_id++){
			translate_binary_to_unit_list(case_id, &rev_active_unit_num, rev_active_unit_id_list);
			int cont_flag = 0; 
			if(rev_case_id_ctr[rev_active_unit_num-1]) 
				for (int idy = 0; idy < blocked_case_id_num; idy++) if(is_subset(blocked_case_id_list[idy], case_id)) cont_flag = 1;
			if(cont_flag) continue;
			rev_dim = (int) maxDim/sqrt(rev_active_unit_num);
			if(rev_active_unit_num != 1 && rev_active_unit_num != 2 && rev_active_unit_num != 4 && rev_active_unit_num != 8) continue;
			if (0 == buff_idx) perform_microbenchmark_2D(host_good_buffs, dev_buffs, CHL_WORKERS, (((long long) maxDim*maxDim)*elemSize*rev_active_unit_num), lddev, ldhost, 
				-1, NULL, 0, rev_dim, rev_active_unit_id_list, rev_active_unit_num, &temp_dummy, &temp_bw, &temp_dummy);
			if (1 == buff_idx) perform_microbenchmark_2D(host_interleaved_buffs, dev_buffs, CHL_MEMLOCS -1, (((long long) maxDim*maxDim)*elemSize*rev_active_unit_num), lddev, ldhost, 
				-1, NULL, 0, rev_dim, rev_active_unit_id_list, rev_active_unit_num, &temp_dummy, &temp_bw, &temp_dummy);
			if (2 == buff_idx) perform_microbenchmark_2D(host_bad_buffs, dev_buffs, CHL_WORKERS, (((long long) maxDim*maxDim)*elemSize*rev_active_unit_num), lddev, ldhost, 
				-1, NULL, 0, rev_dim, rev_active_unit_id_list, rev_active_unit_num, &temp_dummy, &temp_bw, &temp_dummy);
			int del_flag = 0; 
			if(rev_case_id_ctr[rev_active_unit_num-1]) for (int idz = 0; idz < CHL_WORKERS; idz++) if(idz < rev_active_unit_num && rev_case_id_ctr[idz] 
				&& temp_bw <= rev_bw_per_devnum[idz][rev_peak_bw_idx[idz]][buff_idx]/4) del_flag = 1; 
			if(del_flag){
				blocked_case_id_list[blocked_case_id_num++] = case_id;
				if(buff_idx) remove_case_id(case_id, rev_id_per_devnum[rev_active_unit_num-1], rev_bw_per_devnum[rev_active_unit_num-1], explored_cases);
//#ifdef CLDEBUG
				fprintf(stderr, "Blocking (rev_active_unit_id_list[%d] = %s since D2H bw for buff_idx = %d (%lf Gb/s) too low\n",
				rev_active_unit_num, printlist<int>(rev_active_unit_id_list, rev_active_unit_num), buff_idx, temp_bw);
//#endif
			}
			if(!del_flag){
				if(temp_bw > rev_bw_per_devnum[rev_active_unit_num-1][rev_peak_bw_idx[rev_active_unit_num-1]][buff_idx]) 
					rev_peak_bw_idx[rev_active_unit_num-1] = rev_case_id_ctr[rev_active_unit_num-1];
				rev_bw_per_devnum[rev_active_unit_num-1][rev_case_id_ctr[rev_active_unit_num-1]][buff_idx] = temp_bw;
				rev_id_per_devnum[rev_active_unit_num-1][rev_case_id_ctr[rev_active_unit_num-1]++] = case_id;
			}
#ifdef CLDEBUG
			fprintf(stderr,"------------------------------------------------------------------------------------------------------------------------------------------------------\n");
#endif
		}
	}
	for(int idx = 0; idx < CHL_WORKERS; idx++){
		if(!rev_case_id_ctr[idx]) continue;
		bubble_sort_3_vals(rev_id_per_devnum[idx], rev_bw_per_devnum[idx], weights, rev_case_id_ctr[idx]);
		fprintf(stderr,"------------------------------------------------------------------------------------------------------------------------------------------------------\n\n");
		for(int idx1 = 0; idx1 < rev_case_id_ctr[idx]; idx1++){
			translate_binary_to_unit_list(rev_id_per_devnum[idx][idx1], &rev_active_unit_num, rev_active_unit_id_list);
			fprintf(stderr, "D2H Candidate %d for CHL_WORKERS = %d : %s -> NUMA BW = %.2lf Gb/s, Inter BW = %.2lf Gb/s, Bad NUMA BW = %.2lf Gb/s\n", 
			idx1, idx + 1, printlist<int>(rev_active_unit_id_list, rev_active_unit_num), rev_bw_per_devnum[idx][idx1][0], 
				rev_bw_per_devnum[idx][idx1][1], rev_bw_per_devnum[idx][idx1][2]);
		}
		if(log_result_flag){
			sprintf(filename, "%s/Benchmark-Results/D2H_chl_microbench_devNum-%d_ver-%s.log", DEPLOYDB, idx+1, VERSION);
			log_results(filename, idx + 1, rev_id_per_devnum[rev_active_unit_num-1], rev_bw_per_devnum[rev_active_unit_num-1], rev_case_id_ctr[rev_active_unit_num-1]); 
		}
		fprintf(stderr,"\n------------------------------------------------------------------------------------------------------------------------------------------------------\n");
	}
	
	int bid_id_per_devnum[CHL_WORKERS][explored_cases] = {-42}, bid_case_id_ctr[CHL_WORKERS] = {0};
	int bid_rev_id_per_devnum[CHL_WORKERS][explored_cases] = {-42};
	double bid_simu_bw_per_devnum[CHL_WORKERS][explored_cases][3] = {0};
	double bid_bw_per_devnum[CHL_WORKERS][explored_cases][3] = {0};
	double bid_rev_bw_per_devnum[CHL_WORKERS][explored_cases][3] = {0};
	int bid_simu_peak_bw_idx[CHL_WORKERS];
	double rev_bw, bid_bw;
	for(int idx = 0; idx < CHL_WORKERS; idx++){
		for(int idx1 = 0; idx1 < explored_cases; idx1++){
			bid_id_per_devnum[idx][idx1] = -1;
			bid_rev_id_per_devnum[idx][idx1] = -1;
			for(int idx2 = 0; idx2 < 3; idx2++){
				bid_simu_bw_per_devnum[idx][idx1][idx2] = bid_bw_per_devnum[idx][idx1][idx2] = bid_rev_bw_per_devnum[idx][idx1][idx2] = -1;
			}
		}
		bid_simu_peak_bw_idx[idx] = 0;
	}
	for (int buff_idx = 0; buff_idx < buff_types; buff_idx++){
		for (int idx = 0; idx < CHL_WORKERS; idx++) bid_case_id_ctr[idx] = 0;
		for (int case_id = 1; case_id < explored_cases; case_id++){
			translate_binary_to_unit_list(case_id, &active_unit_num, active_unit_id_list);
			if(!is_in_list(case_id, id_per_devnum[active_unit_num-1], case_id_ctr[active_unit_num-1]) ||
			!is_in_list(case_id, rev_id_per_devnum[active_unit_num-1], rev_case_id_ctr[active_unit_num-1])) continue;
			if(active_unit_num != 1 && active_unit_num != 2 && active_unit_num != 4 && active_unit_num != 8) continue;
			dim = (int) maxDim/sqrt(active_unit_num);
			int test_best_num = MAX_WORKER_CONFIG;
			int rev_id_list[test_best_num + 1];
			rev_id_list[0] = case_id; 
			if(active_unit_num/2-1 > 0) for (int rev_case_idx = 0; rev_case_idx < test_best_num; rev_case_idx++) 
				rev_id_list[rev_case_idx+1] = rev_id_per_devnum[active_unit_num/2-1][rev_case_idx];
			else test_best_num = 0;

			for (int rev_case_idx = 0; rev_case_idx < test_best_num + 1; rev_case_idx++){
				int rev_case_id = rev_id_list[rev_case_idx];
				translate_binary_to_unit_list(rev_case_id, &rev_active_unit_num, rev_active_unit_id_list);
				rev_dim = (int) maxDim/sqrt(rev_active_unit_num);
				if (0 == buff_idx) perform_microbenchmark_2D(host_good_buffs, dev_buffs, CHL_WORKERS, (((long long) maxDim*maxDim)*elemSize*std::max(active_unit_num, rev_active_unit_num)), lddev, ldhost, 
					dim, active_unit_id_list, active_unit_num, rev_dim, rev_active_unit_id_list, rev_active_unit_num, &temp_bw, &rev_bw, &bid_bw);
				if (1 == buff_idx) perform_microbenchmark_2D(host_interleaved_buffs, dev_buffs, CHL_MEMLOCS - 1, (((long long) maxDim*maxDim)*elemSize*std::max(active_unit_num, rev_active_unit_num)), lddev, ldhost, 
					dim, active_unit_id_list, active_unit_num, rev_dim, rev_active_unit_id_list, rev_active_unit_num, &temp_bw, &rev_bw, &bid_bw);
				if (2 == buff_idx) perform_microbenchmark_2D(host_bad_buffs, dev_buffs, CHL_WORKERS, (((long long) maxDim*maxDim)*elemSize*std::max(active_unit_num, rev_active_unit_num)), lddev, ldhost, 
					dim, active_unit_id_list, active_unit_num, rev_dim, rev_active_unit_id_list, rev_active_unit_num, &temp_bw, &rev_bw, &bid_bw);
				if(bid_bw > bid_simu_bw_per_devnum[active_unit_num-1][bid_simu_peak_bw_idx[active_unit_num-1]][buff_idx]) 
					bid_simu_peak_bw_idx[active_unit_num-1] = bid_case_id_ctr[active_unit_num-1];
				bid_simu_bw_per_devnum[active_unit_num-1][bid_case_id_ctr[active_unit_num-1]][buff_idx] = bid_bw;
				bid_bw_per_devnum[active_unit_num-1][bid_case_id_ctr[active_unit_num-1]][buff_idx] = temp_bw;
				bid_rev_bw_per_devnum[active_unit_num-1][bid_case_id_ctr[active_unit_num-1]][buff_idx] = rev_bw;
				bid_id_per_devnum[active_unit_num-1][bid_case_id_ctr[active_unit_num-1]] = case_id;
				bid_rev_id_per_devnum[active_unit_num-1][bid_case_id_ctr[active_unit_num-1]++] = rev_case_id;
	#ifdef CLDEBUG
				fprintf(stderr,"------------------------------------------------------------------------------------------------------------------------------------------------------\n");
	#endif	
			}
		}
	}
	
	for(int idx = 0; idx < CHL_WORKERS; idx++){
		if(!bid_case_id_ctr[idx]) continue;
		bubble_sort_2_idx_3x3_vals(bid_id_per_devnum[idx], bid_rev_id_per_devnum[idx], 
			bid_simu_bw_per_devnum[idx], bid_bw_per_devnum[idx], bid_rev_bw_per_devnum[idx], weights, bid_case_id_ctr[idx]);
		fprintf(stderr,"------------------------------------------------------------------------------------------------------------------------------------------------------\n\n");
		for(int idx1 = 0; idx1 < bid_case_id_ctr[idx]; idx1++){
			translate_binary_to_unit_list(bid_id_per_devnum[idx][idx1], &active_unit_num, active_unit_id_list);
			translate_binary_to_unit_list(bid_rev_id_per_devnum[idx][idx1], &rev_active_unit_num, rev_active_unit_id_list);
			fprintf(stderr, "H<->D (HbidD) candidate %d for CHL_WORKERS = %d : H2D : %s, D2H: %s\n-> Bid NUMA BW = %.2lf Gb/s, Bid inter BW = %.2lf Gb/s, Bid bad NUMA BW = %.2lf Gb/s\n"
			"-> H2D NUMA BW = %.2lf Gb/s, H2D inter BW = %.2lf Gb/s, H2D bad NUMA BW = %.2lf Gb/s\n-> D2H NUMA BW = %.2lf Gb/s, D2H inter BW = %.2lf Gb/s, D2H bad NUMA BW = %.2lf Gb/s\n\n", 
			idx1, idx + 1, printlist<int>(active_unit_id_list, active_unit_num), printlist<int>(rev_active_unit_id_list, rev_active_unit_num),
			bid_simu_bw_per_devnum[idx][idx1][0], bid_simu_bw_per_devnum[idx][idx1][1], 
			bid_simu_bw_per_devnum[idx][idx1][2],
			bid_bw_per_devnum[idx][idx1][0], bid_bw_per_devnum[idx][idx1][1], bid_bw_per_devnum[idx][idx1][2],
			bid_rev_bw_per_devnum[idx][idx1][0], bid_rev_bw_per_devnum[idx][idx1][1], bid_rev_bw_per_devnum[idx][idx1][2]);

		}
		if(log_result_flag){
			sprintf(filename, "%s/Benchmark-Results/HbidD_chl_microbench_devNum-%d_ver-%s.log", DEPLOYDB, idx+1, VERSION);
			log_results_bid(filename, idx + 1, bid_id_per_devnum[active_unit_num-1], bid_rev_id_per_devnum[active_unit_num-1], 
			bid_simu_bw_per_devnum[active_unit_num-1], bid_bw_per_devnum[active_unit_num-1], bid_rev_bw_per_devnum[active_unit_num-1], bid_case_id_ctr[active_unit_num-1]); 
		}
		fprintf(stderr,"\n------------------------------------------------------------------------------------------------------------------------------------------------------\n");
	}*/
	CHLSyncCheckErr();
	fprintf(stderr,"------------------------------------------------------------------------------------------------------------------------------------------------------\n");
	timer = csecond();
	for(int dev_id_idx = 0 ; dev_id_idx < CHL_WORKERS; dev_id_idx++){
  		CHLFree(host_buffs[2*dev_id_idx + 0], ldhost*ldhost*elemSize, -1);
  		CHLFree(host_buffs[2*dev_id_idx + 1], ldhost*ldhost*elemSize, -1);
	}
	timer = csecond() - timer;
	fprintf(stderr, "Free host_buffs size = (%d x %d) x %d complete:\t alloc_timer=%lf ms\n", ldhost, ldhost, elemSize, timer  * 1000);

	timer = csecond();
	for(short dev_id_idx = 0 ; dev_id_idx < CHL_WORKERS; dev_id_idx++){
  		CHLFree(dev_buffs[2*dev_id_idx + 0], lddev*lddev*elemSize, (dev_id_idx));
  		CHLFree(dev_buffs[2*dev_id_idx + 1], lddev*lddev*elemSize, (dev_id_idx));
	}
	timer = csecond() - timer;
	fprintf(stderr, "Free dev_buffs size = (%d x %d) x %d complete:\t alloc_timer=%lf ms\n", lddev, lddev, elemSize, timer  * 1000);
	fprintf(stderr,"------------------------------------------------------------------------------------------------------------------------------------------------------\n");
  	return 0;
}
