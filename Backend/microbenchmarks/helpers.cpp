///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief Some convinient C/C++ utilities for CHL.
///

#include "smart_wrappers.hpp"
#include "grid_amalgamation.hpp"
#include "microbenchmarks.hpp"

#include <numa.h>

int use_square_Tiles = 0; 
int elemSize = sizeof(double);
CQueue_p h2d_queue_list[32] = {NULL}, d2h_queue_list[32] = {NULL};

double latbw_linear_regression(double* y, double* x, int samples, double* ans_b, double* ans_a, int active_unit_num){
 	/*	
	double sumX=0, sumX2=0, sumY=0, sumXY=0;
 	for(int i=0; i<samples; i++){
		sumX = sumX + x[i];
		sumX2 = sumX2 + x[i]*x[i];
		sumY = sumY + y[i];
		sumXY = sumXY + x[i]*y[i];
	}
 	//Calculating a and b 
	//*ans_a = (samples*sumXY-sumX*sumY)/(samples*sumX2-sumX*sumX); /// TODO: These result in negative latencies...
	//*ans_b = (sumY*sumX2-sumX*sumXY)/(samples*sumX2-sumX*sumX); //(sumY - (*ans_b)*sumX)/samples; /// TODO: And NNLS seems like an overkill...
	*/
	*ans_a = y[samples-1]/x[samples-1]; // Assuming the last sample has near-0 latency
	double sumb = 0; 
	int active_samples = samples; 
	for(int i=0; i<samples; i++){
		double rem_lat = y[i] - (*ans_a)*x[i];
		if(rem_lat >= 0.0) sumb += rem_lat;
		else active_samples--;
	}
	*ans_a = (y[samples-1]- sumb/active_samples)/x[samples-1];
	*ans_b = sumb/active_samples;
	//*ans_b = abs(sumY - (*ans_a)*sumX)/samples;
	double t_lat = *ans_b, t_b = *ans_a;
	fprintf(stderr,"Latency estimated by linear regresion :\t t_lat = %e\n", t_lat);
	fprintf(stderr,"BW-inv estimated by linear regresion:\t t_b = %e ( %lf Gb/s)\n", t_b, (1/t_b)/1e9);
	long byte_breakpoint = (int) (t_lat/t_b); 
	if (!use_square_Tiles){
		fprintf(stderr,"Latency-BW estimated breakpoint:\t b_eq = %ld bytes\n", byte_breakpoint/active_unit_num);
		fprintf(stderr,"75%% BW estimated breakpoint:\t\t b_75 = %ld bytes\n", byte_breakpoint/active_unit_num*2);
		fprintf(stderr,"95%% BW estimated breakpoint:\t\t b_95 = %ld bytes\n\n", byte_breakpoint/active_unit_num*20);
	}
	else{
		fprintf(stderr,"Latency-BW estimated breakpoint:\t b_eq = 8-bytes X (%d X %d)\n", (int) sqrt(byte_breakpoint/active_unit_num/8), (int) sqrt(byte_breakpoint/active_unit_num/8));
		fprintf(stderr,"75%% BW estimated breakpoint:\t\t b_75 = 8-bytes X (%d X %d)\n", (int) sqrt(byte_breakpoint/active_unit_num/4), (int) sqrt(byte_breakpoint/active_unit_num/4));
		fprintf(stderr,"95%% BW estimated breakpoint:\t\t b_95 = 8-bytes X (%d X %d)\n\n", (int) sqrt(byte_breakpoint/active_unit_num*20/8), (int) sqrt(byte_breakpoint/active_unit_num*20/8));
	}
	double APE = 0; 
	long long byte_lbw_bp = 0, byte_75pbw_bp = 0, byte_95pbw_bp = 0;
	for(int i=0; i<samples; i++){
		double pred_y =  t_lat + x[i]*t_b; 
#ifdef PDEBUG
		fprintf(stderr, "Point (x = %d):\t y=%lf ms (%lf Gb/s), pred_y = %lf ms (%lf Gb/s) PE = %.1lf\n", 
		int(x[i]), y[i] * 1000, Gval_per_s(int(x[i]), y[i]), pred_y * 1000, Gval_per_s(int(x[i]), pred_y), (pred_y - y[i])/y[i]*100);
#endif
		APE += abs((pred_y - y[i])/y[i]*100);
		if (!byte_lbw_bp && y[i]/x[i] <= (*ans_a)/0.5) byte_lbw_bp = x[i];
		if (!byte_75pbw_bp && y[i]/x[i] <= (*ans_a)/0.75) byte_75pbw_bp = x[i];
		if (!byte_95pbw_bp && y[i]/x[i] <= (*ans_a)/0.95) byte_95pbw_bp = x[i];

	}
	if (!use_square_Tiles){
		fprintf(stderr,"Latency-BW empirical breakpoint:\t b_eq = %lld bytes\n", byte_lbw_bp/active_unit_num);
		fprintf(stderr,"75%% BW empirical breakpoint:\t\t b_75 = %lld bytes\n", byte_75pbw_bp/active_unit_num);
		fprintf(stderr,"95%% BW empirical breakpoint:\t\t b_95 = %lld bytes\n\n", byte_95pbw_bp/active_unit_num);
	}
	else{
		fprintf(stderr,"Latency-BW empirical breakpoint:\t b_eq = 8-bytes X (%d X %d)\n", (int) sqrt(byte_lbw_bp/active_unit_num/8), (int) sqrt(byte_lbw_bp/active_unit_num/8));
		fprintf(stderr,"75%% BW empirical breakpoint:\t\t b_75 = 8-bytes X (%d X %d)\n", (int) sqrt(byte_75pbw_bp/active_unit_num/8), (int) sqrt(byte_75pbw_bp/active_unit_num/8));
		fprintf(stderr,"95%% BW empirical breakpoint:\t\t b_95 = 8-bytes X (%d X %d)\n\n", (int) sqrt(byte_95pbw_bp/active_unit_num/8), (int) sqrt(byte_95pbw_bp/active_unit_num/8));
	}

	return APE/samples; 
}
int confidence_interval_5_percent(long int sample_sz, double cpu_timer, double* transfer_t_vals, double* transfer_t_sum_ptr, double* transfer_t_mean_ptr, double* error_margin_ptr){
	transfer_t_vals[sample_sz-1] = cpu_timer;
	(*transfer_t_sum_ptr) += transfer_t_vals[sample_sz-1];
	(*transfer_t_mean_ptr) = (*transfer_t_sum_ptr)/sample_sz;
	if (sample_sz < 2) return 0;
	double std_dev = 0; 
	for (int i = 0; i < sample_sz; i++) std_dev += pow(transfer_t_vals[i] - (*transfer_t_mean_ptr), 2);
	std_dev /= sample_sz;
		std_dev = sqrt(std_dev);
	boost::math::students_t dist(sample_sz - 1);
	double Td = boost::math::quantile(boost::math::complement(dist, alphaCI / 2));
	(*error_margin_ptr) = Td*std_dev/sqrt(sample_sz);
#ifdef DPDEBUG
	fprintf(stderr, "\tItter %ld:\t mean=%lf, std_dev = %lf, Error margin =%lf\n", sample_sz, (*transfer_t_mean_ptr) , std_dev, (*error_margin_ptr));
#endif
	if ((*error_margin_ptr)/(*transfer_t_mean_ptr)  * 100 <= 0.5) return 1;
	else return 0;
}

/// @brief Initialize and benchmark the affinity of the system GPUs with all numa nodes.
/// @param microbench_bytes (input): The ammount of bytes that will be used for numa microbenchmarks 
/// @param numa_high_affinity (output): The highest-affinity score node for each GPU 
/// @param numa_low_affinity (output): The lowest-affinity score node for each GPU 
/// @param nuMap_weights_h2d (output): (A pointer to) The affinity score map for each GPU<->numa_node for h2d trasfers
/// @param nuMap_weights_d2h (output): (A pointer to) The affinity score map for each GPU<->numa_node for d2h trasfers
/// @param nuMap_weights_bid (output): (A pointer to) The affinity score map for each GPU<->numa_node for bidirectional trasfers
void find_numa_connectivity(long long microbench_bytes, int* numa_high_affinity, int* numa_low_affinity, double *** nuMap_weights_h2d, double *** nuMap_weights_d2h, double *** nuMap_weights_bid){
	int numa_nodes = numa_max_node() + 1;
	fprintf(stderr,"------------------------------------------------------------------------------------------------------------------------------------------------------\n");
	fprintf(stderr, "find_numa_connectivity() : GPUs = %d, Numa nodes = %d\n", CHL_WORKERS, numa_nodes);
	// TODO: Not clear how to do this automatically...manually for now. Hints for obtaining manual:
	// TODO: 1) nvidia-smi topo --matrix -> [NUMA affinity]
	// TODO: 2) If NUMA affinity from (1) is not defined combine: 
	// TODO: 	a) nvidia-smi topo --matrix -> [GPU<->CPU affinity]
	// TODO:	b) numactl -H -> [NUMA<->CPU affinity]
	for (int idx = 0; idx < CHL_WORKERS; idx++){
		numa_high_affinity[idx] = idx/((CHL_WORKERS)/numa_nodes)/2*2 + 1;
		numa_low_affinity[idx] = (CHL_WORKERS - 1) - numa_high_affinity[idx];
	}

	/*int dim = sqrt(microbench_bytes/elemSize); 
	*nuMap_weights_h2d = (double **) malloc((CHL_WORKERS)*sizeof(double*));
	*nuMap_weights_d2h = (double **) malloc((CHL_WORKERS)*sizeof(double*));
	*nuMap_weights_bid = (double **) malloc((CHL_WORKERS)*sizeof(double*));
	for (int idx = 0; idx < CHL_WORKERS; idx++){
		(*nuMap_weights_h2d)[idx] = (double *) malloc(numa_nodes*sizeof(double));
		(*nuMap_weights_d2h)[idx] = (double *) malloc(numa_nodes*sizeof(double));
		(*nuMap_weights_bid)[idx] = (double *) malloc(numa_nodes*sizeof(double));
	}
	void* loc_buffs[numa_nodes][2];
	void* worker_buffs[CHL_WORKERS][2];
	for(int numa_id_idx = 0 ; numa_id_idx < numa_nodes; numa_id_idx++){
		loc_buffs[numa_id_idx][0] = numa_bind_pin_malloc(microbench_bytes, numa_id_idx, 1);
		loc_buffs[numa_id_idx][1] = numa_bind_pin_malloc(microbench_bytes, numa_id_idx, 1);
	}
	for(int dev_id_idx = 0 ; dev_id_idx < CHL_WORKERS; dev_id_idx++){
		worker_buffs[dev_id_idx][0] = CHLMalloc(microbench_bytes, deidxize(dev_id_idx), 1);
		worker_buffs[dev_id_idx][1] = CHLMalloc(microbench_bytes, deidxize(dev_id_idx), 1);
	}

	for(int dev_id_idx = 0 ; dev_id_idx < CHL_WORKERS; dev_id_idx++){
		int loc_dest = dev_id_idx;
		CHLSelectDevice(loc_dest);
		for(int numa_id_idx = 0 ; numa_id_idx < numa_nodes; numa_id_idx++){
			double cpu_timer, transfer_t_vals[MICRO_MAX_ITER] = {0}, transfer_t_sum = 0, transfer_t_mean = 0, error_margin = 0;
			double error_margin_bid;		
			int sample_sz = 0;
			for (sample_sz = 1; sample_sz < MICRO_MAX_ITER + 1; sample_sz++){
				cpu_timer = csecond();
				for (int reps = 0; reps < INTERCONNECT_LOAD_TRANSFERS; reps++)
				if (use_square_Tiles) CHLMemcpy2DAsync(worker_buffs[idxize(loc_dest)][0], dim,
					loc_buffs[numa_id_idx][0], dim, dim, dim, elemSize,
					loc_dest, -1, h2d_queue_list[idxize(loc_dest)]);
				else CHLMemcpyAsync(worker_buffs[idxize(loc_dest)][0], loc_buffs[numa_id_idx][0],
					microbench_bytes, loc_dest, -1, h2d_queue_list[idxize(loc_dest)]);
				h2d_queue_list[idxize(loc_dest)]->sync_barrier();
				//CHLSyncCheckErr();
				cpu_timer  = csecond() - cpu_timer;
				int complete_flag = confidence_interval_5_percent(sample_sz, cpu_timer, transfer_t_vals, &transfer_t_sum, &transfer_t_mean, &error_margin);
				if (sample_sz > MICRO_MIN_ITER && complete_flag) break;
			}
			transfer_t_mean/=INTERCONNECT_LOAD_TRANSFERS;
			(*nuMap_weights_h2d)[dev_id_idx][numa_id_idx] = Gval_per_s(dim*dim*elemSize, transfer_t_mean);
			fprintf(stderr,"------------------------------------------------------------------------------------------------------------------------------------------------------\n");
			fprintf(stderr, "NUMA microbenchmark ( NuNode -> dev [%d -> %d], dim = %d, total_bytes = %ld) complete: error Margin (percentage of mean) = %lf %%, Itter = %ld"
				"\n\tH2D mean_exec_t=%lf ms  ( %lf Gb/s)\n", numa_id_idx, loc_dest, dim, dim*dim*elemSize,
				error_margin/(transfer_t_mean*(INTERCONNECT_LOAD_TRANSFERS)) * 100, sample_sz, transfer_t_mean  * 1000, Gval_per_s(dim*dim*elemSize, transfer_t_mean));
			sample_sz = 0;
			transfer_t_sum = transfer_t_mean = error_margin = 0;
			for (sample_sz = 1; sample_sz < MICRO_MAX_ITER + 1; sample_sz++){
				cpu_timer = csecond();
				for (int reps = 0; reps < INTERCONNECT_LOAD_TRANSFERS; reps++)
				if (use_square_Tiles) CHLMemcpy2DAsync(loc_buffs[numa_id_idx][1], dim, worker_buffs[idxize(loc_dest)][1], dim,
					dim, dim, elemSize,
					-1, loc_dest, d2h_queue_list[idxize(loc_dest)]);
				else CHLMemcpyAsync(loc_buffs[numa_id_idx][1], worker_buffs[idxize(loc_dest)][1],
					microbench_bytes, -1, loc_dest, d2h_queue_list[idxize(loc_dest)]);
				d2h_queue_list[idxize(loc_dest)]->sync_barrier();
				//CHLSyncCheckErr();
				cpu_timer  = csecond() - cpu_timer;
				int complete_flag = confidence_interval_5_percent(sample_sz, cpu_timer, transfer_t_vals, &transfer_t_sum, &transfer_t_mean, &error_margin);
				if (sample_sz > MICRO_MIN_ITER && complete_flag) break;
			}
			transfer_t_mean/=INTERCONNECT_LOAD_TRANSFERS;
			(*nuMap_weights_d2h)[dev_id_idx][numa_id_idx] = Gval_per_s(dim*dim*elemSize, transfer_t_mean);
			fprintf(stderr,"------------------------------------------------------------------------------------------------------------------------------------------------------\n");
			fprintf(stderr, "NUMA microbenchmark ( dev -> NuNode [%d -> %d], dim = %d, total_bytes = %ld) complete: error Margin (percentage of mean) = %lf %%, Itter = %ld"
				"\n\tD2H mean_exec_t=%lf ms  ( %lf Gb/s)\n", loc_dest, numa_id_idx,  dim, dim*dim*elemSize,
				error_margin/(transfer_t_mean*(INTERCONNECT_LOAD_TRANSFERS)) * 100, sample_sz, transfer_t_mean  * 1000, Gval_per_s(dim*dim*elemSize, transfer_t_mean));
			sample_sz = 0;
			transfer_t_sum = transfer_t_mean = error_margin = 0;
			for (sample_sz = 1; sample_sz < MICRO_MAX_ITER + 1; sample_sz++){
				cpu_timer = csecond();
				for (int reps = 0; reps < INTERCONNECT_LOAD_TRANSFERS; reps++)
				if (use_square_Tiles) CHLMemcpy2DAsync(worker_buffs[idxize(loc_dest)][0], dim,
					loc_buffs[numa_id_idx][0], dim, dim, dim, elemSize,
					loc_dest, -1, h2d_queue_list[idxize(loc_dest)]);
				else CHLMemcpyAsync(worker_buffs[idxize(loc_dest)][0], loc_buffs[numa_id_idx][0],
					microbench_bytes, loc_dest, -1, h2d_queue_list[idxize(loc_dest)]);
				for (int reps = 0; reps < INTERCONNECT_LOAD_TRANSFERS; reps++)
				if (use_square_Tiles) CHLMemcpy2DAsync(loc_buffs[numa_id_idx][1], dim, worker_buffs[idxize(loc_dest)][1], dim,
					dim, dim, elemSize,
					-1, loc_dest, d2h_queue_list[idxize(loc_dest)]);
				else CHLMemcpyAsync(loc_buffs[numa_id_idx][1], worker_buffs[idxize(loc_dest)][1],
					microbench_bytes, -1, loc_dest, d2h_queue_list[idxize(loc_dest)]);
				h2d_queue_list[idxize(loc_dest)]->sync_barrier();
				d2h_queue_list[idxize(loc_dest)]->sync_barrier();
				//CHLSyncCheckErr();
				cpu_timer  = csecond() - cpu_timer;
				int complete_flag = confidence_interval_5_percent(sample_sz, cpu_timer, transfer_t_vals, &transfer_t_sum, &transfer_t_mean, &error_margin);
				if (sample_sz > MICRO_MIN_ITER && complete_flag) break;
			}
			transfer_t_mean/=INTERCONNECT_LOAD_TRANSFERS;
			(*nuMap_weights_bid)[dev_id_idx][numa_id_idx] = Gval_per_s(2*dim*dim*elemSize, transfer_t_mean);
			fprintf(stderr,"------------------------------------------------------------------------------------------------------------------------------------------------------\n");
			fprintf(stderr, "NUMA microbenchmark ( dev <-> NuNode [%d <-> %d], dim = %d, total_bytes = %ld) complete: error Margin (percentage of mean) = %lf %%, Itter = %ld"
				"\n\tD2H mean_exec_t=%lf ms  ( %lf Gb/s)\n", loc_dest, numa_id_idx,  dim, 2*dim*dim*elemSize,
				error_margin/(transfer_t_mean*(INTERCONNECT_LOAD_TRANSFERS)) * 100, sample_sz, transfer_t_mean  * 1000, Gval_per_s(2*dim*dim*elemSize, transfer_t_mean));
		}
	}*/
	for (int idx = 0; idx < CHL_WORKERS; idx++) fprintf(stderr, "->\tDevice %d: High affinity NuNo = %d, Low affinity NuNo = %d\n", idx, numa_high_affinity[idx], numa_low_affinity[idx]);
	fprintf(stderr,"------------------------------------------------------------------------------------------------------------------------------------------------------\n");
}


double perform_microbenchmark_2D(void** loc_buffs, void** worker_buffs, int loc, long long total_bytes_per_side, int lddev, int ldhost, 
	int dim, int* unit_ids, int unit_num, double* worker_wise_bws, int rev_dim, int* rev_unit_ids, int rev_unit_num, double* worker_wise_rev_bws, double* h2d_bw, double* d2h_bw, double* bid_bw){
	Event_timer_p device_timer[CHL_WORKERS], rev_device_timer[CHL_WORKERS];
	for(int idx = 0; idx < CHL_WORKERS; idx++){
		device_timer[idx] = new Event_timer(idx);
		rev_device_timer[idx] = new Event_timer(idx);
	}
#ifdef PDEBUG
	fprintf(stderr,"------------------------------------------------------------------------------------------------------------------------------------------------------\n");
	fprintf(stderr, "perform_microbenchmark_2D: loc = %d, total_bytes_per_side = %lld MB, lddev = %d, ldhost =%d -> ", loc, (long long) (total_bytes_per_side/1e6), lddev, ldhost);
#endif

	double dev_t[CHL_WORKERS], rev_dev_t[CHL_WORKERS];
	double transfer_t_vals[2][CHL_WORKERS][MICRO_MAX_ITER] = {0}, transfer_t_sum[2][CHL_WORKERS] = {0}, transfer_t_mean[2][CHL_WORKERS] = {0}, error_margin[2][CHL_WORKERS] = {0}, bench_t = csecond();
	long long trasfer_mult = (long long) (1.0*total_bytes_per_side/(dim*dim*elemSize*unit_num) + 1), 
		trasfer_mult_rev = (long long) (1.0*total_bytes_per_side/(rev_dim*rev_dim*elemSize*rev_unit_num) + 1);
	int sample_sz = 0;
	for (sample_sz = 1; sample_sz < MICRO_MAX_ITER + 1; sample_sz++){
		for(int active_unit_idx = 0; active_unit_idx < unit_num; active_unit_idx++){
			int loc_dest = unit_ids[active_unit_idx];
#ifdef PDEBUG
			//fprintf(stderr, "perform_microbenchmark_2D: Running microbench itter %d\n", sample_sz-1);
			int memloc = CHLGetPtrLoc(loc_buffs[2*loc_dest]);
			if(memloc != loc) warning("perform_microbenchmark_2D: Suspicious input buffer %p for loc = %d (%s), pointed to %d (%s) instead\n", 
				loc_buffs[2*loc_dest], loc, mem_name(loc), memloc, mem_name(memloc));
#endif
			CHLSelectDevice(loc_dest);
			device_timer[loc_dest]->start_point(h2d_queue_list[(loc_dest)]);
			if(loc != loc_dest) for (long long reps = 0; reps < trasfer_mult*INTERCONNECT_LOAD_TRANSFERS; reps++)
			h2d_queue_list[loc_dest]->memcpy2DAsync(worker_buffs[2*loc_dest], lddev,
				loc_buffs[2*loc_dest], ldhost, dim, dim, elemSize,
				loc_dest, loc, 1);
			device_timer[loc_dest]->stop_point(h2d_queue_list[(loc_dest)]);
		}
		for(int active_unit_idx = 0; active_unit_idx < rev_unit_num; active_unit_idx++){
			int loc_rev_src = rev_unit_ids[active_unit_idx];
			CHLSelectDevice(loc_rev_src);
			rev_device_timer[loc_rev_src]->start_point(d2h_queue_list[(loc_rev_src)]);
			if(loc != loc_rev_src) for (long long reps = 0; reps < trasfer_mult_rev*INTERCONNECT_LOAD_TRANSFERS; reps++)
			d2h_queue_list[(loc_rev_src)]->memcpy2DAsync(loc_buffs[2*loc_rev_src + 1], ldhost, 
				worker_buffs[2*loc_rev_src + 1], lddev, rev_dim, rev_dim, elemSize,
				loc, loc_rev_src, 1);
			rev_device_timer[loc_rev_src]->stop_point(d2h_queue_list[(loc_rev_src)]);
		}
		for(int active_unit_idx = 0; active_unit_idx < CHL_WORKERS; active_unit_idx++){
			h2d_queue_list[active_unit_idx]->sync_barrier();
			d2h_queue_list[active_unit_idx]->sync_barrier();
		}
		for(int active_unit_idx = 0; active_unit_idx < unit_num; active_unit_idx++){
			int loc_dest = unit_ids[active_unit_idx];
			CHLSelectDevice(loc_dest);
			dev_t[(loc_dest)] = device_timer[(loc_dest)]->sync_get_time()/1000/(trasfer_mult*INTERCONNECT_LOAD_TRANSFERS);
		}
		for(int active_unit_idx = 0; active_unit_idx < rev_unit_num; active_unit_idx++){
			int loc_rev_src = rev_unit_ids[active_unit_idx];
			CHLSelectDevice(loc_rev_src);
			rev_dev_t[(loc_rev_src)] = rev_device_timer[(loc_rev_src)]->sync_get_time()/1000/(trasfer_mult_rev*INTERCONNECT_LOAD_TRANSFERS);
		}
		CHLSyncCheckErr();
		int complete_flag = 1;
		for(int active_unit_idx = 0; active_unit_idx < unit_num; active_unit_idx++){
			int loc_dest = unit_ids[active_unit_idx];
			if(!confidence_interval_5_percent(sample_sz, dev_t[(loc_dest)], transfer_t_vals[0][(loc_dest)], &transfer_t_sum[0][(loc_dest)], 
				&transfer_t_mean[0][(loc_dest)], &error_margin[0][(loc_dest)])) complete_flag = 0;
		}
		for(int active_unit_idx = 0; active_unit_idx < rev_unit_num; active_unit_idx++){
			int loc_rev_src = rev_unit_ids[active_unit_idx];
			if(!confidence_interval_5_percent(sample_sz, rev_dev_t[(loc_rev_src)], transfer_t_vals[1][(loc_rev_src)], &transfer_t_sum[1][(loc_rev_src)], 
				&transfer_t_mean[1][(loc_rev_src)], &error_margin[1][(loc_rev_src)])) complete_flag = 0;
		}
		if (sample_sz > MICRO_MIN_ITER && complete_flag) break;
	}
	(*h2d_bw) = (*d2h_bw) = -1;
	for(int active_unit_idx = 0; active_unit_idx < unit_num; active_unit_idx++){
		int loc_dest = unit_ids[active_unit_idx];
		(*h2d_bw) = fmax((*h2d_bw), transfer_t_mean[0][loc_dest]);
		if(worker_wise_bws && transfer_t_mean[0][loc_dest]) worker_wise_bws[active_unit_idx] = Gval_per_s(dim*dim*elemSize, transfer_t_mean[0][loc_dest]);
		else if(worker_wise_bws) worker_wise_bws[active_unit_idx] = -1;
	}
	for(int active_unit_idx = 0; active_unit_idx < rev_unit_num; active_unit_idx++){
		int loc_rev_src = rev_unit_ids[active_unit_idx];
		(*d2h_bw) = fmax((*d2h_bw), transfer_t_mean[1][loc_rev_src]);
		if(worker_wise_rev_bws && transfer_t_mean[1][loc_rev_src]) worker_wise_rev_bws[active_unit_idx] = Gval_per_s(dim*dim*elemSize, transfer_t_mean[1][loc_rev_src]);
		else if(worker_wise_rev_bws) worker_wise_rev_bws[active_unit_idx] = -1;


	}
	if (unit_num && rev_unit_num) (*bid_bw) = Gval_per_s(((long long)dim*dim*elemSize*unit_num) + ((long long)rev_dim*rev_dim*elemSize*rev_unit_num), fmax((*h2d_bw), (*d2h_bw)));
	if (unit_num) (*h2d_bw) = Gval_per_s(dim*dim*elemSize*unit_num, (*h2d_bw));
	if (rev_unit_num) (*d2h_bw) = Gval_per_s(rev_dim*rev_dim*elemSize*rev_unit_num, (*d2h_bw));
#ifdef PDEBUG
	fprintf(stderr, "Ran %d itterations for convergence.\n"
		"-> dim = %d, unit_ids = %s :\t\t h2d_bw(total) = %.2lf Gb/s\n"
		"-> rev_dim = %d, rev_unit_ids = %s :\t\t d2h_bw(total) = %.2lf Gb/s\n-> bidirectional_bw(total) = %.2lf Gb/s\n",
		sample_sz,
		dim, printlist<int>(unit_ids, unit_num), (*h2d_bw),
		rev_dim, printlist<int>(rev_unit_ids, rev_unit_num), (*d2h_bw), (*bid_bw));
	fprintf(stderr,"------------------------------------------------------------------------------------------------------------------------------------------------------\n");
#endif
	return csecond() - bench_t;
}
