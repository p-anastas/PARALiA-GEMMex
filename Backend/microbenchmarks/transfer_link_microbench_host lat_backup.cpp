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

/// TODO: This is the number of other-link transfers used for overlaping.
/// The benchmark results are correct ONLY if Link_t > MAX_ASSUMED_OTHER_LINK_TIMES_FASTER * Other_link_t.
/// For current systems 10 is sufficient - larger multipliers increase total benchmark time.
#define MAX_ASSUMED_OTHER_LINK_TIMES_FASTER 10

int use_square_Tiles = 0; 

double linear_regression(double* y, double* x, int samples, double* ans_b, double* ans_a, int active_unit_num){
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
		fprintf(stderr,"Latency-BW estimated breakpoint:\t b_eq = %d bytes\n", byte_breakpoint/active_unit_num);
		fprintf(stderr,"75% BW estimated breakpoint:\t\t b_75 = %d bytes\n", byte_breakpoint/active_unit_num*2);
		fprintf(stderr,"95% BW estimated breakpoint:\t\t b_95 = %d bytes\n\n", byte_breakpoint/active_unit_num*20);
	}
	else{
		fprintf(stderr,"Latency-BW estimated breakpoint:\t b_eq = 8-bytes X (%d X %d)\n", (int) sqrt(byte_breakpoint/active_unit_num/8), (int) sqrt(byte_breakpoint/active_unit_num/8));
		fprintf(stderr,"75% BW estimated breakpoint:\t\t b_75 = 8-bytes X (%d X %d)\n", (int) sqrt(byte_breakpoint/active_unit_num/4), (int) sqrt(byte_breakpoint/active_unit_num/4));
		fprintf(stderr,"95% BW estimated breakpoint:\t\t b_95 = 8-bytes X (%d X %d)\n\n", (int) sqrt(byte_breakpoint/active_unit_num*20/8), (int) sqrt(byte_breakpoint/active_unit_num*20/8));
	}
	double APE = 0; 
	long long byte_lbw_bp = 0, byte_75pbw_bp = 0, byte_95pbw_bp = 0;
	for(int i=0; i<samples; i++){
		double pred_y =  t_lat + x[i]*t_b; 
#ifdef CLDEBUG
		fprintf(stderr, "Point (x = %d):\t y=%lf ms (%lf Gb/s), pred_y = %lf ms (%lf Gb/s) PE = %.1lf\n", 
		int(x[i]), y[i] * 1000, Gval_per_s(int(x[i]), y[i]), pred_y * 1000, Gval_per_s(int(x[i]), pred_y), (pred_y - y[i])/y[i]*100);
#endif
		APE += abs((pred_y - y[i])/y[i]*100);
		if (!byte_lbw_bp && y[i]/x[i] <= (*ans_a)/0.5) byte_lbw_bp = x[i];
		if (!byte_75pbw_bp && y[i]/x[i] <= (*ans_a)/0.75) byte_75pbw_bp = x[i];
		if (!byte_95pbw_bp && y[i]/x[i] <= (*ans_a)/0.95) byte_95pbw_bp = x[i];

	}
	if (!use_square_Tiles){
		fprintf(stderr,"Latency-BW empirical breakpoint:\t b_eq = %d bytes\n", byte_lbw_bp/active_unit_num);
		fprintf(stderr,"75% BW empirical breakpoint:\t\t b_75 = %d bytes\n", byte_75pbw_bp/active_unit_num);
		fprintf(stderr,"95% BW empirical breakpoint:\t\t b_95 = %d bytes\n\n", byte_95pbw_bp/active_unit_num);
	}
	else{
		fprintf(stderr,"Latency-BW empirical breakpoint:\t b_eq = 8-bytes X (%d X %d)\n", (int) sqrt(byte_lbw_bp/active_unit_num/8), (int) sqrt(byte_lbw_bp/active_unit_num/8));
		fprintf(stderr,"75% BW empirical breakpoint:\t\t b_75 = 8-bytes X (%d X %d)\n", (int) sqrt(byte_75pbw_bp/active_unit_num/8), (int) sqrt(byte_75pbw_bp/active_unit_num/8));
		fprintf(stderr,"95% BW empirical breakpoint:\t\t b_95 = 8-bytes X (%d X %d)\n\n", (int) sqrt(byte_95pbw_bp/active_unit_num/8), (int) sqrt(byte_95pbw_bp/active_unit_num/8));
	}

	return APE/samples; 
}
int confidence_interval_5_percent(long int sample_sz, double cpu_timer, double transfer_t_vals[MICRO_MAX_ITER], double* transfer_t_sum_ptr, double* transfer_t_mean_ptr, double* error_margin_ptr){
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
	//fprintf(stderr, "\tItter %d:\t mean=%lf, std_dev = %lf, Error margin =%lf\n", sample_sz, (*transfer_t_mean_ptr) , std_dev, (*error_margin_ptr));
	if ((*error_margin_ptr)/(*transfer_t_mean_ptr)  * 100 <= 5) return 1;
	else return 0;
}
int main(const int argc, const char *argv[]) {

  int ctr = 1, samples, dev_id, dev_count;

	short numa_aware_flag = 0;
	int minDim = MIN_DIM_TRANS, maxDim = 0, step = STEP_TRANS;
	if (minDim < 1) error("Transfer Microbench: Bytes must be > 0");

	switch (argc) {
	case (3):
		numa_aware_flag = atoi(argv[ctr++]);
		use_square_Tiles = atoi(argv[ctr++]);
		break;
	default:
		error("Incorrect input arguments. Usage: ./correct_run numa_aware_flag(=0 or 1) use_square_Tiles(=0 or 1)\n");
  	}

	char *filename = (char *) malloc(1024 * sizeof(char));
	sprintf(filename, "%s/Benchmark-Results/transfer_link_microbenchmark_host_%s.log", DEPLOYDB, VERSION);
	//check_benchmark(filename);

	short elemSize = sizeof(double);

	// Define the max size of a benchmark kernel to run on this machine.
	maxDim = std::min(MAX_DIM_TRANS, (int) CoCoGetMaxDimSqAsset2D(2*(LOC_NUM-1), sizeof(double), STEP_TRANS, -1));
	
	for(short dev_id_idx = 0 ; dev_id_idx < LOC_NUM; dev_id_idx++){
		maxDim = std::min(maxDim, (int)CoCoGetMaxDimSqAsset2D(2, sizeof(double), STEP_TRANS, deidxize(dev_id_idx)));
	}
  	int ldhost = std::min((int) (1.4*MAX_DIM_TRANS), (int) CoCoGetMaxDimSqAsset2D(2*2*(LOC_NUM-1), sizeof(double), STEP_TRANS, -1)),
		lddev = maxDim;
	fprintf(stderr,"\ntransfer_link_microbench_host: \nSystem = %s\nminDim = %d, maxDim = %d, step = %d(adaptive)\nldhost = %d, lddev = %d\n", 
		TESTBED, minDim, maxDim, step, ldhost, lddev);
	fprintf(stderr,"----------------------------------------------------------------------\n");

	void* host_buffs[(LOC_NUM-1)][2];
	void* dev_buffs[(LOC_NUM-1)][2];
	double timer = csecond();
	for(short dev_id_idx = 0 ; dev_id_idx < LOC_NUM -1; dev_id_idx++){
  		host_buffs[dev_id_idx][0] = CoCoMalloc(ldhost*ldhost*elemSize, -1, 1);
  		host_buffs[dev_id_idx][1] = CoCoMalloc(ldhost*ldhost*elemSize, -1, 1);
	}
	timer = csecond() - timer;
	fprintf(stderr, "Allocation host_buffs size = (%d x %d) x %d complete:\t alloc_timer=%lf ms\n", ldhost, ldhost, elemSize, timer  * 1000);

	timer = csecond();
	for(short dev_id_idx = 0 ; dev_id_idx < LOC_NUM -1; dev_id_idx++){
		dev_buffs[dev_id_idx][0] = CoCoMalloc(lddev*lddev*elemSize, deidxize(dev_id_idx), 1);
		dev_buffs[dev_id_idx][1] = CoCoMalloc(lddev*lddev*elemSize, deidxize(dev_id_idx), 1);
	}
	timer = csecond() - timer;
	fprintf(stderr, "Allocation dev_buffs size = (%d x %d) x %d complete:\t alloc_timer=%lf ms\n", lddev, lddev, elemSize, timer  * 1000);
	fprintf(stderr,"----------------------------------------------------------------------\n");

	CQueue_p h2d_queue_list[LOC_NUM-1] = {NULL}, d2h_queue_list[LOC_NUM-1] = {NULL};
	for(short dev_id_idx = 0 ; dev_id_idx < LOC_NUM-1; dev_id_idx++){
		//printf("dev_id = %d, dev_id_idx = %d, dev_id_idy = %d, LOC_NUM = %d\n", dev_id, dev_id_idx, dev_id_idy, LOC_NUM);
		short queue_id = deidxize(dev_id_idx);
		h2d_queue_list[dev_id_idx] = new CommandQueue(queue_id, 0);
		d2h_queue_list[dev_id_idx] = new CommandQueue(queue_id, 0);
	}

	fprintf(stderr, "Warming up");
	/// Warmup.
	for (int it = 0; it < 3; it++){
		fprintf(stderr, ".");
		for(short dev_id_idx = 0 ; dev_id_idx < LOC_NUM-1; dev_id_idx++){
		CoCoMemcpy2DAsync(dev_buffs[dev_id_idx][0], lddev,
			host_buffs[dev_id_idx][0], ldhost,
			maxDim, maxDim, elemSize,
			deidxize(dev_id_idx), -1, h2d_queue_list[dev_id_idx]);
			h2d_queue_list[dev_id_idx]->sync_barrier();
		}
		for(short dev_id_idx = 0 ; dev_id_idx < LOC_NUM-1; dev_id_idx++){
			CoCoMemcpy2DAsync(host_buffs[dev_id_idx][1], ldhost,
				dev_buffs[dev_id_idx][1], lddev,
				maxDim, maxDim, elemSize,
				-1, deidxize(dev_id_idx), d2h_queue_list[dev_id_idx]);
			d2h_queue_list[dev_id_idx]->sync_barrier();
		}
	}
	fprintf(stderr, " complete.\n");
	CoCoSyncCheckErr();

	int loc_src = -1;
	int active_unit_num, active_unit_id_list[LOC_NUM-1] = {-42}, explored_cases = pow(2,(LOC_NUM-1));
	int best_case_id_per_devnum[LOC_NUM-1] = {-42};
	double best_bw_per_devnum[LOC_NUM-1] = {-42}, best_lat_per_devnum[LOC_NUM-1] = {-42};
	for (int case_id = 1; case_id < explored_cases; case_id++){
		translate_binary_to_unit_list(case_id, &active_unit_num, active_unit_id_list);
		int dim;
		step = STEP_TRANS;
		double MB_t[int((maxDim-minDim)/step) +1], MB_dim[int((maxDim-minDim)/step) + 1];
		int MB_num = 0;
		double cpu_timer, transfer_t_vals[MICRO_MAX_ITER] = {0}, transfer_t_sum = 0, transfer_t_mean = 0, error_margin = 0, bench_t = csecond();
		double transfer_t_bid_sum, transfer_t_bid_mean, error_margin_bid;
		long int sample_sz, sample_sz_bid;
		for (dim = minDim; dim < maxDim; dim+=step){ // maxDim+1
			if (dim >= step * 8) step*=2;
			if (dim > maxDim) break;
			transfer_t_sum = transfer_t_mean = bench_t = error_margin = 0;
			//fprintf(stderr, "Cublas-chunk Link %d->%d (Chunk %dx%d):\n", loc_src, loc_dest, dim, dim);
			sample_sz = 0;
			for (sample_sz = 1; sample_sz < MICRO_MAX_ITER + 1; sample_sz++){
				cpu_timer = csecond();
					for(int active_unit_idx = 0; active_unit_idx < active_unit_num; active_unit_idx++){
					int loc_dest = active_unit_id_list[active_unit_idx];
					CoCoPeLiaSelectDevice(loc_dest);
					if (use_square_Tiles) CoCoMemcpy2DAsync(dev_buffs[idxize(loc_dest)][0], lddev,
						host_buffs[idxize(loc_dest)][0], ldhost, dim, dim, elemSize,
						loc_dest, loc_src, h2d_queue_list[idxize(loc_dest)]);
					else CoCoMemcpyAsync(dev_buffs[idxize(loc_dest)][0], host_buffs[idxize(loc_dest)][0],
						dim*dim*elemSize, loc_dest, loc_src, h2d_queue_list[idxize(loc_dest)]);
					}
					for(int active_unit_idx = 0; active_unit_idx < active_unit_num; active_unit_idx++) 
						h2d_queue_list[idxize(active_unit_id_list[active_unit_idx])]->sync_barrier();


				cpu_timer  = csecond() - cpu_timer;
				int complete_flag = confidence_interval_5_percent(sample_sz, cpu_timer, transfer_t_vals, &transfer_t_sum, &transfer_t_mean, &error_margin);
				if (sample_sz > MICRO_MIN_ITER && complete_flag) break;
			}
			MB_dim[MB_num] = 1.0*dim*dim*elemSize*active_unit_num;
			MB_t[MB_num++] = transfer_t_mean;
			CoCoSyncCheckErr();
			//fprintf(stderr, "Microbenchmark (dim1 = dim2 = %d) complete:\t mean_exec_t=%lf ms  ( %lf Gb/s), Error Margin (percentage of mean) = %lf %%"
			//", Itter = %ld\n\n", dim, transfer_t_mean  * 1000, Gval_per_s(dim*dim*8, transfer_t_mean), error_margin/transfer_t_mean  * 100, sample_sz);
		}
		double t_lat = 0, t_b = 0; 
		fprintf(stderr,"----------------------------------------------------------------------\n");
		double MAPE = linear_regression(MB_t, MB_dim, MB_num, &t_lat, &t_b, active_unit_num);
		bench_t = csecond() - bench_t;
		fprintf(stderr,"\nlinear_regression complete on %d samples for Host -> %s : t_lat = %e, t_b = %e, MAPE = %.1lf, bench_t=%lf\n", 
		MB_num, printlist<int>(active_unit_id_list, active_unit_num), t_lat, t_b, MAPE, bench_t);
		if (1/t_b >= best_bw_per_devnum[active_unit_num-1]*(1+NORMALIZE_NEAR_SPLIT_LIMIT) ||
			(1/t_b >= best_bw_per_devnum[active_unit_num-1]*(1-NORMALIZE_NEAR_SPLIT_LIMIT) && 
			t_lat < best_lat_per_devnum[active_unit_num-1])){
				best_bw_per_devnum[active_unit_num-1] = 1/t_b;
				best_lat_per_devnum[active_unit_num-1] = t_lat;
				best_case_id_per_devnum[active_unit_num-1] = case_id;
			}
		fprintf(stderr,"----------------------------------------------------------------------\n");
	}
	for(int idx = 0; idx < LOC_NUM - 1; idx++){
		fprintf(stderr,"----------------------------------------------------------------------\n\n");
		translate_binary_to_unit_list(best_case_id_per_devnum[idx], &active_unit_num, active_unit_id_list);
		fprintf(stderr, "Best device set for dev_num = %d -> %s : t_lat = %e, t_b = %e ( %.2lf Gb/s)\n", active_unit_num, 
		printlist<int>(active_unit_id_list, active_unit_num), best_lat_per_devnum[idx], 1/best_bw_per_devnum[idx], best_bw_per_devnum[idx]/1e9);
		fprintf(stderr,"----------------------------------------------------------------------\n\n");

	}
	/*
	double transfer_timer, ETA_timer;
	int share_matrix[LOC_NUM][LOC_NUM], block_matrix[LOC_NUM][LOC_NUM], share_matrix_idx = 0, block_matrix_idx = 0; 
	for(short dev_id_idx = 0 ; dev_id_idx < LOC_NUM; dev_id_idx++){
		for(short dev_id_idy = 0 ; dev_id_idy < LOC_NUM; dev_id_idy++){
			if (dev_id_idx == dev_id_idy ||
				(dev_id_idx == idxize(loc_dest) && dev_id_idy == idxize(loc_src))) continue;
			bench_t = csecond();
			for(int itt = 0 ; itt < MAX_ASSUMED_OTHER_LINK_TIMES_FASTER; itt++) CoCoMemcpy2DAsync(unit_buffs[2*dev_id_idx+1], ldest,
										unit_buffs[2*dev_id_idy+1], ldsrc,
										dim, dim, elemSize,
										deidxize(dev_id_idx), deidxize(dev_id_idy), transfer_queue_list[dev_id_idx][dev_id_idy]);
			device_timer->start_point(transfer_queue_list[idxize(loc_dest)][idxize(loc_src)]);
			cpu_timer = csecond();
			CoCoMemcpy2DAsync(unit_buffs[2*idxize(loc_dest)], ldest,
										unit_buffs[2*idxize(loc_src)], ldsrc,
										dim, dim, elemSize,
										loc_dest, loc_src, transfer_queue_list[idxize(loc_dest)][idxize(loc_src)]);
			device_timer->stop_point(transfer_queue_list[idxize(loc_dest)][idxize(loc_src)]);
			transfer_queue_list[idxize(loc_dest)][idxize(loc_src)]->sync_barrier();
			cpu_timer = csecond() - cpu_timer;
			ETA_timer = cpu_timer;
			transfer_timer = device_timer->sync_get_time()/1000;
			transfer_queue_list[dev_id_idx][dev_id_idy]->sync_barrier();
			CoCoSyncCheckErr();
			bench_t = csecond() - bench_t;
			//fprintf(stderr, "Shared Link (%d->%d) transfer complete:\t shared_timer=%lf ms  ( %lf Gb/s)\n\n",
			//	deidxize(dev_id_idy), deidxize(dev_id_idx), shared_timer  * 1000, Gval_per_s((dim)*(dim)*elemSize, shared_timer));

			if (transfer_t_mean < ETA_timer*(1-NORMALIZE_NEAR_SPLIT_LIMIT)){
				if (!(transfer_timer*(1-NORMALIZE_NEAR_SPLIT_LIMIT) < ETA_timer && ETA_timer > transfer_timer*(1+NORMALIZE_NEAR_SPLIT_LIMIT))){
					fprintf(stderr, "Link(%2d->%2d) & Link(%2d->%2d) partially shared: Shared_BW: %1.2lf %%\n\n",
					loc_src, loc_dest, deidxize(dev_id_idy), deidxize(dev_id_idx), 100*transfer_t_mean/transfer_timer);
					share_matrix[dev_id_idx][dev_id_idy] = 1; 
					share_matrix_idx++;
				}
				else{
					fprintf(stderr, "Link(%2d->%2d) blocked by Link(%2d->%2d) for %lf ms ( %lf %% of its transfer_t) and then ran with BW_ratio: %1.2lf %%\n\n",
					loc_src, loc_dest, deidxize(dev_id_idy), deidxize(dev_id_idx), 1000*(ETA_timer-transfer_timer), 100*(ETA_timer-transfer_timer)/transfer_timer, 100*transfer_t_mean/transfer_timer);
					block_matrix[dev_id_idx][dev_id_idy] = 1; 
					block_matrix_idx++;		
				}
			}
		}
	}

	fprintf(stderr,"----------------------------------------------------------------------\n");

	fprintf(stderr, "Link %d->%d shares resources with %d other links... %d combinations:\n", loc_src, loc_dest, share_matrix_idx, (int) pow(2, share_matrix_idx));

	fprintf(stderr,"----------------------------------------------------------------------\n");
*/
	CoCoSyncCheckErr();
	fprintf(stderr,"----------------------------------------------------------------------\n");
	timer = csecond();
	for(short dev_id_idx = 0 ; dev_id_idx < LOC_NUM -1; dev_id_idx++){
  		CoCoFree(host_buffs[dev_id_idx][0], ldhost*ldhost*elemSize, -1);
  		CoCoFree(host_buffs[dev_id_idx][1], ldhost*ldhost*elemSize, -1);
	}
	timer = csecond() - timer;
	fprintf(stderr, "Free host_buffs size = (%d x %d) x %d complete:\t alloc_timer=%lf ms\n", ldhost, ldhost, elemSize, timer  * 1000);

	timer = csecond();
	for(short dev_id_idx = 0 ; dev_id_idx < LOC_NUM -1; dev_id_idx++){
  		CoCoFree(dev_buffs[dev_id_idx][0], lddev*lddev*elemSize, deidxize(dev_id_idx));
  		CoCoFree(dev_buffs[dev_id_idx][1], lddev*lddev*elemSize, deidxize(dev_id_idx));
	}
	timer = csecond() - timer;
	fprintf(stderr, "Free dev_buffs size = (%d x %d) x %d complete:\t alloc_timer=%lf ms\n", lddev, lddev, elemSize, timer  * 1000);
	fprintf(stderr,"----------------------------------------------------------------------\n");
  	return 0;
}
