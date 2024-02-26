///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief The header for the CoCoPeLia Microbenchmarks (Deployment_phase).
///

#ifndef MICROBENCH_H
#define MICROBENCH_H

/// Check for BOOST use (enabled from CmakeFiles) or number of ITER for microbenchmarks.
#ifdef AUTO_BENCH_USE_BOOST
#include <boost/math/distributions/students_t.hpp>
#define alphaCI 0.05
#define MICRO_MIN_ITER 10
#define MICRO_MAX_ITER 100000
#else
#ifndef ITER
#error
#endif
#endif

/// The minimum dimensions and step sizes (of Dim*Dim chunks) for microbenchmarks
#define MIN_DIM_TRANS 256
#define STEP_TRANS 256
#define MAX_DIM_TRANS 8192
// Define how many simultaneous transfers should be used to emulate interconnect load. 
#define INTERCONNECT_LOAD_TRANSFERS 3

#define MIN_DIM_BLAS3 256
#define STEP_BLAS3 256
#define MAX_DIM_BLAS3 10000

#define MIN_DIM_BLAS2 256
#define STEP_BLAS2 256
#define MAX_DIM_BLAS2 10000

#define MIN_DIM_BLAS1 4096
#define STEP_BLAS1 4096

#include <stdio.h>
#include <cstring>
#include <cuda.h>
#include "cublas_v2.h"

extern int use_square_Tiles;
extern int elemSize;
extern CQueue_p h2d_queue_list[32], d2h_queue_list[32];

int confidence_interval_5_percent(long int sample_sz, double cpu_timer, double* transfer_t_vals, double* transfer_t_sum_ptr, double* transfer_t_mean_ptr, double* error_margin_ptr);
double latbw_linear_regression(double* y, double* x, int samples, double* ans_b, double* ans_a, int active_unit_num);
void find_numa_connectivity(long long microbench_bytes, int* numa_high_affinity, int* numa_low_affinity, double *** nuMap_weights_h2d, double *** nuMap_weights_d2h, double *** nuMap_weights_bid);
double perform_microbenchmark_2D(void** loc_buffs, void** worker_buffs, int loc, long long total_bytes_per_side, int lddev, int ldhost, 
	int dim, int* unit_ids, int unit_num, double* worker_wise_bws, int rev_dim, int* rev_unit_ids, int rev_unit_num, double* worker_wise_rev_bws, double* h2d_bw, double* d2h_bw, double* bid_bw);

#endif
