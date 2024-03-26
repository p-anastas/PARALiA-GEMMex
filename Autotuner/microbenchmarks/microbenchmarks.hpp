///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief The header for the CoCoPeLia Microbenchmarks (Deployment_phase).
///

#ifndef MICROBENCH_H
#define MICROBENCH_H

/// Check for BOOST use (enabled from CmakeFiles) or number of ITER for microbenchmarks.
#include <boost/math/distributions/students_t.hpp>
#define alphaCI 0.05
#define MICRO_MIN_ITER 10
#define MICRO_MAX_ITER 10000

//#define MIN_DIM_BLAS3 256
//#define STEP_BLAS3 256
//#define MAX_DIM_BLAS3 10000

#include <stdio.h>
#include <cstring>
#include <cuda.h>
#include "cublas_v2.h"

int confidence_interval_5_percent(long int sample_sz, double cpu_timer, double* transfer_t_vals, double* transfer_t_sum_ptr, double* transfer_t_mean_ptr, double* error_margin_ptr);

#endif
