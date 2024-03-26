///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief Some CUDA function calls with added error-checking
///

#include "chl_smart_wrappers.hpp"
#include <omp.h>

template<typename VALUETYPE>
void CHLParallelVecInitHost(VALUETYPE *vec, long long length, int seed)
{
	srand(seed);
	//#pragma omp parallel for
	for (long long i = 0; i < length; i++) vec[i] = (VALUETYPE) Drandom();
}

template void CHLParallelVecInitHost<double>(double *vec, long long length, int seed);
template void CHLParallelVecInitHost<float>(float *vec, long long length, int seed);

long PAGE_sz = sysconf(_SC_PAGE_SIZE);

template<typename VALUETYPE>
void CHLTouche(VALUETYPE *vec, long long vec_length, int vec_elemSize)
{
	int elem_offset = PAGE_sz/vec_elemSize;
#ifndef PRODUCTION
	#pragma omp parallel
	{
		/* Obtain thread number */
 		int tid = omp_get_thread_num();
 		/* Only master thread does this */
 		if (tid == 0){
			int nthreads = omp_get_num_threads();
			fprintf(stderr, "CHLTouche(%p, %lld, %d): Using %d threads for touching memory\n", vec, vec_length, vec_elemSize, nthreads);
		}
	}
#endif
	#pragma omp parallel for
	for (long long i = 0; i < vec_length/elem_offset; i++) vec[i*elem_offset] = (VALUETYPE) Drandom();
}

template void CHLTouche<double>(double *vec, long long vec_length, int vec_elemSize);
template void CHLTouche<float>(float *vec, long long vec_length, int vec_elemSize);
