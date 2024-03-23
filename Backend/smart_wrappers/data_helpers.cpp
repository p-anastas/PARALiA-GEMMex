///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief Some CUDA function calls with added error-checking
///

#include "smart_wrappers.hpp"
#include <omp.h>

template<typename VALUETYPE>
void CHLParallelVecInitHost(VALUETYPE *vec, long long length, int seed)
{
	srand(seed);
	#pragma omp parallel for
	for (long long i = 0; i < length; i++) vec[i] = (VALUETYPE) Drandom();
}

template void CHLParallelVecInitHost<double>(double *vec, long long length, int seed);
template void CHLParallelVecInitHost<float>(float *vec, long long length, int seed);