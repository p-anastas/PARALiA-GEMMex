///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief The DGEMM CoCopeLia implementation.
///

#include <cblas.h>

#include "Decomposer.hpp"

#include "backend_wrappers.hpp"

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <pthread.h>

/*
typedef struct pthread_data_in{
	void* adrs;
	long long pin_bytes;
	short pin_internally;
}* pthread_data_in_p;

void* prepareAsync_backend(void* compressed_data){
	pthread_data_in_p prep_data = (pthread_data_in_p)compressed_data;
	if (prep_data->pin_internally) cudaHostRegister(prep_data->adrs,prep_data->pin_bytes,cudaHostRegisterPortable);
#ifdef DEBUG
	fprintf(stderr,"pin asset= %d\n", prep_data->pin_internally);
#endif
	return NULL;
}
*/

void Decomposer::prepareAsync(){
	if(loc >= CHL_WORKERS) // && cudaErrorInvalidValue==cudaPointerGetAttributes(&attributes, adrs)) 
		pin_internally = 1;
	if (pin_internally){
		fprintf(stderr, "Pinning Matrix  %d\n", pin_internally);
		cudaError_t err;
		err = cudaHostRegister(adrs, get_mem_size(),cudaHostRegisterPortable);
		if(cudaSuccess != err) fprintf(stderr, "cudaHostRegister returned %s\n", cudaGetErrorString(err));
		//massert(cudaSuccess == err, "Decomposer::prepareAsync failed to pin memory with error %s\n", cudaGetErrorString(err));
		cudaGetLastError();
	}
}

void Decomposer::resetProperties(){
	if (pin_internally){
		cudaHostUnregister(adrs);
		cudaGetLastError();
	}
}
