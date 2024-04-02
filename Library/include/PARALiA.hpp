///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief The external wrapper for PARALiA + wrapped cuBLASXt
///
#ifndef PARALiA_H
#define PARALiA_H

#include "Autotuner.hpp"
#include "Decomposer.hpp"
#include "DataCaching.hpp"
#include "Resource_manager.hpp"
#include <cuda_fp16.h>

typedef class ProblemMetadata{
public:	
	ATC_p autotuner;
	const char* problem_name;
	void* problem_wrap; 
	int decom_num;
	Decom2D_p decom[10];
	Buffer_p SAB[64] = {NULL}; 
	void print();
	~ProblemMetadata();
}* PMD_p; 

extern PMD_p PMD_cache[PROBLEM_MD_CACHE]; 
extern int PMD_cache_entries; 

/// The PARALiA Dgemm implementation.
ATC_p PARALiADgemm(char TransA,  char TransB, long int M, long int N, long int K,
	double alpha, double* A, long int ldA, double* B, long int ldB, double beta, double* C, long int ldC);

/// A modification of PARALiADgemm but with a given T (mainly for performance/debug purposes)
ATC_p PARALiADgemmControled(char TransA,  char TransB, long int M, long int N, long int K,
	double alpha, double* A, long int ldA, double* B, long int ldB, double beta, double* C, long int ldC, ATC_p predef_control_values);

/// The PARALiA Dgemm implementation.
ATC_p PARALiASgemm(char TransA,  char TransB, long int M, long int N, long int K,
	float alpha, float* A, long int ldA, float* B, long int ldB, float beta, float* C, long int ldC);

/// A modification of PARALiADgemm but with a given T (mainly for performance/debug purposes)
ATC_p PARALiASgemmControled(char TransA,  char TransB, long int M, long int N, long int K,
	float alpha, float* A, long int ldA, float* B, long int ldB, float beta, float* C, long int ldC, ATC_p predef_control_values);

/// The PARALiA Dgemm implementation.
ATC_p PARALiAHgemm(char TransA,  char TransB, long int M, long int N, long int K,
	__half alpha, __half* A, long int ldA, __half* B, long int ldB, __half beta, __half* C, long int ldC);

/// A modification of PARALiADgemm but with a given T (mainly for performance/debug purposes)
ATC_p PARALiAHgemmControled(char TransA,  char TransB, long int M, long int N, long int K,
	__half alpha, __half* A, long int ldA, __half* B, long int ldB, __half beta, __half* C, long int ldC, ATC_p predef_control_values);


///Deallocates the GPU-allocated cache buffer at target device
void PARALiADevCacheFree(int dev_id);

#endif
