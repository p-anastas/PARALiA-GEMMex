///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief BLAS lvl 3 wrappers for benchmarks.
///
#include <cassert>
#include <cublasXt.h>
#include <cblas.h>
#include <cuda_runtime.h>

#include "PARALiA.hpp"
#include "backend_wrappers.hpp"

cublasXtHandle_t handle0 = NULL;
int handle_T = 0; 

void cblas_dgemm_wrap_for_cublasXt(char* gpu_op_A, char* gpu_op_B, int* M, int* N, int* K, double* alpha, double* A, int* ldA, double* B, int* ldB, double* beta, double* C, int* ldC){
  CBLAS_TRANSPOSE cpu_op_A = CblasNoTrans, cpu_op_B = CblasNoTrans;    // CblasNoTrans, CblasTrans

    //fprintf(stderr, "%d %d %d %lf %d %d %lf %d\n",*M, *N, *K, *alpha, *ldA, *ldB, *beta, *ldC);

    if(*gpu_op_A == 'N') cpu_op_A = CblasNoTrans;
    else if(*gpu_op_A == 'T') cpu_op_A = CblasTrans;
    else error("cblas_dgemm_wrap -> Invalid CUBLAS_OP for A");
    if(*gpu_op_B == 'N') cpu_op_B = CblasNoTrans;
    else if(*gpu_op_B == 'T') cpu_op_B = CblasTrans;
    else error("cblas_dgemm_wrap -> Invalid CUBLAS_OP for B");

    cblas_dgemm(CblasColMajor, cpu_op_A, cpu_op_B, *M, *N, *K, *alpha, A, *ldA, B, *ldB, *beta, C, *ldC);
}

double cuBLASXtDgemmWrap(char TransA, char TransB, long int M, long int N, long int K, double alpha, double* A, long int ldA, double* B, long int ldB, double beta, double* C, long int ldC, long int T, double cpu_ratio, int dev_num, int dev_ids[] ){
	int lvl = 1;
	double total_t = csecond();
#ifdef DEBUG
	fprintf(stderr, "|-----> cuBLASXtDgemmWrap(%c,%c,%zu,%zu,%zu,%lf,A(%d),%zu,B(%d),%zu,%lf,C(%d),%zu)\n",
		TransA, TransB, M, N, K, alpha, CHLGetPtrLoc(A), ldA,
		CHLGetPtrLoc(B), ldB, beta, CHLGetPtrLoc(C), ldC);
#endif

#ifdef TEST
	fprintf(stderr, "|-----> cuBLASXtDgemmWrap\n");
	double cpu_timer = csecond();
#endif

	cublasOperation_t gpu_op_A = OpCharToCublas(TransA), gpu_op_B = OpCharToCublas(TransB);
	cublasStatus_t stat;

	/// Required allocations for device
	if(!handle0 || handle_T!=T)
	{
		if(handle0) cublasXtDestroy(handle0);
		assert(CUBLAS_STATUS_SUCCESS == cublasXtCreate(&handle0));
		assert(CUBLAS_STATUS_SUCCESS == cublasXtDeviceSelect(handle0, dev_num, dev_ids));
		assert(CUBLAS_STATUS_SUCCESS == cublasXtSetBlockDim(handle0, T));
		assert(CUBLAS_STATUS_SUCCESS == cublasXtSetCpuRoutine(handle0, CUBLASXT_GEMM, CUBLASXT_DOUBLE, (void*) &cblas_dgemm_wrap_for_cublasXt));
		assert(CUBLAS_STATUS_SUCCESS == cublasXtSetCpuRatio(handle0, CUBLASXT_GEMM, CUBLASXT_DOUBLE, cpu_ratio));
		//assert(CUBLAS_STATUS_SUCCESS == cublasXtSetPinningMemMode(handle0, CUBLASXT_PINNING_ENABLED));
	}
#ifdef TEST
	cpu_timer = csecond() - cpu_timer;
	fprintf(stderr, "cuBLASXt initialization/pinning -> t_init = %lf ms\n", cpu_timer*1000);
    	cpu_timer = csecond();
#endif
	assert(CUBLAS_STATUS_SUCCESS == cublasXtDgemm(handle0, gpu_op_A, gpu_op_B, M, N, K, &alpha, A, ldA, B, ldB, &beta, C, ldC));
	//CHLSyncCheckErr();
#ifdef TEST
	cpu_timer = csecond() - cpu_timer;
	fprintf(stderr, "cuBLASXt execution time -> t_kernel = %lf ms\n", cpu_timer*1000);
#endif

	/// Free local buffers
	//cublasXtDestroy(handle0);
	//CHLSyncCheckErr();
	total_t = csecond() - total_t;
	return total_t;

}

void cblas_sgemm_wrap_for_cublasXt(char* gpu_op_A, char* gpu_op_B, int* M, int* N, int* K, float* alpha, float* A, int* ldA, float* B, int* ldB, float* beta, float* C, int* ldC){
  CBLAS_TRANSPOSE cpu_op_A = CblasNoTrans, cpu_op_B = CblasNoTrans;    // CblasNoTrans, CblasTrans

    //fprintf(stderr, "%d %d %d %lf %d %d %lf %d\n",*M, *N, *K, *alpha, *ldA, *ldB, *beta, *ldC);

    if(*gpu_op_A == 'N') cpu_op_A = CblasNoTrans;
    else if(*gpu_op_A == 'T') cpu_op_A = CblasTrans;
    else error("cblas_dgemm_wrap -> Invalid CUBLAS_OP for A");
    if(*gpu_op_B == 'N') cpu_op_B = CblasNoTrans;
    else if(*gpu_op_B == 'T') cpu_op_B = CblasTrans;
    else error("cblas_dgemm_wrap -> Invalid CUBLAS_OP for B");

cblas_sgemm(CblasColMajor, cpu_op_A, cpu_op_B, *M, *N, *K, *alpha, A, *ldA, B, *ldB, *beta, C, *ldC);
}

double cuBLASXtSgemmWrap(char TransA, char TransB, long int M, long int N, long int K, float alpha, float* A, long int ldA, float* B, long int ldB, float beta, float* C, long int ldC, long int T, double cpu_ratio, int dev_num, int dev_ids[] ){
	int lvl = 1;
	double total_t = csecond();
#ifdef DEBUG
	fprintf(stderr, "|-----> cuBLASXtSgemmWrap(%c,%c,%zu,%zu,%zu,%lf,A(%d),%zu,B(%d),%zu,%lf,C(%d),%zu)\n",
		TransA, TransB, M, N, K, alpha, CHLGetPtrLoc(A), ldA,
		CHLGetPtrLoc(B), ldB, beta, CHLGetPtrLoc(C), ldC);
#endif

#ifdef TEST
	fprintf(stderr, "|-----> cuBLASXtSgemmWrap\n");
	double cpu_timer = csecond();
#endif

	cublasOperation_t gpu_op_A = OpCharToCublas(TransA), gpu_op_B = OpCharToCublas(TransB);
	cublasStatus_t stat;
	cublasXtHandle_t handle0;

	/// Required allocations for device
	assert(CUBLAS_STATUS_SUCCESS == cublasXtCreate(&handle0));
	assert(CUBLAS_STATUS_SUCCESS == cublasXtDeviceSelect(handle0, dev_num, dev_ids));
	assert(CUBLAS_STATUS_SUCCESS == cublasXtSetBlockDim(handle0, T));
	assert(CUBLAS_STATUS_SUCCESS == cublasXtSetCpuRoutine(handle0, CUBLASXT_GEMM, CUBLASXT_DOUBLE, (void*) &cblas_dgemm_wrap_for_cublasXt));
	assert(CUBLAS_STATUS_SUCCESS == cublasXtSetCpuRatio(handle0, CUBLASXT_GEMM, CUBLASXT_DOUBLE, cpu_ratio));
	assert(CUBLAS_STATUS_SUCCESS == cublasXtSetPinningMemMode(handle0, CUBLASXT_PINNING_ENABLED));
#ifdef TEST
	cpu_timer = csecond() - cpu_timer;
	fprintf(stderr, "cuBLASXt initialization/pinning -> t_init = %lf ms\n", cpu_timer*1000);
    	cpu_timer = csecond();
#endif
	assert(CUBLAS_STATUS_SUCCESS == cublasXtSgemm(handle0, gpu_op_A, gpu_op_B, M, N, K, &alpha, A, ldA, B, ldB, &beta, C, ldC));
	CHLSyncCheckErr();
#ifdef TEST
	cpu_timer = csecond() - cpu_timer;
	fprintf(stderr, "cuBLASXt execution time -> t_kernel = %lf ms\n", cpu_timer*1000);
#endif

	/// Free local buffers
	cublasXtDestroy(handle0);
	CHLSyncCheckErr();
	total_t = csecond() - total_t;
	return total_t;
}

double cuBLASHgemmWrap(char TransA,  char TransB, long int M, long int N, long int K,
  __half alpha, __half* A, long int ldA, __half* B, long int ldB, 
  __half beta, __half* C, long int ldC, int dev_id){
	double total_t = csecond();
#ifdef DEBUG
	fprintf(stderr, "|-----> cuBLASHgemmWrap(%c,%c,%zu,%zu,%zu,%lf,A(%d),%zu,B(%d),%zu,%lf,C(%d),%zu)\n",
		TransA, TransB, M, N, K, alpha, CHLGetPtrLoc(A), ldA,
		CHLGetPtrLoc(B), ldB, beta, CHLGetPtrLoc(C), ldC);
#endif

#ifdef TEST
	fprintf(stderr, "|-----> cuBLASHgemmWrap\n");
	double cpu_timer = csecond();
#endif

	cublasOperation_t gpu_op_A = OpCharToCublas(TransA), gpu_op_B = OpCharToCublas(TransB);
	cublasStatus_t stat;
	cublasHandle_t handle0;

	/// Required allocations for device
	assert(CUBLAS_STATUS_SUCCESS == cublasCreate(&handle0));
	__half* A_dev, *B_dev, *C_dev;
	A_dev = (__half*) CHLMalloc(M*K*sizeof(__half), dev_id, 0);
	B_dev = (__half*) CHLMalloc(N*K*sizeof(__half), dev_id, 0);
	C_dev = (__half*) CHLMalloc(M*N*sizeof(__half), dev_id, 1);
	CHLMemcpy(A_dev, A,  M * K *sizeof(__half), dev_id, CHLGetPtrLoc(A));
	CHLMemcpy(B_dev, B,  K * N *sizeof(__half), dev_id, CHLGetPtrLoc(B));
	CHLMemcpy(C_dev, C,  M * N *sizeof(__half), dev_id, CHLGetPtrLoc(C));

#ifdef TEST
	cpu_timer = csecond() - cpu_timer;
	fprintf(stderr, "cuBLASHgemmWrap initialization/copying -> t_init = %lf ms\n", cpu_timer*1000);
    	cpu_timer = csecond();
#endif
	assert(CUBLAS_STATUS_SUCCESS == cublasHgemm(handle0, gpu_op_A, gpu_op_B, M, N, K, &alpha, A_dev, ldA, B_dev, ldB, &beta, C_dev, ldC));
	CHLSyncCheckErr();
#ifdef TEST
	cpu_timer = csecond() - cpu_timer;
	fprintf(stderr, "cuBLASHgemmWrap execution time -> t_kernel = %lf ms\n", cpu_timer*1000);
#endif

	CHLMemcpy(C, C_dev,  M * N *sizeof(__half), CHLGetPtrLoc(C), dev_id);
	CHLFree(A_dev, M*K*sizeof(__half), dev_id);
	CHLFree(B_dev, N*K*sizeof(__half), dev_id);
	CHLFree(C_dev, M*N*sizeof(__half), dev_id);
	/// Free local buffers
	cublasDestroy(handle0);
	CHLSyncCheckErr();
	total_t = csecond() - total_t;
	return total_t;

}