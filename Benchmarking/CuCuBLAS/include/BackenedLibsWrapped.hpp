///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief The external wrapper for CoCoPelia + wrapped cuBLASXt
///
#ifndef BACKENDLIBSWRAPPED_H
#define BACKENDLIBSWRAPPED_H
#include <cuda_fp16.h>

/// cuBLASXt wrappers of Dgemm
double cuBLASXtDgemmWrap(char TransA,  char TransB, long int M, long int N, long int K,
  double alpha, double* A, long int ldA, double* B, long int ldB, double beta, double* C,
  long int ldC, long int T, double cpu_ratio, int dev_num, int dev_ids[]);

/// cuBLASXt wrappers of Sgemm
double cuBLASXtSgemmWrap(char TransA,  char TransB, long int M, long int N, long int K,
  float alpha, float* A, long int ldA, float* B, long int ldB, float beta, float* C,
  long int ldC, long int T, double cpu_ratio, int dev_num, int dev_ids[]);

/// cuBLAS wrappers of Hgemm -> NOT for performance comparisson, only one GPU + transfers!!!!
double cuBLASHgemmWrap(char TransA,  char TransB, long int M, long int N, long int K,
  __half alpha, __half* A, long int ldA, __half* B, long int ldB, __half beta, __half* C, long int ldC, int dev_id);

/// cuBLAS wrappers of Daxpy
double cuBLASDaxpyWrap(long int N, double alpha, double* x, long int incx,
  double* y, long int incy, double cpu_ratio, int dev_num, int dev_ids[]);

/// cuBLAS wrappers of Saxpy
double cuBLASDaxpyWrap(long int N, float alpha, float* x, long int incx,
  float* y, long int incy, double cpu_ratio, int dev_num, int dev_ids[]);

/// cuBLAS wrappers of Ddot
double cuBLASDdotWrap(long int N, double* x, long int incx,
  double* y, long int incy, double* result, double cpu_ratio, int dev_num, int dev_ids[]);

double cuBLASDgemvWrap(char TransA, long int M, long int N, double alpha, double *A, long int lda,
		double* x, long int incx, double beta, double* y, long int incy, double cpu_ratio, int dev_num, int dev_ids[] );

#endif
