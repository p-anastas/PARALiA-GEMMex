///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief The headers for functions for general use throught SKPeLia
///

#ifndef UNIHELPERS_BACK_H
#define UNIHELPERS_BACK_H

#include <cuda.h>
#include <cublas_v2.h>

#include "chl_smart_wrappers.hpp"
#include <atomic>

template<typename VALUETYPE> class gemm_backend_in{
public:
	char TransA,  TransB;
	int M, N, K, ldA, ldB, ldC;
	VALUETYPE alpha,beta;
	void **A, **B, **C;
	void* A_tile_v, *B_tile_v, *C_tile_v;
	short dev_id;
};

template<typename VALUETYPE> class gemv_backend_in{
public:
	char TransA,  incx, incy;
	int M, N, ldA;
	VALUETYPE alpha,beta;
	void **A, **x, **y;
	short dev_id;
};

template<typename VALUETYPE> class axpy_backend_in{
public:
		int N, incx, incy;
	VALUETYPE alpha;
	void **x, **y;
	short dev_id;
};

template<typename VALUETYPE> class slaxpby_backend_in{
public:
	int N, incx, incy;
	int slide_x, slide_y;
	VALUETYPE alpha, beta;
	void **x, **y;
	short dev_id;
};

template<typename VALUETYPE> class axpby_backend_in{
public:
	int N, incx, incy;
	VALUETYPE alpha, beta;
	void **x, **y;
	short dev_id;
};

template<typename VALUETYPE> class dot_backend_in{
public:
	int N, incx, incy;
	void **x, **y;
	VALUETYPE* result;
	short dev_id;
};

#include <cblas.h>

void TransposeTranslate(char TransChar, CBLAS_TRANSPOSE* cblasFlag, cublasOperation_t* cuBLASFlag, long int* ldim, long int dim1, long int dim2);

cublasOperation_t OpCblasToCublas(CBLAS_TRANSPOSE src);
CBLAS_TRANSPOSE OpCublasToCblas(cublasOperation_t src);
cublasOperation_t OpCharToCublas(char src);
CBLAS_TRANSPOSE OpCharToCblas(char src);
char PrintCublasOp(cublasOperation_t src);


/// Internally used utils TODO: Is this the correct way softeng wise?
void cudaCheckErrors();

// Lock wrapped_lock. This functions is fired in a queue to lock when it reaches that point.
void SKQueueLock(void* wrapped_lock);
// Unlock wrapped_lock. This functions is fired in a queue to unlock when it reaches that point.
void SKQueueUnlock(void* wrapped_lock);

// Struct containing an int pointer
typedef struct Ptr_atomic_int{
	std::atomic<int>* ato_int_ptr;
}* Ptr_atomic_int_p;
void SKIncAsync(void* wrapped_ptr_int);
void SKDecAsync(void* wrapped_ptr_int);

// Struct containing an int pointer and an int for Asynchronous set
typedef struct Ptr_and_int{
	int* int_ptr;
	int val;
}* Ptr_and_int_p;
void SKSetInt(void* wrapped_ptr_and_val);

// Struct containing a void pointer and a void for Asynchronous set
typedef struct Ptr_and_parent{
	void** ptr_parent;
	void* ptr_val;
}* Ptr_and_parent_p;
void SKSetPtr(void* wrapped_ptr_and_parent);

void SKSetTimerAsync(void* wrapped_timer_Ptr);

void SKFreeAllocAsync(void* backend_data);

void cublas_wrap_ddot(void* backend_data, void* queue_wrap_p);
void cublas_wrap_daxpy(void* backend_data, void* queue_wrap_p);
void cublas_wrap_saxpy(void* backend_data, void* queue_wrap_p);
void cublas_wrap_daxpby(void* backend_data, void* queue_wrap_p);

void cublas_wrap_dgemv(void* backend_data, void* queue_wrap_p);

void cublas_wrap_dgemm(void* backend_data, void* queue_wrap_p);
void cublas_wrap_sgemm(void* backend_data, void* queue_wrap_p);

void cblas_wrap_ddot(void* backend_data);
void cblas_wrap_daxpy(void* backend_data);
void cblas_wrap_saxpy(void* backend_data);
void cblas_wrap_daxpby(void* backend_data);


void cblas_wrap_dgemv(void* backend_data);

void cblas_wrap_dgemm(void* backend_data);
void cblas_wrap_sgemm(void* backend_data);

//void custom_cpu_wrap_dslaxpby(void* backend_data);
//void custom_avx2_cpu_wrap_dslaxpby(void* backend_data);

void custom_gpu_wrap_dslaxpby(void* backend_data, CQueue_p run_queue);

#endif
