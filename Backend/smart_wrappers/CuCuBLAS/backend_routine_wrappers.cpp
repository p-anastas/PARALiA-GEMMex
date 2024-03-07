///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief Backened wrapped host functions for cuda queue firing -> void func (void*)
///

#include <cblas.h>
#include <omp.h>

#include "backend_wrappers.hpp"

void SKQueueLock(void* wrapped_lock){
#ifdef ENABLE_MUTEX_LOCKING
  (*(std::mutex*)wrapped_lock).lock();
#else
  while(__sync_lock_test_and_set ((&(*((int*)wrapped_lock))), 1));
#endif
#ifdef DEBUG
  fprintf(stderr, "SKQueueLock(%p) ran succesfully.\n", wrapped_lock);
#endif
}

void SKQueueUnlock(void* wrapped_lock){
#ifdef ENABLE_MUTEX_LOCKING
	(*(std::mutex*)wrapped_lock).unlock();
#else
  //int* intptr = (int*) wrapped_lock;
  //*intptr = 0;
  __sync_lock_release((&(*((int*) wrapped_lock))));
#endif

#ifdef DEBUG
  fprintf(stderr, "SKQueueUnlock(%p) ran succesfully.\n", wrapped_lock);
#endif
}

void SKIncAsync(void* wrapped_ptr_int){
  Ptr_atomic_int_p unwrapped = (Ptr_atomic_int_p) wrapped_ptr_int;
  *(unwrapped->ato_int_ptr)++;
  free(unwrapped);
#ifdef DEBUG
  fprintf(stderr, "SKIncAsync(%p, new_val=%d) ran succesfully.\n", unwrapped->ato_int_ptr, (*(unwrapped->ato_int_ptr)).load());
#endif
}

void SKDecAsync(void* wrapped_ptr_int){
  Ptr_atomic_int_p unwrapped = (Ptr_atomic_int_p) wrapped_ptr_int;
  (*(unwrapped->ato_int_ptr))--;
  free(unwrapped);
#ifdef DEBUG
  fprintf(stderr, "SKDecAsync(%p, new_val=%d) ran succesfully.\n", unwrapped->ato_int_ptr, (*(unwrapped->ato_int_ptr)).load());
#endif
}

void SKSetInt(void* wrapped_ptr_and_val){
  Ptr_and_int_p unwrapped = (Ptr_and_int_p) wrapped_ptr_and_val;
  *(unwrapped->int_ptr) = unwrapped->val;
  free(unwrapped);
#ifdef DEBUG
  fprintf(stderr, "SKSetVal(%p, %d) ran succesfully.\n", unwrapped->int_ptr, unwrapped->val);
#endif
}

void SKSetPtr(void* wrapped_ptr_and_parent){
  Ptr_and_parent_p unwrapped = (Ptr_and_parent_p) wrapped_ptr_and_parent;
  void* prev_ptr = *(unwrapped->ptr_parent);
  *(unwrapped->ptr_parent) = unwrapped->ptr_val;
  free(unwrapped);
#ifdef DEBUG
  fprintf(stderr, "SKSetPtr(prev=%p, %p) ran succesfully.\n", prev_ptr, unwrapped->ptr_val);
#endif
}

void cblas_wrap_daxpy(void* backend_data){
  axpy_backend_in<double>* ptr_ker_translate = (axpy_backend_in<double>*) backend_data;
#ifdef DEBUG
  fprintf(stderr,"cblas_wrap_daxpy: cublasDaxpy(dev_id = %d,\
    N = %d, alpha = %lf, x = %p, incx = %d, y = %p, incy = %d)\n",
    ptr_ker_translate->dev_id, ptr_ker_translate->N, ptr_ker_translate->alpha,
    (double*) *ptr_ker_translate->x, ptr_ker_translate->incx,
    (double*) *ptr_ker_translate->y, ptr_ker_translate->incy);
#endif
  cblas_daxpy(ptr_ker_translate->N, ptr_ker_translate->alpha,
    (double*) *ptr_ker_translate->x, ptr_ker_translate->incx, (double*)
    *ptr_ker_translate->y, ptr_ker_translate->incy);
}

void cblas_wrap_daxpby(void* backend_data){
  axpby_backend_in<double>* ptr_ker_translate = (axpby_backend_in<double>*) backend_data;
#ifdef DEBUG
  fprintf(stderr,"cblas_wrap_daxpby: cblas_daxpby(dev_id = %d,\
    N = %d, alpha = %lf, x = %p, incx = %d, b = %lf, y = %p, incy = %d)\n",
    ptr_ker_translate->dev_id, ptr_ker_translate->N, ptr_ker_translate->alpha,
    (double*) *ptr_ker_translate->x, ptr_ker_translate->incx, ptr_ker_translate->beta,
    (double*) *ptr_ker_translate->y, ptr_ker_translate->incy);
#endif
  cblas_daxpby(ptr_ker_translate->N, ptr_ker_translate->alpha,
    (double*) *ptr_ker_translate->x, ptr_ker_translate->incx, ptr_ker_translate->beta,
    (double*) *ptr_ker_translate->y, ptr_ker_translate->incy);
}

void* custom_cpu_wrap_dslaxpby_pthread_wrap(void* backend_data){
  slaxpby_backend_in<double>* ptr_ker_translate = (slaxpby_backend_in<double>*) backend_data;
#ifdef DEBUG
  fprintf(stderr,"custom_cpu_wrap_dslaxpby_pthread_wrap(dev_id = %d,\
    N = %d, alpha = %lf, x = %p, incx = %d, b = %lf, y = %p, incy = %d, slide_x = %d, slide_y = %d)\n",
    ptr_ker_translate->dev_id, ptr_ker_translate->N, ptr_ker_translate->alpha,
    (double*) *ptr_ker_translate->x, ptr_ker_translate->incx, ptr_ker_translate->beta,
    (double*) *ptr_ker_translate->y, ptr_ker_translate->incy, ptr_ker_translate->slide_x, ptr_ker_translate->slide_y);
#endif
	int cpu_aff = CHL_HWNUMA_AT_MEMLOC[ptr_ker_translate->dev_id];//myself->dev_id; //
	cpu_set_t cpuset;
	CPU_ZERO(&cpuset);
	// Set thread affinity to the first core in its 'friend' CPU.
	// For now, always use core 0 from each CPU for queue worker threads.
	int core_aff = cpu_aff*HW_THREADS/NUMA_HW_NUM + 0;
	//for(int idx = 0; idx < HW_THREADS/NUMA_HW_NUM; idx++){
	//	core_aff = cpu_aff*HW_THREADS/NUMA_HW_NUM + idx;
		CPU_SET(core_aff, &cpuset);
	//}
	pthread_t curr_thread = pthread_self();
	pthread_setaffinity_np(curr_thread, sizeof(cpu_set_t), &cpuset);
  double* y = (double*) *ptr_ker_translate->y, *x = (double*) *ptr_ker_translate->x;
	int i, j, N = ptr_ker_translate->N, offset_x = ptr_ker_translate->slide_x, offset_y = ptr_ker_translate->slide_y; 
	double alpha = ptr_ker_translate->alpha, beta = ptr_ker_translate->beta;
  if(ptr_ker_translate->alpha != 1.0) error("custom_avx2_cpu_wrap_dslaxpby: not implemented for alpha = %lf\n", alpha);
	//fprintf(stderr,"custom_cpu_wrap_dslaxpby: using %d openmp workers\n", omp_get_num_threads());
	#pragma omp parallel for proc_bind(close) num_threads(HW_THREADS/(NUMA_HW_NUM))
	for (i = 0; i < offset_x; i++){
		for (j = 0; j < N; j++){
		//fprintf(stderr, "y[%d] = ax[%d] + by[%d] , a = %lf, b = %lf\n", i*ptr_ker_translate->slide_y + j, i*ptr_ker_translate->N + j, 
		//  i*ptr_ker_translate->slide_y + j, ptr_ker_translate->alpha, ptr_ker_translate->beta);
			y[i*offset_y + j] = x[i*N + j] + beta*y[i*offset_y + j];
		}
	}
  return NULL;
}

void custom_cpu_wrap_dslaxpby_v02(void* backend_data){
  slaxpby_backend_in<double>* ptr_ker_translate = (slaxpby_backend_in<double>*) backend_data;
  pthread_t core_binder; 
	if( pthread_create(&core_binder, NULL, custom_cpu_wrap_dslaxpby_pthread_wrap, backend_data) != 0)
		error("pthread_create failed for custom_cpu_wrap_dslaxpby\n");
  if(pthread_join(core_binder, NULL) != 0 ) error("pthread_join failed for custom_cpu_wrap_dslaxpby\n");
}

//#pragma omp declare simd
void inline simp_d(double x, double beta, double* y){
  *y =  x + beta*(*y);
}

void custom_cpu_wrap_dslaxpby(void* backend_data){
  slaxpby_backend_in<double>* ptr_ker_translate = (slaxpby_backend_in<double>*) backend_data;
#ifdef DEBUG
  fprintf(stderr,"custom_cpu_wrap_dslaxpby(dev_id = %d,\
    N = %d, alpha = %lf, x = %p, incx = %d, b = %lf, y = %p, incy = %d, slide_x = %d, slide_y = %d)\n",
    ptr_ker_translate->dev_id, ptr_ker_translate->N, ptr_ker_translate->alpha,
    (double*) *ptr_ker_translate->x, ptr_ker_translate->incx, ptr_ker_translate->beta,
    (double*) *ptr_ker_translate->y, ptr_ker_translate->incy, ptr_ker_translate->slide_x, ptr_ker_translate->slide_y);
#endif
  double* y = (double*) *ptr_ker_translate->y, *x = (double*) *ptr_ker_translate->x;
	int i, j, N = ptr_ker_translate->N, offset_x = ptr_ker_translate->slide_x, offset_y = ptr_ker_translate->slide_y; 
	double alpha = ptr_ker_translate->alpha, beta = ptr_ker_translate->beta;
  if(ptr_ker_translate->alpha != 1.0) error("custom_avx2_cpu_wrap_dslaxpby: not implemented for alpha = %lf\n", alpha);
	//fprintf(stderr,"custom_cpu_wrap_dslaxpby: using %d openmp workers\n", omp_get_num_threads());
  int numa_memlocs = CHL_MEMLOCS-1-CHL_WORKERS;
	long long thread_offset_chunks = offset_x/(HW_THREADS/numa_memlocs) + ((offset_x%(HW_THREADS/numa_memlocs) ? 1 : 0));
	#pragma omp parallel for schedule(static,1)
	for(int thread_id = 0; thread_id < HW_THREADS; thread_id++){
		int thread_idx = omp_get_thread_num(); // This *should* be the same with thread_id for OMP_PROC_BIND=TRUE
		int thread_core = thread_idx/(HW_THREADS/numa_memlocs);
		int thread_ctr = thread_idx%(HW_THREADS/numa_memlocs);
    if (thread_core == ptr_ker_translate->dev_id - CHL_WORKERS){
      for (i = thread_ctr*thread_offset_chunks; i < (thread_ctr + 1)*thread_offset_chunks; i++){
        if(i < offset_x) for (j = 0; j < N; j++){
          //fprintf(stderr, "y[%d] = ax[%d] + by[%d] , a = %lf, b = %lf\n", i*ptr_ker_translate->slide_y + j, i*ptr_ker_translate->N + j, 
          //  i*ptr_ker_translate->slide_y + j, ptr_ker_translate->alpha, ptr_ker_translate->beta);
          simp_d(x[i*N + j], beta, &(y[i*offset_y + j]));
        }
      }
    }
  }
}

#include <immintrin.h>

void custom_avx2_cpu_wrap_dslaxpby(void* backend_data){
  slaxpby_backend_in<double>* ptr_ker_translate = (slaxpby_backend_in<double>*) backend_data;
#ifdef DEBUG
  fprintf(stderr,"custom_cpu_wrap_dslaxpby(dev_id = %d,\
    N = %d, alpha = %lf, x = %p, incx = %d, b = %lf, y = %p, incy = %d, slide_x = %d, slide_y = %d)\n",
    ptr_ker_translate->dev_id, ptr_ker_translate->N, ptr_ker_translate->alpha,
    (double*) *ptr_ker_translate->x, ptr_ker_translate->incx, ptr_ker_translate->beta,
    (double*) *ptr_ker_translate->y, ptr_ker_translate->incy, ptr_ker_translate->slide_x, ptr_ker_translate->slide_y);
#endif
  double* y = (double*) *ptr_ker_translate->y, *x = (double*) *ptr_ker_translate->x;
	int i, j, N = ptr_ker_translate->N, offset_x = ptr_ker_translate->slide_x, offset_y = ptr_ker_translate->slide_y; 
	double alpha = ptr_ker_translate->alpha, beta = ptr_ker_translate->beta;
  if(ptr_ker_translate->alpha != 1.0) error("custom_avx2_cpu_wrap_dslaxpby: not implemented for alpha = %lf\n", alpha);
	//fprintf(stderr,"custom_cpu_wrap_dslaxpby: using %d openmp workers\n", omp_get_num_threads());
  int numa_memlocs = CHL_MEMLOCS-1-CHL_WORKERS;
	long long thread_offset_chunks = offset_x/(HW_THREADS/numa_memlocs) + ((offset_x%(HW_THREADS/numa_memlocs) ? 1 : 0));
	#pragma omp parallel for schedule(static,1)
	for(int thread_id = 0; thread_id < HW_THREADS; thread_id++){
		int thread_idx = omp_get_thread_num(); // This *should* be the same with thread_id for OMP_PROC_BIND=TRUE
		if(thread_idx!=thread_id) warning("custom_cpu_wrap_dslaxpby: Will not work as intended"
    "since OMP threads are not pinned to HW correctly (?)\n");
		int thread_core = thread_idx/(HW_THREADS/numa_memlocs);
		int thread_ctr = thread_idx%(HW_THREADS/numa_memlocs);
    if (thread_core == ptr_ker_translate->dev_id - CHL_WORKERS){
      for (i = thread_ctr*thread_offset_chunks; i < (thread_ctr + 1)*thread_offset_chunks; i++){
        if(i < offset_x) for (j = 0; j < N; j+=4){
          __m256d v_y, v_x, v_beta;
          //fprintf(stderr, "y[%d] = ax[%d] + by[%d] , a = %lf, b = %lf\n", i*ptr_ker_translate->slide_y + j, i*ptr_ker_translate->N + j, 
          //  i*ptr_ker_translate->slide_y + j, ptr_ker_translate->alpha, ptr_ker_translate->beta);
          v_y =_mm256_load_pd(&y[i*offset_y + j]);
          v_x =_mm256_load_pd(&x[i*N + j]);
	        v_y =_mm256_add_pd(v_y,v_x);
	        _mm256_store_pd(&y[i*offset_y + j], v_y);
          //v_y = _mm256_set_pd(beta*y[i*offset_y + j], beta*y[i*offset_y + j + 1], 
          //  beta*y[i*offset_y + j + 2], beta*y[i*offset_y + j + 3]);
			    //v_sum = _mm256_fmadd_pd(v_a, v_x, v_sum);
          //y[i*offset_y + j] = alpha*x[i*N + j] + beta*y[i*offset_y + j];
        }
      }
    }
  }
}


void cblas_wrap_dgemm(void* backend_data){
  gemm_backend_in<double>* ptr_ker_translate = (gemm_backend_in<double>*) backend_data;
#ifdef DDEBUG
  fprintf(stderr,"cblas_wrap_dgemm: cblas_dgemm(dev_id = %d, TransA = %c, TransB = %c,\
    M = %d, N = %d, K = %d, alpha = %lf, A = %p, lda = %d, \n\
    B = %p, ldb = %d, beta = %lf, C = %p, ldC = %d)\n",
    ptr_ker_translate->dev_id, ptr_ker_translate->TransA, ptr_ker_translate->TransB,
    ptr_ker_translate->M, ptr_ker_translate->N, ptr_ker_translate->K, ptr_ker_translate->alpha,
    (double*) *ptr_ker_translate->A, ptr_ker_translate->ldA,
    (double*) *ptr_ker_translate->B, ptr_ker_translate->ldB,
    ptr_ker_translate->beta, (double*) *ptr_ker_translate->C, ptr_ker_translate->ldC);
#endif
  cblas_dgemm(CblasColMajor,
    OpCharToCblas(ptr_ker_translate->TransA), OpCharToCblas(ptr_ker_translate->TransB),
    ptr_ker_translate->M, ptr_ker_translate->N, ptr_ker_translate->K, ptr_ker_translate->alpha,
    (double*) *ptr_ker_translate->A, ptr_ker_translate->ldA,
    (double*) *ptr_ker_translate->B, ptr_ker_translate->ldB,
    ptr_ker_translate->beta, (double*) *ptr_ker_translate->C, ptr_ker_translate->ldC);
}

void cublas_wrap_daxpy(void* backend_data, void* queue_wrap_p){
  axpy_backend_in<double>* ptr_ker_translate = (axpy_backend_in<double>*) backend_data;
  //CHLSelectDevice(ptr_ker_translate->dev_id);
#ifdef DEBUG
  fprintf(stderr,"cublas_wrap_daxpy: cublasDaxpy(dev_id = %d,\
    N = %d, alpha = %lf, x = %p, incx = %d, y = %p, incy = %d)\n",
    ptr_ker_translate->dev_id, ptr_ker_translate->N, ptr_ker_translate->alpha,
    (double*) *ptr_ker_translate->x, ptr_ker_translate->incx,
    (double*) *ptr_ker_translate->y, ptr_ker_translate->incy);
#endif
  cublasHandle_t temp_handle = *((cublasHandle_t*)((CQueue_p)queue_wrap_p)->backend_comp_md);

  cublasStatus_t stat = cublasDaxpy(temp_handle,
    ptr_ker_translate->N, (double*) &ptr_ker_translate->alpha, (double*) *ptr_ker_translate->x,
    ptr_ker_translate->incx, (double*) *ptr_ker_translate->y, ptr_ker_translate->incy);
  if(stat != CUBLAS_STATUS_SUCCESS) error("cublas_wrap_daxpy failed: cublasDaxpy(temp_handle = %p, dev_id = %d,\
    N = %d, alpha = %lf, x = %p, incx = %d, y = %p, incy = %d)\n",
    temp_handle, ptr_ker_translate->dev_id, ptr_ker_translate->N, ptr_ker_translate->alpha,
    (double*) *ptr_ker_translate->x, ptr_ker_translate->incx,
    (double*) *ptr_ker_translate->y, ptr_ker_translate->incy);
}

void cublas_wrap_daxpby(void* backend_data, void* queue_wrap_p){
  axpby_backend_in<double>* ptr_ker_translate = (axpby_backend_in<double>*) backend_data;
  //CHLSelectDevice(ptr_ker_translate->dev_id);
#ifdef DEBUG
  fprintf(stderr,"cublas_wrap_daxpby:\n cublasDscal(dev_id = %d,\
    N = %d, alpha = %lf, x = %p, incx = %d)\ncublasDaxpy(dev_id = %d,\
    N = %d, alpha = %lf, x = %p, incx = %d, y = %p, incy = %d)\n",
    ptr_ker_translate->dev_id, ptr_ker_translate->N, ptr_ker_translate->beta, (double*) *ptr_ker_translate->y, ptr_ker_translate->incy,
    ptr_ker_translate->dev_id, ptr_ker_translate->N, ptr_ker_translate->alpha, 
    (double*) *ptr_ker_translate->x, ptr_ker_translate->incx,
    (double*) *ptr_ker_translate->y, ptr_ker_translate->incy);
#endif
    cublasHandle_t temp_handle = *((cublasHandle_t*)((CQueue_p)queue_wrap_p)->backend_comp_md);

  massert(CUBLAS_STATUS_SUCCESS == cublasDscal(temp_handle,
    ptr_ker_translate->N, (double*) &ptr_ker_translate->beta, (double*) *ptr_ker_translate->y,
    ptr_ker_translate->incy), "cublasDscal failed\n");
  massert(CUBLAS_STATUS_SUCCESS == cublasDaxpy(temp_handle,
    ptr_ker_translate->N, (double*) &ptr_ker_translate->alpha, (double*) *ptr_ker_translate->x,
    ptr_ker_translate->incx, (double*) *ptr_ker_translate->y, ptr_ker_translate->incy),
    "cublasDaxpy failed\n");
}

void cublas_wrap_dgemm(void* backend_data, void* queue_wrap_p){
  gemm_backend_in<double>* ptr_ker_translate = (gemm_backend_in<double>*) backend_data;
  //CHLSelectDevice(ptr_ker_translate->dev_id);
  cublasHandle_t temp_handle = *((cublasHandle_t*)((CQueue_p)queue_wrap_p)->backend_comp_md);
#ifdef DEBUG
  fprintf(stderr,"cublas_wrap_dgemm: cublasDgemm(temp_handle = %p, dev_id = %d, TransA = %c, TransB = %c,\
    M = %d, N = %d, K = %d, alpha = %lf, A = %p, lda = %d, \n\
    B = %p, ldb = %d, beta = %lf, C = %p, ldC = %d)\n",
    temp_handle, ptr_ker_translate->dev_id, ptr_ker_translate->TransA, ptr_ker_translate->TransB,
    ptr_ker_translate->M, ptr_ker_translate->N, ptr_ker_translate->K, ptr_ker_translate->alpha,
    (double*) *ptr_ker_translate->A, ptr_ker_translate->ldA,
    (double*) *ptr_ker_translate->B, ptr_ker_translate->ldB,
    ptr_ker_translate->beta, (double*) *ptr_ker_translate->C, ptr_ker_translate->ldC);
#endif

  massert(CUBLAS_STATUS_SUCCESS == cublasDgemm(temp_handle,
    OpCharToCublas(ptr_ker_translate->TransA), OpCharToCublas(ptr_ker_translate->TransB),
    ptr_ker_translate->M, ptr_ker_translate->N, ptr_ker_translate->K, &ptr_ker_translate->alpha,
    (double*) *ptr_ker_translate->A, ptr_ker_translate->ldA,
    (double*) *ptr_ker_translate->B, ptr_ker_translate->ldB,
    &ptr_ker_translate->beta, (double*) *ptr_ker_translate->C, ptr_ker_translate->ldC),
    "cublas_wrap_dgemm: cublasDgemm failed\n");
}