///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief Some CUDA function calls with added error-checking
///

#include <cstdio>
#include <typeinfo>
#include <float.h>
#include <curand.h>
#include <cuda.h>

#include "chl_smart_wrappers.hpp"
#include "backend_wrappers.hpp"

int CHLGetDevice(){
  int dev_id = -1;
  cudaError_t err = cudaGetDevice(&dev_id);
  massert(cudaSuccess == err,
    "CHLGetDevice: cudaGetDevice failed - %s\n", cudaGetErrorString(err));
  return dev_id;
}

void CHLSelectDevice(short dev_id){
  if(dev_id >= 0 && dev_id < CHL_WORKERS){
  cudaError_t err = cudaSetDevice(dev_id);
  massert(cudaSuccess == err,
    "CHLSelectDevice(%d): cudaSetDevice(%d) failed - %s\n", dev_id, dev_id, cudaGetErrorString(err));
  }
  else if(dev_id >= CHL_WORKERS && dev_id < CHL_MEMLOCS){  /// "Host" device loc id used by CHL
    cudaSetDevice((dev_id - CHL_WORKERS)*CHL_WORKERS/(CHL_MEMLOCS-CHL_WORKERS));
  }
  else error("CHLSelectDevice(%d): invalid dev_id\n", dev_id);
}

void CHLDevGetMemInfo(long long* free_dev_mem, long long* max_dev_mem){
  size_t free_dev_mem_tmp, max_dev_mem_tmp;
    int tmp_dev_id;
    cudaError_t err = cudaGetDevice(&tmp_dev_id);
    // TODO: For the CPU this function returns device 0 memory availability. Its a feature not a bug.
    massert(cudaSuccess == err,
      "CHLDevGetMemInfo: cudaGetDevice failed - %s\n", cudaGetErrorString(err));
    err = cudaMemGetInfo(&free_dev_mem_tmp, &max_dev_mem_tmp);
  	massert(cudaSuccess == err,
      "CHLDevGetMemInfo: cudaMemGetInfo failed - %s\n", cudaGetErrorString(err));
    *free_dev_mem = (long long) free_dev_mem_tmp;
    *max_dev_mem = (long long) max_dev_mem_tmp;
}

void CHLSyncCheckErr(){
  int prev_loc = CHLGetDevice();
  for(int dev_idx = 0; dev_idx < CHL_WORKERS; dev_idx++)
  {
    CHLSelectDevice(dev_idx); 
    cudaError_t errSync = cudaDeviceSynchronize();
    if (errSync != cudaSuccess)
      printf("Sync kernel error: %s\n", cudaGetErrorString(errSync));
  }
  CHLSelectDevice(prev_loc); 
}

void CHLASyncCheckErr(){
  cudaError_t errAsync = cudaGetLastError();
  if (errAsync != cudaSuccess)
    printf("Async kernel error: %s\n", cudaGetErrorString(errAsync));
}

void cudaCheckErrors(){
	//CHLASyncCheckErr();
	CHLSyncCheckErr();
}


void CHLEnableLinks(short target_dev_i, short num_devices){
	short lvl = 2;
#ifdef DEBUG
	lprintf(lvl-1, "|-----> CHLEnableGPUPeer(%d,%d)\n", target_dev_i, num_devices);
#endif
#ifdef TEST
	lprintf(lvl-1, "|-----> CHLEnableGPUPeer\n");
	double cpu_timer = csecond();
#endif
	int dev_id_target = target_dev_i;
	CHLSelectDevice(dev_id_target);
	for(int j=0; j<num_devices;j++){
		int dev_id_current = j;
		if (dev_id_target == dev_id_current || dev_id_target >= CHL_WORKERS || dev_id_current >= CHL_WORKERS) continue;
		int can_access_peer;
		massert(cudaSuccess == cudaDeviceCanAccessPeer(&can_access_peer, dev_id_target, dev_id_current), "PARALiaDgemm: cudaDeviceCanAccessPeer failed\n");
		if(can_access_peer){
			cudaError_t check_peer = cudaDeviceEnablePeerAccess(dev_id_current, 0);
			if(check_peer == cudaSuccess){ ;
#ifdef DEBUG
				lprintf(lvl, "Enabled Peer access for dev %d to dev %d\n", dev_id_target, dev_id_current);
#endif
			}
			else if (check_peer == cudaErrorPeerAccessAlreadyEnabled){
				cudaGetLastError();
#ifdef DEBUG
				lprintf(lvl, "Peer access already enabled for dev %d to dev %d\n", dev_id_target, dev_id_current);
#endif
			}
			else error("Enabling Peer access failed for %d to dev %d\n", dev_id_target, dev_id_current);
		}
	}
#ifdef TEST
	cpu_timer = csecond() - cpu_timer;
	lprintf(lvl, "Utiilizing Peer access for dev %d -> t_enable =%lf ms\n", dev_id_target, 1000*cpu_timer);
	cpu_timer = csecond();
	lprintf(lvl-1, "<-----|\n");
#endif
#ifdef DEBUG
	lprintf(lvl-1, "<-----|\n");
#endif
}

__global__ void cuda_hello(){
    printf("\n\n<<<<<<<<<<<<<<<<<<<Hello World from GPU!>>>>>>>>>>>>>>>>>>>>>\n\n");
}

__global__ void dslaxpby(long long N, double alpha, double* x, long long offset_x, double beta, double* y, long long offset_y)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  int j = blockIdx.y*blockDim.y + threadIdx.y;
  if (i >= offset_x || j >= N) return;
  y[i*offset_y + j] = /*alpha**/x[i*offset_x + j] + beta * y[i*offset_y + j];
  //printf("Here %d, %d\n", i, j);
}

__global__ void dslaxpby_gridstride(long long N, double alpha, 
  double* x,  long long offset_x, double beta, double* y, long long offset_y)
{
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < offset_x; i += blockDim.x * gridDim.x) 
    {
      for (int j = blockIdx.y * blockDim.y + threadIdx.y; j < N; j += blockDim.y * gridDim.y) 
      {
          y[i*offset_y + j] = /*alpha**/x[i*offset_x + j] + beta * y[i*offset_y + j];
      }
    }
}

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

void custom_gpu_wrap_dslaxpby(void* backend_data, CQueue_p run_queue){
  slaxpby_backend_in* ptr_ker_translate = (slaxpby_backend_in*) backend_data;
#ifdef DEBUG
  fprintf(stderr,"custom_gpu_wrap_dslaxpby(dev_id = %d,\
    N = %d, alpha = %lf, x = %p, incx = %d, b = %lf, y = %p, incy = %d, slide_x = %d, slide_y = %d)\n",
    ptr_ker_translate->dev_id, ptr_ker_translate->N, *((double*)ptr_ker_translate->alpha),
    (double*) *ptr_ker_translate->x, ptr_ker_translate->incx, *((double*)ptr_ker_translate->beta),
    (double*) *ptr_ker_translate->y, ptr_ker_translate->incy, ptr_ker_translate->slide_x, ptr_ker_translate->slide_y);
#endif
  /*
  dim3 grid_sz((ptr_ker_translate->slide_x), (ptr_ker_translate->N + 1023)/1024);
  dim3 block_sz(1, 1024);
  dslaxpby<<<grid_sz, block_sz, 0, *((cudaStream_t*)run_queue->backend_queue_ptr) >>>
  (ptr_ker_translate->N, *((double*)ptr_ker_translate->alpha), (double*) *ptr_ker_translate->x, ptr_ker_translate->slide_x, 
  ptr_ker_translate->beta, (double*) *ptr_ker_translate->y, ptr_ker_translate->slide_y);
  */
  int numSMs;
  cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, ptr_ker_translate->dev_id);
  dim3 grid_sz(numSMs, 1), block_sz(1, 1024);
  cudaStream_t stream = *((cudaStream_t*) run_queue->backend_queue_ptr);
  dslaxpby_gridstride<<<grid_sz, block_sz, 0, stream>>>
  (ptr_ker_translate->N, *((double*)ptr_ker_translate->alpha), (double*) *ptr_ker_translate->x, ptr_ker_translate->slide_x, 
  *((double*)ptr_ker_translate->beta), (double*) *ptr_ker_translate->y, ptr_ker_translate->slide_y);
  
  //gpuErrchk( cudaPeekAtLastError() );
  //gpuErrchk( cudaDeviceSynchronize() );

}