///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief Some CUDA function calls with added error-checking
///

#include <cuda.h>
#include <cublas_v2.h>
#include <curand.h>

#include "sys/types.h"
#include "sys/sysinfo.h"
#include <numa.h>

#include "smart_wrappers.hpp"

void CHLGetMaxCPUmem(size_t* free_mem, size_t* max_mem){
	struct sysinfo memInfo;
	sysinfo (&memInfo);
	long long virtualMemFree = memInfo.freeram;
	//Add other values in next statement to avoid int overflow on right hand side...
	virtualMemFree += memInfo.freeswap;
	virtualMemFree *= memInfo.mem_unit;
	*free_mem = (size_t) virtualMemFree;
	long long totalVirtualMem = memInfo.totalram;
	//Add other values in next statement to avoid int overflow on right hand side...
	totalVirtualMem += memInfo.totalswap;
	totalVirtualMem *= memInfo.mem_unit;
	*max_mem = (size_t) totalVirtualMem;
	return;
}

long int CHLGetMaxDimSqAsset2D(short Asset2DNum, short dsize, long int step, int loc){
	size_t free_mem, max_mem;
	if (loc >= 0){
		int prev_loc; cudaGetDevice(&prev_loc);
		cudaSetDevice(loc);
		massert(cudaSuccess == cudaMemGetInfo(&free_mem, &max_mem), "backend_get_max_dim_sq_Asset2D: cudaMemGetInfo failed");
		cudaSetDevice(prev_loc);
	} else CHLGetMaxCPUmem(&free_mem, &max_mem);

	// Define the max size of a benchmark kernel to run on this machine.
	long int maxDim = (( (long int) sqrt((free_mem*PROBLEM_GPU_PERCENTAGE/100.0)/(Asset2DNum*dsize))) / step) * step;
	return maxDim;
}

long int CHLGetMaxDimAsset1D(short Asset1DNum, short dsize, long int step, int loc){
	size_t free_mem, max_mem;
	if (loc >= 0){
		int prev_loc; cudaGetDevice(&prev_loc);
		cudaSetDevice(loc);
		massert(cudaSuccess == cudaMemGetInfo(&free_mem, &max_mem), "backend_get_max_dim_sq_Asset2D: cudaMemGetInfo failed");
		cudaSetDevice(prev_loc);
	} else CHLGetMaxCPUmem(&free_mem, &max_mem);

	long int maxDim = (( (long int) (free_mem*PROBLEM_GPU_PERCENTAGE/100.0)/(Asset1DNum*dsize)) / step) * step;
	return maxDim;
}

void *gpu_malloc(long long count) {
	void *ret;
	massert(cudaMalloc(&ret, count) == cudaSuccess,
		cudaGetErrorString(cudaGetLastError()));
	//long long chunk_count = 1;
	//int* slice_memloc_in = (int*) malloc (chunk_count*sizeof(int));
	//slice_memloc_in[0] = CHLGetDevice(); 
	//allocated_data[allocated_data_num++] = new MEMMetadata(ret, slice_memloc_in, count, 1);
	return ret;
}

void pin_mem_wrap(void** ptr, long long bytes){
	cudaHostRegister(*ptr,bytes,cudaHostRegisterPortable);
}

void *pin_malloc(long long count, short W_flag) {
	void *ret;
	massert(cudaHostAlloc ( &ret, count, cudaHostAllocDefault)  == cudaSuccess,
	cudaGetErrorString(cudaGetLastError()));
	return ret;
}

void *numa_inter_pin_malloc(long long count, short W_flag){
	void *ret;
	ret = numa_alloc_interleaved(count);//numa_alloc_onnode(count, 0);//

	//long long chunk_count = 1;
	//int* slice_memloc_in = (int*) malloc (chunk_count*sizeof(int));
	//slice_memloc_in[0] = CHL_MEMLOCS - 1; 
	//allocated_data[allocated_data_num++] = new MEMMetadata(ret, slice_memloc_in, count, 1);
	//if (!W_flag) cudaHostRegister(ret,count,cudaHostRegisterPortable | cudaHostRegisterReadOnly);
	//else 
	cudaHostRegister(ret,count,cudaHostRegisterPortable);
	return ret;
}

void *numa_bind_pin_malloc(long long count, int node_num, short W_flag){
	void *ret;
	ret = numa_alloc_onnode(count, translate_mem_idx_to_hw(node_num));
	//long long chunk_count = 1;
	//int* slice_memloc_in = (int*) malloc (chunk_count*sizeof(int));
	//slice_memloc_in[0] = node_num + CHL_WORKERS; 
	//allocated_data[allocated_data_num++] = new MEMMetadata(ret, slice_memloc_in, count, 1);
	//if (!W_flag) cudaHostRegister(ret,count,cudaHostRegisterPortable | cudaHostRegisterReadOnly);
	//else 
	cudaHostRegister(ret,count,cudaHostRegisterPortable);
	return ret;
}

void* CHLMalloc(long long bytes, int loc, short W_flag){
  void *ptr = NULL;
  if (loc == CHL_MEMLOCS - 1) {
#ifdef CLDEBUG
    fprintf(stderr, "Allocating %lld bytes to interleaved NUMA alloc...\n", bytes);
#endif
	ptr = numa_inter_pin_malloc(bytes, W_flag);
  }
  else if (loc >= CHL_WORKERS && loc < CHL_MEMLOCS - 1) {
#ifdef CLDEBUG
    fprintf(stderr, "Allocating %lld bytes to NUM MEM %d [%s]...\n", bytes, loc, mem_name(loc));
#endif
	ptr = numa_bind_pin_malloc(bytes, loc, W_flag);
  }
  else if (loc == -1){
#ifdef CLDEBUG
    fprintf(stderr, "Allocating %lld bytes and touching in close numa nodes with CHLMallocHostTouchSmart\n", bytes);
#endif
	ptr = CHLMallocHostTouchSerial(bytes);
		//CHLMallocHostTouchSmart(std::sqrt(bytes/sizeof(double)),
		 //std::sqrt(bytes/sizeof(double)), sizeof(double), 'N');
  }
  else if (loc >= 0 && loc < CHL_WORKERS){
    int prev_loc; cudaGetDevice(&prev_loc);
#ifdef CLDEBUG
    fprintf(stderr, "Allocating %lld bytes to device(%d)...\n", bytes, loc);
    //if (prev_loc != loc) warning("CHLMalloc: Malloc'ed memory in other device (Previous device: %d, Malloc in: %d)\n", prev_loc, loc);
#endif
    cudaSetDevice(loc);
    ptr = gpu_malloc(bytes);


	CHLSyncCheckErr();
    if (prev_loc != loc)cudaSetDevice(prev_loc);
  }
  else error("CHLMalloc: Invalid device id/location %d\n", loc);
  CHLSyncCheckErr();
  return ptr;
}

void gpu_free(void *gpuptr) {
  massert(cudaFree(gpuptr) == cudaSuccess,
          cudaGetErrorString(cudaGetLastError()));
}

void pin_free(void *gpuptr) {
  massert(cudaFreeHost(gpuptr) == cudaSuccess,
          cudaGetErrorString(cudaGetLastError()));
}

void numa_pin_free(void *gpuptr, long long bytes) {
  cudaHostUnregister(gpuptr);
  numa_free(gpuptr,bytes);
}


void CHLFree(void * ptr, long long bytes, int loc){
	//if (??? == loc) free(ptr);
	if ((loc >= CHL_WORKERS && loc < CHL_MEMLOCS) || loc == -1) numa_pin_free(ptr, bytes);
	else if (loc >= 0 && loc < CHL_WORKERS){
		int prev_loc; cudaGetDevice(&prev_loc);
		//if (prev_loc != loc) warning("CHLFree: Freed memory in other device (Previous device: %d, Free in: %d)\n", prev_loc, loc);
		cudaSetDevice(loc);
		gpu_free(ptr);
		CHLSyncCheckErr();
		if (prev_loc != loc){
		//warning("CHLFree: Reseting device to previous: %d\n", prev_loc);
			cudaSetDevice(prev_loc);
		}
	}
	else error("CHLFree: Invalid device id/location %d\n", loc);
	CHLSyncCheckErr();
}


short CHLGetPtrLoc(void * in_ptr)
{
// This is legacy code for CUDA 9.2 <<. It should not be used due to CUDA ptr_att back-end struct changes in latest versions
#ifdef CUDA_9_WRAPPER_MESS
  error("CHLGetPtrLoc(9.2 version, ptr =%p): Not supported anymore\n", in_ptr);	
#else
	int loc = -42;
	cudaPointerAttributes ptr_att;
	if (cudaSuccess != cudaPointerGetAttributes(&ptr_att, in_ptr)) error("CHLGetPtrLoc(cuda 10+ version, ptr =%p):\
	//Pointer not visible to CUDA, host alloc or error\n", in_ptr);
	if (ptr_att.type == cudaMemoryTypeHost){
		loc = translate_hw_numa_to_mem_idx(get_hw_numa_idx(in_ptr));
		if (loc == -1) error("CHLGetPtrLoc(cuda 10+ version, ptr =%p): scan_allocated_data_for_ptr_loc(ptr) did not find pointer\n", in_ptr);
	}
	else if (ptr_att.type == cudaMemoryTypeDevice) loc = ptr_att.device;
	// TODO: Unified memory is considered available in the GPU as cuBLASXt ( not bad, not great)
	else if (ptr_att.type == cudaMemoryTypeManaged) loc = ptr_att.device;
	else error("CHLGetPtrLoc(cuda 10+ version, loc = %d, ptr =%p): Invalid memory type\n", ptr_att.device, in_ptr);
	return loc;
#endif
}

short CHLGetPtrAdvLoc(void * in_ptr, long long dim1, long long dim2, int elemSize)
{
	int loc = -42;
	cudaPointerAttributes ptr_att;
	if (cudaSuccess != cudaPointerGetAttributes(&ptr_att, in_ptr)) error("CHLGetPtrAdvLoc(cuda 10+ version, ptr =%p):"
	"Pointer not visible to CUDA, host alloc or error\n", in_ptr);
	if (ptr_att.type == cudaMemoryTypeHost){
		loc = translate_hw_numa_to_mem_idx(get_hw_numa_idx(in_ptr));
		if (loc == -1) error("CHLGetPtrLoc(cuda 10+ version, ptr =%p): scan_allocated_data_for_ptr_loc(ptr)"
		"did not find pointer\n", in_ptr);
		int advanced_allocation = 0;
		//Check for custom chunk-interleaved matrix at host memlocs.
		if(loc == CHL_WORKERS){
			advanced_allocation = 1; 
			long long addr_chunk_sz = dim1*dim2*elemSize / (CHL_MEMLOCS - CHL_WORKERS - 1);
			for (int idx = 1; idx < CHL_MEMLOCS - CHL_WORKERS - 1; idx++) if (CHL_WORKERS + idx != 
				translate_hw_numa_to_mem_idx(get_hw_numa_idx(in_ptr +idx*addr_chunk_sz))) advanced_allocation = 0; 
		}
		if (advanced_allocation){
			fprintf(stderr, "CHLGetPtrAdvLoc(%p, %lld, %lld, %d):"
			"Advanced allocation identified (type = chunk-interleaved matrix)\n", in_ptr, dim1, dim2, elemSize); 
			return -1;
		}
		else return loc; 
	}
	else return CHLGetPtrLoc(in_ptr);
}

void CHLMemcpy(void* dest, void* src, long long bytes, int loc_dest, int loc_src)
{
	massert(loc_dest >= -1 && loc_dest < CHL_MEMLOCS, "CHLMemcpy: Invalid destination device: %d\n", loc_dest);
	massert(loc_src >= -1 && loc_src < CHL_MEMLOCS, "CHLMemcpy: Invalid source device: %d\n", loc_src);

	enum cudaMemcpyKind kind;
	if ((loc_src >= CHL_WORKERS  || loc_src == -1 ) && (loc_dest >= CHL_WORKERS   || loc_dest == -1 )) 
		kind = cudaMemcpyHostToHost;
	else if (loc_dest >= CHL_WORKERS  || loc_dest == -1) kind = cudaMemcpyDeviceToHost;
	else if (loc_src >= CHL_WORKERS || loc_src == -1) kind = cudaMemcpyHostToDevice;
	else kind = cudaMemcpyDeviceToDevice;

#ifdef DEBUG
	if (loc_src == loc_dest) warning("CHLMemcpy(dest=%p, src=%p, bytes=%lld, loc_dest=%d, loc_src=%d): Source location matches destination\n",
	dest, src, bytes, loc_dest, loc_src);
#endif
	massert(CUBLAS_STATUS_SUCCESS == cudaMemcpy(dest, src, bytes, kind), "CHLMemcpy: cudaMemcpy from device src=%d to dest=%d failed\n", loc_src, loc_dest);
	CHLSyncCheckErr();
}

/*void CHLMemcpyAsync(void* dest, void* src, long long bytes, int loc_dest, int loc_src, CQueue_p transfer_queue)
{

	cudaStream_t stream = *((cudaStream_t*) transfer_queue->backend_queue_ptr);

	massert(loc_dest >= -1 && loc_dest < CHL_MEMLOCS, "CHLMemcpyAsync: Invalid destination device: %d\n", loc_dest);
	massert(loc_src >= -1 && loc_src < CHL_MEMLOCS, "CHLMemcpyAsync: Invalid source device: %d\n", loc_src);

	enum cudaMemcpyKind kind;
	if ((loc_src >= CHL_WORKERS  || loc_src == -1 ) && (loc_dest >= CHL_WORKERS   || loc_dest == -1 )) 
		kind = cudaMemcpyHostToHost;
	else if (loc_dest >= CHL_WORKERS  || loc_dest == -1) kind = cudaMemcpyDeviceToHost;
	else if (loc_src >= CHL_WORKERS || loc_src == -1) kind = cudaMemcpyHostToDevice;
	else kind = cudaMemcpyDeviceToDevice;

	if (loc_src == loc_dest) warning("CHLMemcpyAsync(dest=%p, src=%p, bytes=%lld, loc_dest=%d, loc_src=%d): Source location matches destination\n",
	dest, src, bytes, loc_dest, loc_src);
	massert(cudaSuccess == cudaMemcpyAsync(dest, src, bytes, kind, stream),
	"CHLMemcpy2D: cudaMemcpyAsync failed\n");
	//CHLSyncCheckErr();
}*/

void CHLMemcpy2D(void* dest, long int ldest, void* src, long int ldsrc, long int rows, long int cols, short elemSize, int loc_dest, int loc_src){
	short lvl = 6;
#ifdef DDEBUG
	lprintf(lvl, "CHLMemcpy2D(dest=%p, ldest =%zu, src=%p, ldsrc = %zu, rows = %zu, cols = %zu, elemsize = %d, loc_dest = %d, loc_src = %d)\n",
		dest, ldest, src, ldsrc, rows, cols, elemSize, loc_dest, loc_src);
#endif
	massert(loc_dest >= -1 && loc_dest < CHL_MEMLOCS, "CHLMemcpy2D: Invalid destination device: %d\n", loc_dest);
	massert(loc_src >= -1 && loc_src < CHL_MEMLOCS, "CHLMemcpy2D: Invalid source device: %d\n", loc_src);

	enum cudaMemcpyKind kind;
	if ((loc_src >= CHL_WORKERS  || loc_src == -1 ) && (loc_dest >= CHL_WORKERS   || loc_dest == -1 )) 
		kind = cudaMemcpyHostToHost;
	else if (loc_dest >= CHL_WORKERS  || loc_dest == -1) kind = cudaMemcpyDeviceToHost;
	else if (loc_src >= CHL_WORKERS || loc_src == -1) kind = cudaMemcpyHostToDevice;
	else kind = cudaMemcpyDeviceToDevice;

	if (loc_src == loc_dest) warning("CHLMemcpy2D(dest=%p, ldest =%zu, src=%p, ldsrc = %zu, rows=%zu, cols=%zu, elemSize =%d, loc_dest=%d, loc_src=%d): Source location matches destination\n",
	dest, ldest, src, ldsrc, rows, cols, elemSize, loc_dest, loc_src);
	//if (loc_src == -1 && loc_dest >=0) massert(CUBLAS_STATUS_SUCCESS == cublasSetMatrix(rows, cols, elemSize, src, ldsrc, dest, ldest), "CHLMemcpy2DAsync: cublasSetMatrix failed\n");
	//else if (loc_src >=0 && loc_dest == -1) massert(CUBLAS_STATUS_SUCCESS == cublasGetMatrix(rows, cols, elemSize, src, ldsrc, dest, ldest),  "CHLMemcpy2DAsync: cublasGetMatrix failed");
	//else 
	massert(cudaSuccess == cudaMemcpy2D(dest, ldest*elemSize, src, ldsrc*elemSize, rows*elemSize, cols, kind),
	"CHLMemcpy2D: cudaMemcpy2D failed\n");


}

/*
void CoCMempy2DAsyncWrap3D(void* dest, long int ldest, void* src, long int ldsrc, long int rows, long int cols, short elemSize, int loc_dest, int loc_src, CQueue_p transfer_queue){
	// Convert 2d input (as CHLMemcpy2DAsync) to 3D for ...reasons.
	enum cudaMemcpyKind kind = cudaMemcpyDefault;
	cudaStream_t stream = *((cudaStream_t*) transfer_queue->backend_queue_ptr);
	cudaMemcpy3DParms* cudaMemcpy3DParms_p = (cudaMemcpy3DParms*) calloc(1, sizeof(cudaMemcpy3DParms));
	cudaMemcpy3DParms_p->extent = make_cudaExtent(rows*elemSize, cols, 1);
	cudaMemcpy3DParms_p->srcPtr = make_cudaPitchedPtr (src, ldsrc*elemSize, rows, cols );
	cudaMemcpy3DParms_p->dstPtr = make_cudaPitchedPtr (dest, ldest*elemSize, rows, cols );
	massert(cudaSuccess == cudaMemcpy3DAsync ( cudaMemcpy3DParms_p, stream) , "cudaMemcpy3DAsync failed\n");
}*/

template<typename VALUETYPE>
void CHLVecInit(VALUETYPE *vec, long long length, int seed, int loc)
{
  if (!vec) error("CHLVecInit: vec is not allocated (correctly)\n");
  if (loc < -1  || loc >= CHL_MEMLOCS) error("CHLVecInit: Invalid device id/location %d\n", loc);
  else if (loc >= CHL_WORKERS || loc == -1) CHLParallelVecInitHost(vec, length, seed);
  else {
	int prev_loc; cudaGetDevice(&prev_loc);

	//if (prev_loc != loc) warning("CHLVecInit: Initialized vector in other device (Previous device: %d, init in: %d)\n", prev_loc, loc);
    	cudaSetDevice(loc);
	curandGenerator_t gen;
	/* Create pseudo-random number generator */
	massert(curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT) == cudaSuccess,
          cudaGetErrorString(cudaGetLastError()));
	/* Set seed */
	massert(curandSetPseudoRandomGeneratorSeed(gen, seed) == cudaSuccess,
          cudaGetErrorString(cudaGetLastError()));
	if (typeid(VALUETYPE) == typeid(float))
	  massert(curandGenerateUniform(gen, (float*) vec, length) == cudaSuccess,
            cudaGetErrorString(cudaGetLastError()));
	else if (typeid(VALUETYPE) == typeid(double))
	  massert(curandGenerateUniformDouble(gen, (double*) vec, length) == cudaSuccess,
            cudaGetErrorString(cudaGetLastError()));
	CHLSyncCheckErr();
    	if (prev_loc != loc){
		//warning("CHLVecInit: Reseting device to previous: %d\n", prev_loc);
		cudaSetDevice(prev_loc);
	}
  }
  CHLSyncCheckErr();
}

template void CHLVecInit<double>(double *vec, long long length, int seed, int loc);
template void CHLVecInit<float>(float *vec, long long length, int seed, int loc);