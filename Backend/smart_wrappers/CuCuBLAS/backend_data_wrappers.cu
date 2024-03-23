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

int CHL_MEMLOCS = CHL_WORKERS + 2; 
//long PAGE_sz = sysconf (_SC_PAGESIZE);

char* mem_name(int idx){
	char* ans = (char*) malloc (10*sizeof(char));
	if(idx < 0) error("mem_name(%d) unsupported\n", idx);
	else if(idx < CHL_WORKERS) sprintf(ans, "Dev-%d", translate_mem_idx_to_hw(idx));
	else if (idx < CHL_MEMLOCS - 1) sprintf(ans, "Host");
	else if (idx == CHL_MEMLOCS - 1){
		sprintf(ans, "Numa-inter");
	}
	else error("mem_name(%d) unsupported\n", idx);
	return ans;
}

int get_hostmem_idx(void* addr){
	int is_interleaved = -1;
	get_mempolicy(&is_interleaved, NULL, 0, addr, MPOL_F_ADDR);
	if(is_interleaved == MPOL_INTERLEAVE) return CHL_MEMLOCS -1;
	//int numa_node = -1;
	//get_mempolicy(&numa_node, NULL, 0, addr, MPOL_F_NODE | MPOL_F_ADDR);
	//return numa_node; 
	return CHL_WORKERS;
}

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

void *gpu_malloc(long long count) {
	void *ret;
	massert(cudaMalloc(&ret, count) == cudaSuccess,
		cudaGetErrorString(cudaGetLastError()));
	return ret;
}

void pin_mem_wrap(void** ptr, long long bytes){
	cudaHostRegister(*ptr,bytes,cudaHostRegisterPortable);
}

void *pin_malloc(long long count, short W_flag) {
	void *ret;
	ret = malloc(count);
	cudaHostRegister(ret,count,cudaHostRegisterPortable);
	//massert(cudaHostAlloc ( &ret, count, cudaHostAllocDefault)  == cudaSuccess,
	//cudaGetErrorString(cudaGetLastError()));
	return ret;
}

void *numa_inter_pin_malloc(long long count, short W_flag){
	void *ret;
	ret = numa_alloc_interleaved(count);
	cudaHostRegister(ret,count,cudaHostRegisterPortable);
	return ret;
}

void *numa_bind_pin_malloc(long long count, int node_num, short W_flag){
	void *ret;
	ret = numa_alloc_onnode(count, translate_mem_idx_to_hw(node_num));
	cudaHostRegister(ret,count,cudaHostRegisterPortable);
	return ret;
}

void* CHLMalloc(long long bytes, int loc, short W_flag){
	void *ptr = NULL;
	if (loc == CHL_MEMLOCS - 1) {
#ifdef CLDEBUG
		fprintf(stderr, "Allocating %lld bytes with interleaved NUMA alloc and pinning...\n", bytes);
#endif
		ptr = numa_inter_pin_malloc(bytes, W_flag);
	}
	else if (loc == CHL_WORKERS) {
#ifdef CLDEBUG
		fprintf(stderr, "Allocating %lld bytes with simple malloc and pinning...\n", bytes);
#endif
		ptr = numa_inter_pin_malloc(bytes, W_flag);		
	}
	else if (loc > CHL_WORKERS) {
//#ifdef CLDEBUG
		warning("Allocating %lld bytes to NUM MEM %d [%s]...\n", bytes, loc, mem_name(loc));
//#endif
		ptr = numa_bind_pin_malloc(bytes, loc, W_flag);
	}
	else if (loc == -1){
#ifdef CLDEBUG
    	fprintf(stderr, "Allocating %lld bytes with interleaved NUMA alloc without pinning...\n", bytes);
#endif
		ptr = numa_alloc_interleaved(bytes);
	}
	else if (loc == -2){
#ifdef CLDEBUG
    	fprintf(stderr, "Allocating %lld bytes with simple malloc without pinning...\n", bytes);
#endif
		ptr = malloc(bytes);
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
		if (prev_loc != loc) cudaSetDevice(prev_loc);
	}
	else error("CHLMalloc: Invalid device id/location %d\n", loc);
	return ptr;
}

void gpu_free(void *gpuptr) {
	massert(cudaFree(gpuptr) == cudaSuccess,
			cudaGetErrorString(cudaGetLastError()));
}

void pin_free(void *ptr) {
	//massert(cudaFreeHost(gpuptr) == cudaSuccess,
	//        cudaGetErrorString(cudaGetLastError()));
	cudaHostUnregister(ptr);
	free(ptr);
}

void numa_pin_free(void *ptr, long long bytes) {
	cudaHostUnregister(ptr);
	numa_free(ptr,bytes);
}

void CHLFree(void * ptr, long long bytes, int loc){
	//if (??? == loc) free(ptr);
	if (loc > CHL_WORKERS && loc < CHL_MEMLOCS) numa_pin_free(ptr, bytes);
	else if (loc == CHL_WORKERS) pin_free(ptr);
	else if (loc == -1) numa_free(ptr,bytes);
	else if (loc == -2 ) free(ptr);
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
}

short CHLGetPtrLoc(void * in_ptr)
{
// This is legacy code for CUDA 9.2 <<. It should not be used due to CUDA ptr_att back-end struct changes in latest versions
#ifdef CUDA_9_WRAPPER_MESS
  error("CHLGetPtrLoc(9.2 version, ptr =%p): Not supported anymore\n", in_ptr);	
#else
	int loc = -42;
	cudaPointerAttributes ptr_att;
	if (cudaSuccess != cudaPointerGetAttributes(&ptr_att, in_ptr)) error("CHLGetPtrLoc(cuda 10+ version, ptr =%p): "
	"Pointer not visible to CUDA, host alloc or error\n", in_ptr);
	else if (ptr_att.type == cudaMemoryTypeDevice) loc = ptr_att.device;
	// TODO: Unified memory is considered available in the GPU as cuBLASXt ( not bad, not great)
	else if (ptr_att.type == cudaMemoryTypeManaged) loc = ptr_att.device;
	else{ // if (ptr_att.type == cudaMemoryTypeHost){
		// Note: This does not take into account if memory is pinned or not because it is not needed!
		loc = get_hostmem_idx();
	}
	return loc;
#endif
}

void CHLMemcpy(void* dest, void* src, long long bytes, int loc_dest, int loc_src)
{
	massert(loc_dest >= -2 && loc_dest < CHL_MEMLOCS, "CHLMemcpy: Invalid destination device: %d\n", loc_dest);
	massert(loc_src >= -2 && loc_src < CHL_MEMLOCS, "CHLMemcpy: Invalid source device: %d\n", loc_src);

#ifdef DEBUG
	if (loc_src == loc_dest) warning("CHLMemcpy(dest=%p, src=%p, bytes=%lld, loc_dest=%d, loc_src=%d): Source location matches destination\n",
	dest, src, bytes, loc_dest, loc_src);
#endif
	massert(CUBLAS_STATUS_SUCCESS == cudaMemcpy(dest, src, bytes, cudaMemcpyDefault), "CHLMemcpy: cudaMemcpy from device src=%d to dest=%d failed\n", loc_src, loc_dest);
	CHLSyncCheckErr();
}

/*void CHLMemcpyAsync(void* dest, void* src, long long bytes, int loc_dest, int loc_src, CQueue_p transfer_queue)
{
	massert(loc_dest >= -2 && loc_dest < CHL_MEMLOCS, "CHLMemcpy: Invalid destination device: %d\n", loc_dest);
	massert(loc_src >= -2 && loc_src < CHL_MEMLOCS, "CHLMemcpy: Invalid source device: %d\n", loc_src);

	cudaStream_t stream = *((cudaStream_t*) transfer_queue->backend_queue_ptr);

	massert(loc_dest >= -1 && loc_dest < CHL_MEMLOCS, "CHLMemcpyAsync: Invalid destination device: %d\n", loc_dest);
	massert(loc_src >= -1 && loc_src < CHL_MEMLOCS, "CHLMemcpyAsync: Invalid source device: %d\n", loc_src);


	if (loc_src == loc_dest) warning("CHLMemcpyAsync(dest=%p, src=%p, bytes=%lld, loc_dest=%d, loc_src=%d): Source location matches destination\n",
	dest, src, bytes, loc_dest, loc_src);
	massert(cudaSuccess == cudaMemcpyAsync(dest, src, bytes, cudaMemcpyDefault, stream),
	"CHLMemcpy2D: cudaMemcpyAsync failed\n");
	//CHLSyncCheckErr();
}*/

void CHLMemcpy2D(void* dest, long int ldest, void* src, long int ldsrc, long int rows, long int cols, short elemSize, int loc_dest, int loc_src){
	short lvl = 6;
#ifdef DDEBUG
	lprintf(lvl, "CHLMemcpy2D(dest=%p, ldest =%zu, src=%p, ldsrc = %zu, rows = %zu, cols = %zu, elemsize = %d, loc_dest = %d, loc_src = %d)\n",
		dest, ldest, src, ldsrc, rows, cols, elemSize, loc_dest, loc_src);
#endif
	massert(loc_dest >= -2 && loc_dest < CHL_MEMLOCS, "CHLMemcpy2D: Invalid destination device: %d\n", loc_dest);
	massert(loc_src >= -2 && loc_src < CHL_MEMLOCS, "CHLMemcpy2D: Invalid source device: %d\n", loc_src);

	if (loc_src == loc_dest) warning("CHLMemcpy2D(dest=%p, ldest =%zu, src=%p, ldsrc = %zu, rows=%zu, cols=%zu, elemSize =%d, loc_dest=%d, loc_src=%d): Source location matches destination\n",
	dest, ldest, src, ldsrc, rows, cols, elemSize, loc_dest, loc_src);
	//if (loc_src == -1 && loc_dest >=0) massert(CUBLAS_STATUS_SUCCESS == cublasSetMatrix(rows, cols, elemSize, src, ldsrc, dest, ldest), "CHLMemcpy2DAsync: cublasSetMatrix failed\n");
	//else if (loc_src >=0 && loc_dest == -1) massert(CUBLAS_STATUS_SUCCESS == cublasGetMatrix(rows, cols, elemSize, src, ldsrc, dest, ldest),  "CHLMemcpy2DAsync: cublasGetMatrix failed");
	//else 
	massert(cudaSuccess == cudaMemcpy2D(dest, ldest*elemSize, src, ldsrc*elemSize, rows*elemSize, cols, cudaMemcpyDefault),
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
  if (loc < -2  || loc >= CHL_MEMLOCS) error("CHLVecInit: Invalid device id/location %d\n", loc);
  else if (loc >= CHL_WORKERS || loc < 0 ) CHLParallelVecInitHost(vec, length, seed);
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