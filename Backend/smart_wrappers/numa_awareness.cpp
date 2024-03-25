///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief Some convenient C/C++ utilities for CHL.
///

#include "chl_smart_wrappers.hpp"

#include <cstdlib>
#include <cmath>

#include <unistd.h>
#include <numa.h>
#include <numaif.h>
#include <omp.h>


int DEV_NUM = -1, NUMA_HW_NUM = -1, NIC_NUM = -1, HW_THREADS = -1, NIC_AT_DEV[32] = {0}, NUMA_AT_DEV[32] = {0};
int CHL_WORKERS = -1, CHL_MEMLOCS = -1, CHL_HWNUMA_AT_MEMLOC[32] = {0}, 
	CHL_WORKER_CLOSE_TO_MEMLOC[32] = {0}, CHL_MEMLOC_CLOSE_TO_WORKER[64] = {0}, 
	CHL_INPUT_QUEUES_CASE_IDS[32][MAX_WORKER_CONFIG] = {0}, CHL_OUTPUT_QUEUES_CASE_IDS[32][MAX_WORKER_CONFIG] = {0};

int sconfig_load(){
    char *filename = (char *) malloc(1024 * sizeof(char));
    sprintf(filename, "%s", SCONFIG_PATH);
	FILE* fp = fopen(filename,"r");
	if(!fp) error("sconfig_load(): Failed to open file %s\n", filename); 
	fscanf(fp, "DEV_NUM=%d\n", &DEV_NUM);
	fscanf(fp, "NUMA_HW_NUM=%d\n", &NUMA_HW_NUM);
	fscanf(fp, "NIC_NUM=%d\n", &NIC_NUM);
	fscanf(fp, "nic_at_dev=[");
	if(DEV_NUM > 32) error("sconfig_load(): CHL does not support systems with > 32 devices (%d here), needs metadata storage extension\n", DEV_NUM);
	for(int idx = 0; idx < DEV_NUM; idx ++)
		fscanf(fp, " %d", &(NIC_AT_DEV[idx]));
	fscanf(fp, " ]\nnuma_at_dev=[");
	for(int idx = 0; idx < DEV_NUM; idx ++)
		fscanf(fp, " %d", &(NUMA_AT_DEV[idx]));
	fscanf(fp, " ]\nHW_THREADS=%d\n", &HW_THREADS);
	fclose(fp);
	free(filename);
	fprintf(stderr, "----------------------------------------------\n");
	fprintf(stderr, "sconfig_load() -> Loading system config file:\nDEV_NUM = %d\nNUMA_HW_NUM = %d\nNIC_NUM = %d\nHW_THREADS = %d\n"
		"NIC_AT_DEV[%d] = %s\nNUMA_AT_DEV[%d] = %s\n", DEV_NUM, NUMA_HW_NUM, NIC_NUM, HW_THREADS, DEV_NUM, 
		printlist(NIC_AT_DEV, DEV_NUM), DEV_NUM, printlist(NUMA_AT_DEV, DEV_NUM));
	
	int bind_list[NUMA_HW_NUM] = {0};
	for(int idx = 0; idx < NUMA_HW_NUM; idx ++) bind_list[idx] = -1;
	for(int idx = 0; idx < DEV_NUM; idx ++){
		if(bind_list[NUMA_AT_DEV[idx]]!=-1 && bind_list[NUMA_AT_DEV[idx]]!=NIC_AT_DEV[idx]) error("sconfig_load() bind_list[%d] = %d already defined (new val = %d)\n", 
			NUMA_AT_DEV[idx], bind_list[NUMA_AT_DEV[idx]], NIC_AT_DEV[idx]);
		else if(bind_list[NUMA_AT_DEV[idx]] == -1) bind_list[NUMA_AT_DEV[idx]] = NIC_AT_DEV[idx];
	}

	CHL_WORKERS = DEV_NUM;
	CHL_MEMLOCS = CHL_WORKERS + std::min(NUMA_HW_NUM, NIC_NUM) + 1;//NUMA_HW_NUM + 1;//
	int bind_ctr = 0;
	for(int idx = 0; idx < CHL_MEMLOCS; idx ++) 
		if(idx < CHL_WORKERS) CHL_HWNUMA_AT_MEMLOC[idx] = -1; // MEMLOC is a device
		else if(idx == CHL_MEMLOCS - 1) CHL_HWNUMA_AT_MEMLOC[idx] = -42; // MEMLOC is interleaved memory
		else{
			while(bind_list[bind_ctr] == -1) bind_ctr++;
			CHL_HWNUMA_AT_MEMLOC[idx] = bind_ctr++;
		}
	for(int idx = 0; idx < CHL_WORKERS; idx++) 
	CHL_WORKER_CLOSE_TO_MEMLOC[idx] =  translate_hw_numa_to_mem_idx(NUMA_AT_DEV[idx]);
	for(int idx = 0; idx < CHL_MEMLOCS; idx++){
		for(int idx1 = 0; idx1 < CHL_WORKERS; idx1++) if(CHL_WORKER_CLOSE_TO_MEMLOC[idx1] == idx)
			CHL_MEMLOC_CLOSE_TO_WORKER[idx] = idx1;
	}
	fprintf(stderr, "\nInitializing CHL grid amalgamation dimensions:\nCHL_WORKERS = %d\nCHL_MEMLOCS = %d\nCHL_WORKER_CLOSE_TO_MEMLOC[%d] = %s\nCHL_MEMLOC_CLOSE_TO_WORKER[%d] = %s\nCHL_HWNUMA_TO_MEMLOC[%d] = %s\n| ", 
		CHL_WORKERS, CHL_MEMLOCS, CHL_WORKERS, printlist(CHL_WORKER_CLOSE_TO_MEMLOC, CHL_WORKERS), 
			CHL_MEMLOCS, printlist(CHL_MEMLOC_CLOSE_TO_WORKER, CHL_MEMLOCS), CHL_MEMLOCS, printlist(CHL_HWNUMA_AT_MEMLOC, CHL_MEMLOCS));
	for(int idx = 0; idx < CHL_MEMLOCS; idx ++) fprintf(stderr, "%s | ", mem_name(idx));

	int active_unit_num = CHL_WORKERS, active_unit_id_list[CHL_WORKERS], rev_active_unit_num = CHL_WORKERS, rev_active_unit_id_list[CHL_WORKERS];
	for(int idx = 0; idx < std::min(NIC_NUM, CHL_WORKERS); idx++){
		int worker_num = idx + 1;
		for(int idx1 = 0; idx1 < MAX_WORKER_CONFIG; idx1++) CHL_INPUT_QUEUES_CASE_IDS[idx][idx1] = CHL_OUTPUT_QUEUES_CASE_IDS[idx][idx1] = -1;
			int ctr = 0, ctr_glob = 0, used_nic_list[NIC_NUM], prev_used_workers[CHL_WORKERS] = {0}, used_workers[CHL_WORKERS] = {0};
			for(int idx1 = 0; idx1 < CHL_WORKERS; idx1++) prev_used_workers[idx1] = -1;
			for(int idx1 = 0; idx1 < MAX_WORKER_CONFIG; idx1++){
				ctr = 0;
				for(int idx2 = 0; idx2 < CHL_WORKERS; idx2++) used_workers[idx2] = -1;
				for(int idx2 = 0; idx2 < NIC_NUM; idx2++) used_nic_list[idx2] = -1;
				for(int idx2 = 0; idx2 < CHL_WORKERS; idx2++){
					if(ctr == worker_num) break;
					if(!is_in_list(NIC_AT_DEV[idx2], used_nic_list, ctr) && !is_in_list(idx2, prev_used_workers, ctr_glob)){
						used_workers[ctr] = idx2;
						prev_used_workers[ctr_glob++] = idx2;
						used_nic_list[ctr++] = NIC_AT_DEV[idx2];
					}
				}
				if(ctr == worker_num) CHL_INPUT_QUEUES_CASE_IDS[idx][idx1] = CHL_OUTPUT_QUEUES_CASE_IDS[idx][idx1] = translate_unit_list_to_binary(used_workers, ctr);
			}
	}
	translate_binary_to_unit_list(CHL_INPUT_QUEUES_CASE_IDS[NIC_NUM-1][0], &active_unit_num, active_unit_id_list);
	for(int idx = NIC_NUM; idx < CHL_WORKERS; idx++){
		int worker_num = idx + 1, different_sets = MAX_WORKER_CONFIG/2;
		int ctr, ctr_glob = 0, prev_used_workers[CHL_WORKERS] = {0}, used_workers[CHL_WORKERS] = {0};
		for(int idx1 = 0; idx1 < CHL_WORKERS; idx1++) prev_used_workers[idx1] = -1;
		for(int idx1 = 0; idx1 < different_sets; idx1++){
			ctr = 0;
			for(int idx1 = 0; idx1 < CHL_WORKERS; idx1++) used_workers[idx1] = -1;
			for(int idx1 = 0; idx1 < CHL_WORKERS; idx1++) 
				if(is_in_list(idx1, active_unit_id_list, active_unit_num)) used_workers[ctr++] = idx1;
			int cand = 0; 
			while(ctr < worker_num && cand < CHL_WORKERS){
				if(!is_in_list(cand, used_workers, ctr) && !is_in_list(cand, prev_used_workers, ctr_glob)){
						used_workers[ctr++] = cand;
						prev_used_workers[ctr_glob++] = cand;
				}
				cand++;
			}
			if(ctr == worker_num){
				CHL_INPUT_QUEUES_CASE_IDS[idx][idx1*2] = CHL_OUTPUT_QUEUES_CASE_IDS[idx][idx1*2] = translate_unit_list_to_binary(used_workers, ctr);
				CHL_INPUT_QUEUES_CASE_IDS[idx][idx1*2 + 1] =  translate_unit_list_to_binary(used_workers, ctr);
				CHL_OUTPUT_QUEUES_CASE_IDS[idx][idx1*2 + 1] = CHL_OUTPUT_QUEUES_CASE_IDS[NIC_NUM-1][0];
			}
		}
	}
	for(int idx = 0; idx < CHL_WORKERS; idx++){
		fprintf(stderr, "\nGenerated configurations for distribution to %d workers:\n", idx + 1);
		for(int idx2 = 0; idx2 < MAX_WORKER_CONFIG; idx2++)
		if(CHL_OUTPUT_QUEUES_CASE_IDS[idx][idx2] == -1) break;
		else{
			translate_binary_to_unit_list(CHL_INPUT_QUEUES_CASE_IDS[idx][idx2], &active_unit_num, active_unit_id_list);
			translate_binary_to_unit_list(CHL_OUTPUT_QUEUES_CASE_IDS[idx][idx2], &rev_active_unit_num, rev_active_unit_id_list);
			fprintf(stderr, "Input queues = %s, Output queues = %s\n",
			printlist(active_unit_id_list, active_unit_num), 
			printlist(rev_active_unit_id_list, rev_active_unit_num));
		}
	}

	fprintf(stderr, "\n----------------------------------------------\n");

	return 1;
}

int load_flag = sconfig_load();
/*
MEMMetadata::MEMMetadata(void* mem_ptr_in, int* slice_memloc_in, 
	long long slice_size_in, long long slice_num_in){
#ifdef DCLDEBUG
	fprintf(stderr, "MEMMetadata::MEMMetadata(%p, %s, %lld, %lld)\n", mem_ptr_in, 
		printlist<int>(slice_memloc_in, (int)slice_num_in), slice_size_in, slice_num_in);
#endif
	mem_ptr = mem_ptr_in;
	slice_size = slice_size_in;
	slice_num = slice_num_in; 
	slices_memloc = (int*) malloc (slice_num*sizeof(int));
	for (long long itter = 0; itter < slice_num; itter++) 
		slices_memloc[itter] = slice_memloc_in[itter];
}
	
MEMMetadata::~MEMMetadata(){
	free(slices_memloc);
}

void MEMMetadata::get_slices(int** slice_copy_p, int* slice_num_p){
	*slice_num_p = slice_num;
	*slice_copy_p = (int*) malloc (slice_num*sizeof(int));
	for (long long itter = 0; itter < slice_num; itter++) 
		(*slice_copy_p)[itter] = slices_memloc[itter];
}

int MEMMetadata::get_relative_memloc(void* sub_ptr){
	long long offset = ((char*)sub_ptr) - ((char*)mem_ptr);
	if(offset < 0 || offset >= slice_size*slice_num) return -1;
	else return slices_memloc[offset/slice_size];

}

void MEMMetadata::print_slices(){
	fprintf(stderr, "MEMMetadata::print_slices(): mem_ptr = %p allocated in %lld slices of size %lld KB\n-> [ ", mem_ptr, slice_num, slice_size/1024);
	for (long long itter = 0; itter < slice_num; itter++) 
		fprintf(stderr, "%d ", slices_memloc[itter]);
	fprintf(stderr, "]\n");
}
*/

char* mem_name(int idx){
	char* ans = (char*) malloc (10*sizeof(char));
	if(idx < 0) error("mem_name(%d) unsupported\n", idx);
	else if(idx < CHL_WORKERS) sprintf(ans, "Dev-%d", translate_mem_idx_to_hw(idx));
	else if (idx < CHL_MEMLOCS - 1) sprintf(ans, "NuN-%d", translate_mem_idx_to_hw(idx));
	else if (idx == CHL_MEMLOCS - 1) sprintf(ans, "Inter");
	else error("mem_name(%d) unsupported\n", idx);
	return ans;
}

int get_hw_numa_num(){
	return NUMA_HW_NUM; 
}

int get_hw_numa_idx(void* addr){
	int is_interleaved = -1;
	get_mempolicy(&is_interleaved, NULL, 0, addr, MPOL_F_ADDR);
	if(is_interleaved == MPOL_INTERLEAVE) return -42;
	int numa_node = -1;
	get_mempolicy(&numa_node, NULL, 0, addr, MPOL_F_NODE | MPOL_F_ADDR);
	return numa_node; 
}

void* search_sub_addrs_at_memloc(void* addr, int memloc, long long chunk_offset, long long size, long long* found_chunk_idx){
	char* chunk_addrs = NULL;
	long long chunk_num = size/chunk_offset + ((size%chunk_offset)? 1: 0); 
	for(long long idx = 0; idx < chunk_num; idx++){
		chunk_addrs = ((char*) addr) + idx*chunk_offset;
		if (translate_hw_numa_to_mem_idx(get_hw_numa_idx(chunk_addrs)) == memloc){
			fprintf(stderr, "search_sub_addrs_at_memloc(addr=%p,loc=%d,chunk=%lld): Found chunk_addrs = %p at chunk %lld\n", 
				addr, memloc, chunk_offset, chunk_addrs, idx);
			*found_chunk_idx = idx;
			return chunk_addrs;
		}
	}
	return NULL;
}

void get_hw_numa_list(void* addr, long long chunk_offset, long long chunk_num, long long* numa_list_num, int* numa_list){
	char* chunk_addrs = NULL;
	for(int idx = 0; idx < NUMA_HW_NUM; idx++) numa_list_num[idx] = 0;
	for(long long idx = 0; idx < chunk_num; idx++){
		chunk_addrs = ((char*) addr) + idx*chunk_offset;
		numa_list[idx] = //translate_hw_numa_to_mem_idx
			(get_hw_numa_idx(chunk_addrs));
		numa_list_num[get_hw_numa_idx(chunk_addrs)]++;
	}
}

int translate_mem_idx_to_hw(int idx){
	if(idx < CHL_WORKERS) return idx;
	else if(idx == CHL_MEMLOCS - 1) return idx;
	else if(idx < 0 || idx >= CHL_MEMLOCS) error("translate_mem_idx_to_hw(%d): Invalid idx\n", idx);
	return CHL_HWNUMA_AT_MEMLOC[idx]; 
}

int translate_hw_numa_to_mem_idx(int hw_numa_idx){
	if(-42 == hw_numa_idx) return CHL_MEMLOCS-1; 
	for (int idx = CHL_WORKERS; idx < CHL_MEMLOCS-1; idx++)
	 	if(CHL_HWNUMA_AT_MEMLOC[idx] == hw_numa_idx) return idx;
	return -1; 
}

/*int scan_allocated_data_for_ptr_loc(void* ptr){
#ifdef DCLDEBUG
	fprintf(stderr, "scan_allocated_data_for_ptr_loc(%p)\n", ptr);
#endif
	for (long long itter = 0; itter < allocated_data_num; itter++){
		int loc_cand = allocated_data[itter]->get_relative_memloc(ptr); 
		if (loc_cand != -1){
#ifdef DCLDEBUG
	fprintf(stderr, "ptr = %p found in allocated_data[%lld] (mem_ptr = %p): memloc = %d (name: %s)\n", 
		ptr, itter, allocated_data[itter]->mem_ptr, loc_cand, mem_name(loc_cand));
#endif
			return loc_cand;
		}
	}
	return -1; 
}*/

void* CHLMallocHostTouchDefault(long long bytes){
	void *ret;
	posix_memalign(&ret, PAGE_sz*PAGE_sz, bytes);
	long long chunk_num = bytes/PAGE_sz;
	long long chunk_team_num = chunk_num/HW_THREADS; 
	#pragma omp parallel for schedule(static,chunk_team_num)
	for(long chunk_idx = 0; chunk_idx < chunk_num; chunk_idx++)
		((char*)ret)[chunk_idx*PAGE_sz] = 0;
	pin_mem_wrap(&ret, bytes);
	return ret;
}

void* CHLMallocHostTouchSerial(long long bytes){
	void *ret;
	if(bytes%PAGE_sz) error("CHLMallocHostTouchSmart was not designed to allocate non-multiples of PAGE size (%ld)\n", PAGE_sz);
	const char* env_bind = std::getenv("OMP_PROC_BIND");
	if (strcmp(env_bind,"TRUE")) warning("CHLMallocHostTouchSmart: OMP threads are not pinned to HW, run [export OMP_PROC_BIND=TRUE] for Smart touch.\n"
	"Defaulting to interleaved allocation, performance might degrade\n");
	int env_threads = atoi(std::getenv("OMP_NUM_THREADS"));
    if (env_threads!= HW_THREADS) warning("CHLMallocHostTouchSmart: OMP_NUM_THREADS = %d instead of %d that the framework was configured for, run [export OMP_NUM_THREADS=%d] for Smart touch.\n"
	"Defaulting to interleaved allocation, performance might degrade\n", env_threads, HW_THREADS, HW_THREADS);
	//ret = malloc(bytes);
	// Do not ask questions.... using less addr alingment (actually anything >=512 works in Karolina) results in a few runaway pages to other numa nodes. 
	posix_memalign(&ret, PAGE_sz*PAGE_sz, bytes);
	long long chunk_num = bytes/PAGE_sz;
	int numa_memlocs = CHL_MEMLOCS-1-CHL_WORKERS;
	int alloc_per_thread =  NUMA_HW_NUM/(numa_memlocs);
	long long chunk_alloc_per_thread = chunk_num/numa_memlocs;
	#pragma omp parallel for schedule(static,1)
	for(int thread_id = 0; thread_id < HW_THREADS; thread_id++){
		int thread_idx = omp_get_thread_num(); // This *should* be the same with thread_id for OMP_PROC_BIND=TRUE
		if(thread_idx!=thread_id) warning("CHLMallocHostTouchSmart: Will not work as intended since OMP threads are not pinned to HW correctly (?)\n");
		int thread_core = thread_idx/(HW_THREADS/NUMA_HW_NUM);
		int thread_ctr = thread_idx%(HW_THREADS/NUMA_HW_NUM);
		int thread_mem_idx = thread_core/alloc_per_thread; 
		if (is_in_list(thread_core, CHL_HWNUMA_AT_MEMLOC, CHL_MEMLOCS) && thread_ctr == 0){
			fprintf(stderr, "Thread %d (core = %d, ctr = %d, mem_idx = %d) selected as allocator\n", thread_idx, thread_core, thread_ctr, thread_mem_idx);
			for(long long chunk_idx = 0; chunk_idx < chunk_alloc_per_thread; chunk_idx++){
				((unsigned char*)ret)[chunk_idx*PAGE_sz + chunk_alloc_per_thread*thread_mem_idx*PAGE_sz] = 0;
			}
		}
	}
/*
	for(long chunk_idx = 0; chunk_idx < chunk_num; chunk_idx++){
		int thread_idx = omp_get_thread_num();
		int thread_team_numa_node = thread_idx/(HW_THREADS/NUMA_HW_NUM);
		int chunk_idx_area = 
		if (is_in_list(thread_team_numa_node, CHL_HWNUMA_AT_MEMLOC, CHL_MEMLOCS)){
			((char*)ret)[chunk_idx*PAGE_sz] = 0;
			((char*)ret)[chunk_idx*PAGE_sz + chunk_team_num*16*PAGE_sz] = 0;

		}
	}
*/
	pin_mem_wrap(&ret, bytes);
	return ret;
}

void* CHLMallocHostTouchSmart(long long dim1, long long dim2, int elemSize, char transpose){
	long long GridSz1 = dim1/PAGE_sz;
	long long GridSz2 = dim2/PAGE_sz;
	if (dim1*elemSize%PAGE_sz > 0) GridSz1++;
	if (dim2%PAGE_sz > 0) GridSz2++;
	int numa_memlocs = CHL_MEMLOCS-1-CHL_WORKERS;
	int D1_parts = std::sqrt(numa_memlocs);
  	int D2_parts = D1_parts;
  	if (D1_parts ==0) { D2_parts = numa_memlocs; D1_parts = 1; }
  	else {
    	/* find the most square decomposition of numa_memlocs in D1_parts x D2_parts */
		int g;
		for (g = D1_parts+1; g>0; --g)
		if (numa_memlocs % g == 0) break;
		if (g==0) { D1_parts = numa_memlocs; D2_parts = 1; }
		//if (g==0) { D1_parts = 1; D2_parts = numa_memlocs; }
		else { D1_parts = g; D2_parts = numa_memlocs/g; }
	}
	//int temp = D2_parts;
    //D2_parts = D1_parts;
    //D1_parts = temp;
	void *ret;
	const char* env_bind = std::getenv("OMP_PROC_BIND");
	if (env_bind && strcmp(env_bind,"TRUE")) warning("CHLMallocHostTouchSmart: OMP threads are not pinned to HW, run [export OMP_PROC_BIND=TRUE] for Smart touch.\n"
	"Defaulting to interleaved allocation, performance might degrade\n");
	const char* env_threads_c = std::getenv("OMP_NUM_THREADS");
	int env_threads = HW_THREADS;
	if (!env_threads_c) warning("CHLMallocHostTouchSmart: OMP_NUM_THREADS not defined, run [export OMP_NUM_THREADS=%d] for Smart touch.\n", HW_THREADS);
	else env_threads = atoi(env_threads_c);
	if (env_threads!= HW_THREADS) warning("CHLMallocHostTouchSmart: OMP_NUM_THREADS = %d instead of %d that the framework was configured for, run [export OMP_NUM_THREADS=%d] for Smart touch.\n"
	"Defaulting to interleaved allocation, performance might degrade\n", env_threads, HW_THREADS, HW_THREADS);
	//ret = malloc(bytes);
	// Do not ask questions.... using less addr alingment (actually anything >=512 works in Karolina) results in a few runaway pages to other numa nodes. 
	posix_memalign(&ret, PAGE_sz*PAGE_sz, GridSz1*GridSz2*PAGE_sz*PAGE_sz*elemSize);
	#pragma omp parallel for schedule(static,1)
	for(int thread_id = 0; thread_id < HW_THREADS; thread_id++){
		int thread_idx = omp_get_thread_num(); // This *should* be the same with thread_id for OMP_PROC_BIND=TRUE
		if(thread_idx!=thread_id) warning("CHLMallocHostTouchSmart: Will not work as intended since OMP threads are not pinned to HW correctly (?)\n");
		int thread_core = thread_idx/(HW_THREADS/NUMA_HW_NUM);
		int thread_ctr = thread_idx%(HW_THREADS/NUMA_HW_NUM);
		int thread_mem_idx = thread_core/(NUMA_HW_NUM/numa_memlocs);
		int thread_mem_idx_row = thread_mem_idx/D2_parts;
		int thread_mem_idx_col = thread_mem_idx%D2_parts;
		long long tile_offset = 0;
		if (is_in_list(thread_core, CHL_HWNUMA_AT_MEMLOC, CHL_MEMLOCS) && thread_ctr == 0){
			fprintf(stderr, "Thread %d (core = %d, ctr = %d, mem_idx = %d) selected as allocator if area (Ser = %d, Grid[row,col] = [%d,%d]) \n", 
				thread_idx, thread_core, thread_ctr, thread_mem_idx, thread_mem_idx, thread_mem_idx_row, thread_mem_idx_col);
			for (int itt1 = 0; itt1 < GridSz1; itt1++){
				for (int itt2 = 0 ; itt2 < GridSz2; itt2++){
					/// For column major format assumed
					if (transpose == 'N'){
						tile_offset = (itt1*PAGE_sz*elemSize + itt2*dim1*elemSize*PAGE_sz);
					}
					else if (transpose == 'T'){
						tile_offset = (itt1*dim2*elemSize*PAGE_sz + itt2*PAGE_sz);
					}
					else error("CHLMallocHostTouchSmart: Unknown transpose type\n");
					if (thread_mem_idx_row == itt1/(GridSz1/D1_parts)
					&& thread_mem_idx_col == itt2/(GridSz2/D2_parts))
						for(int itter_idx = 0; itter_idx < PAGE_sz; itter_idx++)
							for(int elem_idx = 0; elem_idx < PAGE_sz*elemSize; elem_idx++)
								((unsigned char*)ret)[tile_offset + itter_idx*dim1*elemSize + elem_idx] = 0;
				}
			}
		}
	}
	pin_mem_wrap(&ret, GridSz1*GridSz2*PAGE_sz*PAGE_sz*elemSize);
	return ret;
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