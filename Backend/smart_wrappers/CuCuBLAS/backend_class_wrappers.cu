///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief Some CUDA function calls with added error-checking
///

#include <cuda.h>
#include <cublas_v2.h>

#include "smart_wrappers.hpp"
#include "backend_wrappers.hpp"

int Event_num_loc[64] = {0};
int Queue_num_loc[64] = {0};

/*****************************************************/
/// Event Status-related functions

const char* print_event_status(event_status in_status){
	switch(in_status){
		case(UNRECORDED):
			return "UNRECORDED";
		case(RECORDED):
			return "RECORDED";
		case(COMPLETE):
			return "COMPLETE";
		case(CHECKED):
			return "CHECKED";
		case(GHOST):
			return "GHOST";
		default:
			error("print_event_status: Unknown state\n");
	}
}

/*****************************************************/
/// Command queue class functions

void CommandQueue_comm_init(CQueue_p myself){
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] |-----> CommandQueue_comm_init(Queue(%d))\n", myself->dev_id, myself->id);
#endif
	;
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] <-----| CommandQueue_comm_init(Queue(%d))\n", myself->dev_id, myself->id);
#endif
}

void CommandQueue_comm_delete(CQueue_p myself){
	free(myself->backend_comm_md);
}

void CommandQueue_comp_init(CQueue_p myself){
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] |-----> CommandQueue_comp_init(Queue(%d))\n", myself->dev_id, myself->id);
#endif
	cudaStream_t stream = *((cudaStream_t*) myself->backend_queue_ptr);
	myself->backend_comp_md = malloc(sizeof(cublasHandle_t));
	massert(CUBLAS_STATUS_SUCCESS == cublasCreate((cublasHandle_t*) myself->backend_comp_md),
		"CommandQueue_comp_init(%d): cublasCreate failed\n", myself->dev_id);
	massert(CUBLAS_STATUS_SUCCESS == cublasSetStream(*((cublasHandle_t*) myself->backend_comp_md), stream),
		"CommandQueue_comp_init(%d): cublasSetStream failed\n", myself->dev_id);
	if(WS_SZ >=0){
	void* local_ws = NULL; if(WS_SZ) cudaMalloc(&local_ws, WS_SZ); 
	massert(CUBLAS_STATUS_SUCCESS == cublasSetWorkspace(*((cublasHandle_t*) 
		myself->backend_comp_md), local_ws, WS_SZ), 
		"CommandQueue_comp_init(%d): cublasSetWorkspace failed\n", myself->dev_id);
	}
	//warning("FIXME: Running on limited SMs, custom stuff, beware this in not a drill\n");
	//massert(CUBLAS_STATUS_SUCCESS == cublasSetSmCountTarget(*((cublasHandle_t*) myself->backend_comp_md), 1),
	//	"CommandQueue::CommandQueue(%d): cublasSetSmCountTarget failed\n", myself->dev_id);
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] <-----| CommandQueue_comp_init(Queue(%d))\n", myself->dev_id, myself->id);
#endif
}

void CommandQueue_comp_delete(CQueue_p myself){
	cublasHandle_t handle = *((cublasHandle_t*) myself->backend_comp_md);
	massert(CUBLAS_STATUS_SUCCESS == cublasDestroy(handle),
		"CommandQueue_comp_delete - cublasDestroy(handle) failed\n");
	free(myself->backend_comp_md);
}

void CommandQueue_sync_barrier(CQueue_p myself, int task_complete_flag){
	cudaStream_t stream = *((cudaStream_t*) myself->backend_queue_ptr);
	cudaError_t err = cudaStreamSynchronize(stream);
	massert(cudaSuccess == err, "CommandQueue_sync_barrier - %s\n", cudaGetErrorString(err));
	if(task_complete_flag) myself->task_id.store(0);
}

void CommandQueue_add_host_func(CQueue_p myself){
	cudaStream_t stream = *((cudaStream_t*) myself->backend_queue_ptr);
	cudaError_t err = cudaLaunchHostFunc(stream, (cudaHostFn_t) myself->func_ptr, myself->data_ptr);
	massert(cudaSuccess == err, "CommandQueue::add_host_func - %s\n", cudaGetErrorString(err));
	//myself->task_id.store(0);
}

void CommandQueue_wait_for_event(CQueue_p myself){
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] |-----> CommandQueue_wait_for_event(Queue(%d), Event(%d))\n", myself->dev_id, myself->id, myself->wait_for_it->id);
#endif
	event_status evstat = myself->wait_for_it->query_status();
	if (evstat != CHECKED){
		if (evstat == UNRECORDED) {
		;
#ifdef UDDEBUG
			warning("CommandQueue_wait_for_event - UNRECORDED event\n");
#endif
		}
		else {
			cudaStream_t stream = *((cudaStream_t*) myself->backend_queue_ptr);
			cudaEvent_t cuda_event = *(cudaEvent_t*) myself->wait_for_it->event_backend_ptr;
			cudaError_t err = cudaStreamWaitEvent(stream, cuda_event, 0);
			massert(cudaSuccess == err, "[dev_id=%3d] CommandQueue_wait_for_event(Queue(%d), Event(%d)() - %s\n", 
			myself->dev_id, myself->id, myself->wait_for_it->id, cudaGetErrorString(err));
		}
	}
	//myself->task_id.store(0);
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] <-----| CommandQueue_wait_for_event(Queue(%d), Event(%d))\n", myself->dev_id, myself->id, myself->wait_for_it->id);
#endif
}

void CommandQueue_record_event(CQueue_p myself){
	if (myself->record_it->status != UNRECORDED && 
		myself->record_it->status != CHECKED){ // TODO: previous -> if (status != UNRECORDED)
		;
#ifdef UDEBUG
		warning("Event(%d,dev_id = %d)::record_to_queue(%d): Recording %s event\n", myself->record_it->id, 
		myself->record_it->dev_id, myself->dev_id, print_event_status(myself->record_it->status));
#endif
		if(myself->dev_id != myself->record_it->dev_id)
			error("Event(%d,dev_id = %d)::record_to_queue(%d): Recording %s event in illegal dev\n",
				myself->record_it->id, myself->record_it->dev_id, 
				myself->dev_id, print_event_status(myself->record_it->status));

	}
	//TODO: previous -> if (status == UNRECORDED)
	else if (myself->record_it->status == UNRECORDED || myself->record_it->status == CHECKED){ 
		if(myself->dev_id > -1){
			/// TODO: This used to be an error, but with soft reset it was problematic...is it ok?
			; //warning("Event(%d,dev_id = %d)::record_to_queue(%d) - UNRECORDED event suspicious dev_id\n",
			  //	id, dev_id, Rr->dev_id);
		}
		myself->record_it->dev_id = myself->dev_id;
		cudaError_t err = cudaEventCreate(( cudaEvent_t*) myself->record_it->event_backend_ptr);
		massert(cudaSuccess == err, "Event(%d,dev_id = %d)::record_to_queue(%d): - cudaEventCreate - %s\n",
			myself->record_it->id, myself->record_it->dev_id, myself->dev_id, cudaGetErrorString(err));
	}
	//int curr_dev_id = CHLGetDevice();
	//if (myself->dev_id!=curr_dev_id) warning("Event(%d,dev_id = %d)::record_to_queue(%d): Wrong current device %d\n",  
	//	myself->record_it->id, myself->record_it->dev_id, myself->dev_id, curr_dev_id);
	cudaEvent_t cuda_event= *(cudaEvent_t*) myself->record_it->event_backend_ptr;
	cudaStream_t stream = *((cudaStream_t*) myself->backend_queue_ptr);
	cudaError_t err = cudaEventRecord(cuda_event, stream);
	myself->record_it->status = RECORDED;
	massert(cudaSuccess == err, "Event(%d,dev_id = %d)::record_to_queue(%d, dev_id=%d) - cudaEventRecord - %s\n",  
		myself->record_it->id, myself->record_it->dev_id, myself->id, myself->dev_id, cudaGetErrorString(err));
	//myself->task_id.store(0);
}

void CommandQueue_run_operation(CQueue_p myself){
	if  (!(myself->type > 0)) error("CommandQueue_run_operation(): Invoked on communication queue\n");
	const char* opname = (const char*) myself->func_ptr;
	if(myself->target_dev_id == -1 || myself->target_dev_id >= CHL_WORKERS){
		if (!strcmp(opname, "Dgemm")) myself->func_ptr = (void*) &cblas_wrap_dgemm;
		if (!strcmp(opname, "Daxpy")) myself->func_ptr = (void*) &cblas_wrap_daxpy;
		if (!strcmp(opname, "Daxpby")) myself->func_ptr = (void*) &cblas_wrap_daxpby;
		if (!strcmp(opname, "Dslaxpby")) myself->func_ptr = (void*) &custom_cpu_wrap_dslaxpby;
		CommandQueue_add_host_func(myself);
	}
	else{
		//TODO: Why is this needed here, I wonder!
		CHLSelectDevice(myself->target_dev_id);
		if (!strcmp(opname, "Dgemm")) cublas_wrap_dgemm(myself->data_ptr, myself);
		if (!strcmp(opname, "Daxpy")) cublas_wrap_daxpy(myself->data_ptr, myself);
		if (!strcmp(opname, "Daxpby")) cublas_wrap_daxpby(myself->data_ptr, myself);
		if (!strcmp(opname, "Dslaxpby")) custom_gpu_wrap_dslaxpby(myself->data_ptr, myself);
	}
	//myself->task_id.store(0);
}
void CommandQueue_delete(CQueue_p myself){
	CommandQueue_sync_barrier(myself, 0);
	if(COMMUNICATION == myself->type) CommandQueue_comm_delete(myself);
	else if(COMPUTATION == myself->type) CommandQueue_comp_delete(myself);
	else if (MIXED == myself->type){
		CommandQueue_comm_delete(myself);
		CommandQueue_comp_delete(myself);
	}
	cudaStream_t stream = *((cudaStream_t*) myself->backend_queue_ptr);
	cudaError_t err = cudaStreamDestroy(stream);
	massert(cudaSuccess == err, "CommandQueue_delete - cudaStreamDestroy failed: %s\n", cudaGetErrorString(err));
	free(myself->backend_queue_ptr);
	//myself->task_id.store(0);
}

#ifdef TTEST /// C programmers hate him
int n_trans_ctr = 0;
long long n_bytes[100000] = {0};
int n_locs[10000][2];
double n_timers[100000][3];
int n_timer_ctr[64][64] = {{0}};
double n_link_gbytes_s[64][64] = {{0}};
int my_log_lock = 0; /// This might slow down things, but it is needed. 

void n_reseTTEST(){
	for(int k = 0; k < n_trans_ctr; k++){
		n_bytes[k] = 0;
		for(int m = 0; m < 3; m++) n_timers[k][m] = 0;
	}
	n_trans_ctr = 0;
	for (int d1 = 0; d1 < 64; d1++)
		for (int d2 = 0; d2 < 64; d2++){
			n_timer_ctr[d1][d2] = 0; 
			n_link_gbytes_s[d1][d2] = 0; 
		}
}

void n_HopMemcpyPrint(){
	lprintf(0,"\n Tranfers Full:\n");
	FILE* fp = fopen("temp_n_trans.log", "w+");
	for(int k = 0; k < n_trans_ctr; k++){
		int src = n_locs[k][0], dest = n_locs[k][1];
		n_timer_ctr[dest][src]++;
		double time = (n_timers[k][2] - n_timers[k][1]), pipe_time = (n_timers[k][2] - n_timers[k][0]);
		n_link_gbytes_s[dest][src]+=Gval_per_s(n_bytes[k], time);
		//lprintf(0, "Normal 2D Trasfer %d->%d : total_t = %lf ms ( %.3lf Gb/s ), pipelined_t = %lf ms ( %.3lf Gb/s )\n", 
		//	src, dest, 1000*time, Gval_per_s(n_bytes[k], time), 1000*pipe_time, Gval_per_s(n_bytes[k], pipe_time));
		fprintf(fp, "%d,%d,[ %d %d ],%ld,%lf,%lf,%lf\n", src, dest, src, dest, n_bytes[k], n_timers[k][0], n_timers[k][1], n_timers[k][2]);
	}
		
	lprintf(0,"\n Full Tranfer Map:\n   |");
	for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
		lprintf(0, "  %s  |", mem_name(d2));
	lprintf(0, "\n   |");
	for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
		lprintf(0, "-------");
	lprintf(0, "\n");
	for (int d1 = 0; d1 < CHL_MEMLOCS; d1++){
		lprintf(0, "%s | ", mem_name(d1));
		for (int d2 = 0; d2 < CHL_MEMLOCS; d2++){
			lprintf(0, "%4d | ", n_timer_ctr[d1][d2]);
		}
		lprintf(0, "\n");
	}

	lprintf(0,"\n Full Tranfer Map Achieved Bandwidths (GB/s):\n   |");
	for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
		lprintf(0, "  %s   |", mem_name(d2));
	lprintf(0, "\n   |");
	for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
		lprintf(0, "--------");
	lprintf(0, "\n");
	for (int d1 = 0; d1 < CHL_MEMLOCS; d1++){
		lprintf(0, "%s | ", mem_name(d1));
		for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
			if (n_timer_ctr[d1][d2]) lprintf(0, "%.2lf | ", n_link_gbytes_s[d1][d2]/n_timer_ctr[d1][d2]);
			else lprintf(0, "  -   | ");
		lprintf(0, "\n");
	}
	fclose(fp);
	n_reseTTEST();
}
#endif 

void CommandQueue_memcpy2DAsync(CQueue_p myself){
#ifdef TTEST
	if(myself->comm_wrap_ptr->log_flag){
		while(__sync_lock_test_and_set(&my_log_lock, 1));
		if (n_trans_ctr > 100000) 
			error("CommandQueue_memcpy2DAsync(dest = %d, src = %d) exeeded 100000 transfers in TTEST Mode\n", 
				myself->comm_wrap_ptr->loc_dest, myself->comm_wrap_ptr->loc_src);
		if(!n_trans_ctr) n_reseTTEST();
		n_bytes[n_trans_ctr] = myself->comm_wrap_ptr->rows*myself->comm_wrap_ptr->cols*myself->comm_wrap_ptr->elemSize;
		n_locs[n_trans_ctr][0] = myself->comm_wrap_ptr->loc_src;
		n_locs[n_trans_ctr][1] = myself->comm_wrap_ptr->loc_dest;
	}
#endif
	cudaStream_t stream = *((cudaStream_t*) myself->backend_queue_ptr);

	massert(myself->comm_wrap_ptr->loc_dest >= -1 && myself->comm_wrap_ptr->loc_dest < CHL_MEMLOCS, 
		"CommandQueue(%d)_memcpy2DAsync: Invalid destination device: %d\n", myself->id, myself->comm_wrap_ptr->loc_dest);
	massert(myself->comm_wrap_ptr->loc_src >= -1 && myself->comm_wrap_ptr->loc_src < CHL_MEMLOCS, 
		"CommandQueue(%d)_memcpy2DAsync: Invalid source device: %d\n", myself->id, myself->comm_wrap_ptr->loc_src);

	enum cudaMemcpyKind kind;
	if ((myself->comm_wrap_ptr->loc_src >= CHL_WORKERS  || myself->comm_wrap_ptr->loc_src == -1 ) 
		&& (myself->comm_wrap_ptr->loc_dest >= CHL_WORKERS   || myself->comm_wrap_ptr->loc_dest == -1 )) 
		kind = cudaMemcpyHostToHost;
	else if (myself->comm_wrap_ptr->loc_dest >= CHL_WORKERS  || myself->comm_wrap_ptr->loc_dest == -1) kind = cudaMemcpyDeviceToHost;
	else if (myself->comm_wrap_ptr->loc_src >= CHL_WORKERS || myself->comm_wrap_ptr->loc_src == -1) kind = cudaMemcpyHostToDevice;
	else kind = cudaMemcpyDeviceToDevice;

	if (myself->comm_wrap_ptr->loc_src == myself->comm_wrap_ptr->loc_dest) warning("CommandQueue(%d)_memcpy2DAsync(dest=%p, ldest =%zu, src=%p,"
	"ldsrc = %zu, rows=%zu, cols=%zu, elemSize =%d, loc_dest=%d, loc_src=%d): Source location matches destination\n",
	myself->id, myself->comm_wrap_ptr->dest, myself->comm_wrap_ptr->ldest, myself->comm_wrap_ptr->src, myself->comm_wrap_ptr->ldsrc, myself->comm_wrap_ptr->rows, 
	myself->comm_wrap_ptr->cols, myself->comm_wrap_ptr->elemSize, myself->comm_wrap_ptr->loc_dest, myself->comm_wrap_ptr->loc_src);
#ifdef TTEST
	if(myself->comm_wrap_ptr->log_flag){
		CHLSetTimerAsync(&(n_timers[n_trans_ctr][0]));
		myself->add_host_func((void*)&CHLSetTimerAsync, 
			(void*) &(n_timers[n_trans_ctr][1]));
	}
#endif	
	//if (loc_src == -1 && loc_dest >=0) massert(CUBLAS_STATUS_SUCCESS == cublasSetMatrixAsync(rows, cols, elemSize, src, ldsrc, dest, ldest, stream), "CHLMemcpy2DAsync: cublasSetMatrixAsync failed\n");
	//else if (loc_src >=0 && loc_dest == -1) massert(CUBLAS_STATUS_SUCCESS == cublasGetMatrixAsync(rows, cols, elemSize, src, ldsrc, dest, ldest, stream),  "CHLMemcpy2DAsync: cublasGetMatrixAsync failed");
	//else 
	massert(cudaSuccess == cudaMemcpy2DAsync(myself->comm_wrap_ptr->dest, myself->comm_wrap_ptr->ldest*myself->comm_wrap_ptr->elemSize, 
	myself->comm_wrap_ptr->src, myself->comm_wrap_ptr->ldsrc*myself->comm_wrap_ptr->elemSize,
		myself->comm_wrap_ptr->rows*myself->comm_wrap_ptr->elemSize, myself->comm_wrap_ptr->cols, kind, stream), 
		"CHLMemcpy2DAsync(dest=%p, ldest =%zu, src=%p, ldsrc = %zu\nrows = %zu, cols = %zu, elemsize = %d,"
		"loc_dest = %d, loc_src = %d): cudaMemcpy2DAsync failed\n",
			myself->comm_wrap_ptr->dest, myself->comm_wrap_ptr->ldest, myself->comm_wrap_ptr->src, myself->comm_wrap_ptr->ldsrc, myself->comm_wrap_ptr->rows, 
			myself->comm_wrap_ptr->cols, myself->comm_wrap_ptr->elemSize, myself->comm_wrap_ptr->loc_dest, myself->comm_wrap_ptr->loc_src);
#ifdef TTEST
	if(myself->comm_wrap_ptr->log_flag){
		myself->add_host_func((void*)&CHLSetTimerAsync, 
			(void*) &(n_timers[n_trans_ctr][2]));
		n_trans_ctr++;
		__sync_lock_release(&my_log_lock);
	}
#endif		
	//myself->task_id.store(0);
}

void* CommandQueue_worker(void* wrapped_CQueue_p){
	CQueue_p myself = (CQueue_p) wrapped_CQueue_p;
	int cpu_aff = CHL_HWNUMA_AT_MEMLOC[CHL_WORKER_CLOSE_TO_MEMLOC[myself->dev_id]];//myself->dev_id; //
	cpu_set_t cpuset;
	CPU_ZERO(&cpuset);
	// Set thread affinity to the first core in its 'friend' CPU.
	// For now, always use core 0 from each CPU for queue worker threads.
	int core_aff = cpu_aff*HW_THREADS/NUMA_HW_NUM + 0;
	for(int idx = 0; idx < HW_THREADS/NUMA_HW_NUM; idx++){
		core_aff = cpu_aff*HW_THREADS/NUMA_HW_NUM + idx;
		CPU_SET(core_aff, &cpuset);
	}
	pthread_t curr_thread = pthread_self();
	pthread_setaffinity_np(curr_thread, sizeof(cpu_set_t), &cpuset);
	CHLSelectDevice(myself->dev_id);
	myself->backend_queue_ptr = malloc(sizeof(cudaStream_t));
	cudaError_t err = cudaStreamCreateWithFlags((cudaStream_t*) myself->backend_queue_ptr, cudaStreamNonBlocking);
	massert(cudaSuccess == err, "CommandQueue_worker(%d) - %s\n", myself->dev_id, cudaGetErrorString(err));
	cudaStream_t stream = *((cudaStream_t*) myself->backend_queue_ptr);
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] CommandQueue_worker(%d) - Initialized stream at CPU %d\n", myself->dev_id, myself->id, cpu_aff); 
#endif
	if(COMMUNICATION == myself->type) CommandQueue_comm_init(myself);
	else if(COMPUTATION == myself->type) CommandQueue_comp_init(myself);
	else if (MIXED == myself->type){
		CommandQueue_comm_init(myself);
		CommandQueue_comp_init(myself);
	}
	myself->task_id.store(0);
	pthread_mutex_lock(myself->backend_mutex_ptr);
	while (1){
		//usleep(300);
		int local_task_id = myself->task_id.load();
		if (local_task_id == 0){
			//continue;
			pthread_cond_wait(myself->backend_cond_ptr, myself->backend_mutex_ptr);
		}
		if (local_task_id == 1) error("CommandQueue_worker(%d):"
			"Contructor called in again on existing queue...bug\n", myself->dev_id);
		else if(local_task_id == 2) CommandQueue_sync_barrier(myself, 0);
		else if(local_task_id == 3) CommandQueue_add_host_func(myself);
		else if(local_task_id == 4) CommandQueue_wait_for_event(myself);
		else if(local_task_id == 5) CommandQueue_record_event(myself);
		else if(local_task_id == 6) CommandQueue_run_operation(myself);
		else if(local_task_id == 7) CommandQueue_memcpy2DAsync(myself);
		else if(local_task_id == -42){
			CommandQueue_delete(myself);
			break;
		}
		else if(local_task_id != 0) error("dev_id=%3d] CommandQueue_worker(%d): was singaled with local_task_id = %d\n", 
			myself->dev_id, myself->id, local_task_id);
		myself->task_id.store(0);
	}
	pthread_mutex_unlock(myself->backend_mutex_ptr);
	fprintf(stderr, "[dev_id=%3d] CommandQueue_worker(%d): I am about to return\n", myself->dev_id, myself->id);
	return NULL;
}

CommandQueue::CommandQueue(int dev_id_in, CQueue_type type_in)
{
	id = Queue_num_loc[dev_id_in]++;
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] |-----> CommandQueue(%d)::CommandQueue()\n", dev_id_in, id);
#endif
	dev_id = dev_id_in;
	type = type_in;
	backend_queue_ptr = backend_comp_md = backend_comm_md = NULL;
	queue_ETA = 0;
	// Spawn pthread
	task_id.store(1);
	comm_wrap_ptr = (CO_p) calloc(1, sizeof(struct comm_wrap));
#ifdef ENABLE_PTHREAD_BACKQUEUES
	pthread_attr_t attr;
	if (pthread_attr_init(&attr) != 0) error("pthread_attr_init failed for CommandQueue(%d)\n", id);
	backend_mutex_ptr = (pthread_mutex_t*) calloc(1, sizeof(pthread_mutex_t));
	if (pthread_mutex_init(backend_mutex_ptr, NULL) != 0) 
		error("pthread_mutex_init failed for CommandQueue(%d)\n", id);
	queue_mutex_ptr = (pthread_mutex_t*) calloc(1, sizeof(pthread_mutex_t));
	if (pthread_mutex_init(queue_mutex_ptr, NULL) != 0) 
		error("pthread_mutex_init failed for CommandQueue(%d)\n", id);
	backend_cond_ptr = (pthread_cond_t*) calloc(1, sizeof(pthread_cond_t));
	if(pthread_cond_init(backend_cond_ptr, NULL) != 0)
		error("pthread_cond_init failed for CommandQueue(%d)\n", id);
	backend_worker_ptr = (pthread_t*) calloc(1, sizeof(pthread_t));
	if(pthread_create(backend_worker_ptr, &attr, CommandQueue_worker, (void*)this) != 0)
		error("pthread_create failed for CommandQueue(%d)\n", id);
	// Wait for its completion and return
	while(task_id.load()){ ; } // block  
#else
	CHLSelectDevice(dev_id);
	backend_queue_ptr = malloc(sizeof(cudaStream_t));
	cudaError_t err = cudaStreamCreateWithFlags((cudaStream_t*) backend_queue_ptr, cudaStreamNonBlocking);
	massert(cudaSuccess == err, "CommandQueue_worker(%d) - %s\n", dev_id, cudaGetErrorString(err));
	cudaStream_t stream = *((cudaStream_t*) backend_queue_ptr);
	if(COMMUNICATION == type) CommandQueue_comm_init(this);
	else if(COMPUTATION == type) CommandQueue_comp_init(this);
	else if (MIXED == type){
		CommandQueue_comm_init(this);
		CommandQueue_comp_init(this);
	}
#endif
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] <-----| CommandQueue(%d)::CommandQueue()\n", dev_id, id);
#endif
}

CommandQueue::~CommandQueue()
{
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] |-----> CommandQueue(%d)::~CommandQueue()\n", dev_id, id);
#endif
#ifdef ENABLE_PTHREAD_BACKQUEUES
	pthread_mutex_lock(queue_mutex_ptr);
	pthread_mutex_lock(backend_mutex_ptr);
	task_id.store(-42);
	pthread_cond_signal(backend_cond_ptr);
	pthread_mutex_unlock(backend_mutex_ptr);
	// Join pthread at backend_worker_ptr
	pthread_join(*backend_worker_ptr, NULL);
	pthread_mutex_destroy(backend_mutex_ptr);
	pthread_cond_destroy(backend_cond_ptr);
	pthread_mutex_unlock(queue_mutex_ptr);
	pthread_mutex_destroy(queue_mutex_ptr);
	free(backend_mutex_ptr);
	free(backend_cond_ptr);	
	free(queue_mutex_ptr);
	free(backend_worker_ptr);
#else
	CHLSelectDevice(dev_id);
	CommandQueue_delete(this);
#endif
	free(comm_wrap_ptr);

#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] <-----| CommandQueue(%d)::~CommandQueue()\n", dev_id, id);
#endif
	return;
}

void CommandQueue::sync_barrier()
{
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] |-----> CommandQueue(%d)::sync_barrier()\n", dev_id, id);
#endif
#ifdef ENABLE_PTHREAD_BACKQUEUES
	pthread_mutex_lock(queue_mutex_ptr);
	pthread_mutex_lock(backend_mutex_ptr);
	task_id.store(2); 
	pthread_cond_signal(backend_cond_ptr);
	pthread_mutex_unlock(backend_mutex_ptr);
	while(task_id.load()){ ; } // block 
	pthread_mutex_unlock(queue_mutex_ptr);
#else
	//CHLSelectDevice(dev_id);
	CommandQueue_sync_barrier(this, 0);
#endif
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] <-----| CommandQueue(%d)::sync_barrier()\n", dev_id, id);
#endif
}

void CommandQueue::add_host_func(void* func, void* data){
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] |-----> CommandQueue(%d)::add_host_func()\n", dev_id, id);
#endif
#ifdef ENABLE_PTHREAD_BACKQUEUES
	pthread_mutex_lock(queue_mutex_ptr);
	pthread_mutex_lock(backend_mutex_ptr);
	func_ptr = func;
	data_ptr = data;
	task_id.store(3); 
	pthread_cond_signal(backend_cond_ptr);
	pthread_mutex_unlock(backend_mutex_ptr);
	while(task_id.load()){ ; } // block 
	pthread_mutex_unlock(queue_mutex_ptr);
#else
	//CHLSelectDevice(dev_id);
	func_ptr = func;
	data_ptr = data;
	CommandQueue_add_host_func(this);
#endif
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] <-----| CommandQueue(%d)::add_host_func()\n", dev_id, id);
#endif
}

void CommandQueue::wait_for_event(Event_p Wevent)
{
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] |-----> CommandQueue(%d)::wait_for_event(Event(%d))\n", dev_id, id, Wevent->id);
#endif
#ifdef ENABLE_PTHREAD_BACKQUEUES
	pthread_mutex_lock(queue_mutex_ptr);
	pthread_mutex_lock(backend_mutex_ptr);
	wait_for_it = Wevent;
	task_id.store(4); 
	pthread_cond_signal(backend_cond_ptr);
	pthread_mutex_unlock(backend_mutex_ptr);
	while(task_id.load()){ ; } // block 
	pthread_mutex_unlock(queue_mutex_ptr);
#else
	wait_for_it = Wevent;
	CommandQueue_wait_for_event(this);
#endif
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] <-----| CommandQueue(%d)::wait_for_event(Event(%d))\n", dev_id, id, Wevent->id);
#endif
	return;
}

void CommandQueue::record_event(Event_p Wevent)
{
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] |-----> CommandQueue(%d)::record_event(Event(%d))\n", dev_id, id, Wevent->id);
#endif
#ifdef ENABLE_PTHREAD_BACKQUEUES
	pthread_mutex_lock(queue_mutex_ptr);
	pthread_mutex_lock(backend_mutex_ptr);
	record_it = Wevent;
	task_id.store(5); 
	pthread_cond_signal(backend_cond_ptr);
	pthread_mutex_unlock(backend_mutex_ptr);
	while(task_id.load()){ ; } // block 
	pthread_mutex_unlock(queue_mutex_ptr);
#else
	CHLSelectDevice(dev_id);
	record_it = Wevent;
	CommandQueue_record_event(this);
#endif
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] <-----| CommandQueue(%d)::record_event(Event(%d))\n", dev_id, id, Wevent->id);
#endif
	return;
}

void CommandQueue::run_operation(void* backend_data, const char* opname, int target_dev_id_in){
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] |-----> CommandQueue(%d)::run_operation()\n", dev_id, id);
#endif
#ifdef ENABLE_PTHREAD_BACKQUEUES
	pthread_mutex_lock(queue_mutex_ptr);
	pthread_mutex_lock(backend_mutex_ptr);
	func_ptr = (void*) opname;
	data_ptr = backend_data;
	target_dev_id = target_dev_id_in;
	task_id.store(6); 
	pthread_cond_signal(backend_cond_ptr);
	pthread_mutex_unlock(backend_mutex_ptr);
	while(task_id.load()){ ; } // block 
	pthread_mutex_unlock(queue_mutex_ptr);
#else
	func_ptr = (void*) opname;
	data_ptr = backend_data;
	target_dev_id = target_dev_id_in;
	CommandQueue_run_operation(this);
#endif
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] <-----| CommandQueue(%d)::run_operation()\n", dev_id, id);
#endif
}

void CommandQueue::memcpy2DAsync(void* dest, long int ldest, void* src, long int ldsrc, long int rows, 
	long int cols, int elemSize, int loc_dest, int loc_src, int log_flag){
#ifdef DDEBUG
	fprintf(stderr, "CommandQueue(%d):memcpy2DAsync(dest=%p, ldest =%zu, src=%p, ldsrc = %zu, rows = %zu, cols = %zu, elemsize = %d, loc_dest = %d, loc_src = %d)\n",
		id, dest, ldest, src, ldsrc, rows, cols, elemSize, loc_dest, loc_src);
#endif
#ifdef ENABLE_PTHREAD_BACKQUEUES
	pthread_mutex_lock(queue_mutex_ptr);
	pthread_mutex_lock(backend_mutex_ptr);
	comm_wrap_ptr->dest = dest;
	comm_wrap_ptr->src = src;
	comm_wrap_ptr->ldest = ldest;
	comm_wrap_ptr->ldsrc = ldsrc;
	comm_wrap_ptr->rows = rows;
	comm_wrap_ptr->cols = cols;
	comm_wrap_ptr->elemSize = elemSize;
	comm_wrap_ptr->loc_dest = loc_dest;
	comm_wrap_ptr->loc_src = loc_src;
	comm_wrap_ptr->log_flag = log_flag;
	task_id.store(7); 
	pthread_cond_signal(backend_cond_ptr);
	pthread_mutex_unlock(backend_mutex_ptr);
	while(task_id.load()){ ; } // block 
	pthread_mutex_unlock(queue_mutex_ptr);
#else
	comm_wrap_ptr->dest = dest;
	comm_wrap_ptr->src = src;
	comm_wrap_ptr->ldest = ldest;
	comm_wrap_ptr->ldsrc = ldsrc;
	comm_wrap_ptr->rows = rows;
	comm_wrap_ptr->cols = cols;
	comm_wrap_ptr->elemSize = elemSize;
	comm_wrap_ptr->loc_dest = loc_dest;
	comm_wrap_ptr->loc_src = loc_src;
	comm_wrap_ptr->log_flag = log_flag;	
	CommandQueue_memcpy2DAsync(this);
#endif
}

/*****************************************************/
/// PARALia 2.0 - timed queues

void CommandQueue::ETA_add_task(long double task_fire_t, long double task_duration){
	queue_ETA = fmax(queue_ETA, task_fire_t) + task_duration;
}

void CommandQueue::ETA_set(long double new_ETA){
	queue_ETA = new_ETA; 
}

long double CommandQueue::ETA_get(){
	return queue_ETA;
}

/*****************************************************/
/// Event class functions.
Event::Event(int dev_id_in)
{
#ifdef UDDEBUG
	fprintf(stderr, "[dev_id=%3d] |-----> Event(%d)::Event()\n", dev_id_in, Event_num_loc[dev_id_in]);
#endif
	event_backend_ptr = malloc(sizeof(cudaEvent_t));
	id = Event_num_loc[dev_id_in];
	Event_num_loc[dev_id_in]++;
	dev_id = dev_id_in - 42;

	status = UNRECORDED;
#ifdef UDDEBUG
	fprintf(stderr, "[dev_id=%3d] <-----| Event(%d)::Event()\n", dev_id, id);
#endif
}

Event::~Event()
{
#ifdef UDDEBUG
	fprintf(stderr, "[dev_id=%3d] |-----> Event(%d)::~Event()\n", dev_id, id);
#endif
	sync_barrier();
	if (dev_id < -1) Event_num_loc[dev_id+42]--;
	else{
			Event_num_loc[dev_id]--;
			cudaError_t err = cudaEventDestroy(*(( cudaEvent_t*) event_backend_ptr));
			massert(cudaSuccess == err, "Event(%d)::~Event() - %s\n", id, cudaGetErrorString(err));
	}

	free(event_backend_ptr);
#ifdef UDDEBUG
	fprintf(stderr, "[dev_id=%3d] <-----| Event(%d)::~Event()\n", dev_id, id);
#endif
}

void Event::sync_barrier()
{
#ifdef UDDEBUG
	fprintf(stderr, "[dev_id=%3d] |-----> Event(%d)::sync_barrier()\n", dev_id, id);
#endif
	if (status != CHECKED){
		if (status == UNRECORDED){;
#ifdef UDEBUG
			warning("[dev_id=%3d] |-----> Event(%d)::sync_barrier() - Tried to sync unrecorded event\n", dev_id, id);
#endif
		}
		else{
			cudaEvent_t cuda_event= *(cudaEvent_t*) event_backend_ptr;
			cudaError_t err = cudaEventSynchronize(cuda_event);
			if (status == RECORDED) status = CHECKED;
			massert(cudaSuccess == err, "Event::sync_barrier() - %s\n", cudaGetErrorString(err));
		}
	}
#ifdef UDDEBUG
	fprintf(stderr, "[dev_id=%3d] <-----| Event(%d)::sync_barrier()\n", dev_id, id);
#endif
	return;
}

void Event::record_to_queue(CQueue_p Rr){
	if (Rr == NULL){
#ifdef UDDEBUG
	fprintf(stderr, "[dev_id=%3d] <-----> Event(%d)::record_to_queue(NULL)\n", dev_id, id);
#endif
		status = CHECKED;
		return;
	}
#ifdef UDDEBUG
	fprintf(stderr, "[dev_id=%3d] |-----> Event(%d)::record_to_queue(Queue(dev_id=%d))\n", dev_id, id, Rr->dev_id);
#endif
	Rr->record_event(this);
#ifdef UDDEBUG
	fprintf(stderr, "[dev_id=%3d] <-----| Event(%d)::record_to_queue(Queue(dev_id=%d))\n", dev_id, id, Rr->dev_id);
#endif
}

event_status Event::query_status(){
#ifdef UDDEBUG
	fprintf(stderr, "[dev_id=%3d] |-----> Event(%d)::query_status()\n", dev_id, id);
#endif
	enum event_status local_status = status;
	if (local_status != CHECKED){
		if (local_status == UNRECORDED) return UNRECORDED;
		cudaEvent_t cuda_event= *(cudaEvent_t*) event_backend_ptr;
		cudaError_t err = cudaEventQuery(cuda_event);

		if (err == cudaSuccess && (local_status == UNRECORDED ||  local_status == COMPLETE));
		else if (err == cudaSuccess && local_status == RECORDED) local_status = status = COMPLETE;
		else if (err == cudaErrorNotReady && local_status == RECORDED);
		else if (err == cudaErrorNotReady && local_status == UNRECORDED){
#ifdef UDEBUG
			// this should not happen in a healthy locked update scenario.
			warning("Event::query_status(): cudaErrorNotReady with status == UNRECORDED should not happen\n");
#endif
			local_status = status = RECORDED;
		}
		else if (err == cudaSuccess &&  local_status == CHECKED){
			;
			// TODO: This should not happen in a healthy locked update scenario.
			// But it does since no locking yet. Not sure of its effects.
#ifdef UDEBUG
			warning("[dev_id=%3d] |-----> Event(%d)::query_status(): cudaSuccess with local_status == CHECKED should not happen\n", dev_id, id);
#endif
		}
		else error("[dev_id=%3d] |-----> Event(%d)::query_status() - %s, local_status=%s, status = %s\n", dev_id, id,
		cudaGetErrorString(err), print_event_status(local_status), print_event_status(status));
	}
#ifdef UDDEBUG
	fprintf(stderr, "[dev_id=%3d] <-----| Event(%d)::query_status() = %s\n", dev_id, id, print_event_status(status));
#endif
	return local_status;
}

void Event::checked(){
#ifdef UDDEBUG
	fprintf(stderr, "[dev_id=%3d] |-----> Event(%d)::checked()\n", dev_id, id);
#endif
	if (status == COMPLETE) status = CHECKED;
	else error("Event::checked(): error event was %s,  not COMPLETE()\n", print_event_status(status));
#ifdef UDDEBUG
	fprintf(stderr, "[dev_id=%3d] <-----| Event(%d)::checked()\n", dev_id, id);
#endif
}

void Event::soft_reset(){
#ifdef UDDEBUG
	fprintf(stderr, "[dev_id=%3d] |-----> Event(%d)::soft_reset()\n", dev_id, id);
#endif
	// sync_barrier();
	// event_status prev_status = status;
	status = UNRECORDED;
	if(dev_id >= -1){
		dev_id = dev_id - 42;
		event_backend_ptr = malloc(sizeof(cudaEvent_t));
	}
#ifdef UDDEBUG
	fprintf(stderr, "[dev_id=%3d] <-----| Event(%d)::soft_reset()\n", dev_id, id);
#endif
}

void Event::reset(){
#ifdef UDDEBUG
	fprintf(stderr, "[dev_id=%3d] Event(%d)::reset()\n", dev_id, id);
#endif
	sync_barrier();
	event_status prev_status = status;
	status = UNRECORDED;
	if(dev_id >= -1){
		dev_id = dev_id - 42;
		cudaError_t err = cudaEventDestroy(*(( cudaEvent_t*) event_backend_ptr));
		massert(cudaSuccess == err, "[dev_id=%3d] Event(%d)::reset - %s\n", dev_id + 42, id, cudaGetErrorString(err));
	}
#ifdef UDDEBUG
	fprintf(stderr, "[dev_id=%3d] <-----| Event(%d)::reset()\n", dev_id, id);
#endif
}

/*****************************************************/
/// Event-based timer class functions

Event_timer::Event_timer(int dev_id) {
  Event_start = new Event(dev_id);
  Event_stop = new Event(dev_id);
  time_ms = 0;
}

void Event_timer::start_point(CQueue_p start_queue)
{
	Event_start->record_to_queue(start_queue);
}

void Event_timer::stop_point(CQueue_p stop_queue)
{
	Event_stop->record_to_queue(stop_queue);
}

double Event_timer::sync_get_time()
{
	float temp_t;
	if(Event_start->query_status() != UNRECORDED){
		Event_start->sync_barrier();
		if(Event_stop->query_status() != UNRECORDED) Event_stop->sync_barrier();
		else error("Event_timer::sync_get_time: Event_start is %s but Event_stop still UNRECORDED\n",
			print_event_status(Event_start->query_status()));
		cudaEvent_t cuda_event_start = *(cudaEvent_t*) Event_start->event_backend_ptr;
		cudaEvent_t cuda_event_stop = *(cudaEvent_t*) Event_stop->event_backend_ptr;
		cudaEventElapsedTime(&temp_t, cuda_event_start, cuda_event_stop);
	}
	else temp_t = 0;
	time_ms = (double) temp_t;
	return time_ms;
}

/*****************************************************/
