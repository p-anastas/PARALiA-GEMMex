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
		default:
			error("print_event_status: Unknown state\n");
	}
}

/*****************************************************/
/// Command queue class functions

void CommandQueue_comp_init(CQueue_p myself){
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] |-----> CommandQueue_comp_init(Queue(%d))\n", dev_id, id);
#endif
	cudaStream_t stream = *((cudaStream_t*) backend_queue_ptr);
	backend_comp_md = malloc(sizeof(cublasHandle_t));
	massert(CUBLAS_STATUS_SUCCESS == cublasCreate((cublasHandle_t*) backend_comp_md),
		"CommandQueue_comp_init(%d): cublasCreate failed\n", dev_id);
	massert(CUBLAS_STATUS_SUCCESS == cublasSetStream(*((cublasHandle_t*) backend_comp_md), stream),
		"CommandQueue_comp_init(%d): cublasSetStream failed\n", dev_id);
	if(WS_SZ >=0){
	void* local_ws = NULL; if(WS_SZ) cudaMalloc(&local_ws, WS_SZ); 
	massert(CUBLAS_STATUS_SUCCESS == cublasSetWorkspace(*((cublasHandle_t*) 
		backend_comp_md), local_ws, WS_SZ), 
		"CommandQueue_comp_init(%d): cublasSetWorkspace failed\n", dev_id);
	}
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] <-----| CommandQueue_comp_init(Queue(%d))\n", dev_id, id);
#endif
}

CommandQueue::CommandQueue(int dev_id_in, CQueue_type type_in)
{
	id = Queue_num_loc[dev_id_in]++;
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] |-----> CommandQueue(%d)::CommandQueue()\n", dev_id_in, id);
#endif
	dev_id = dev_id_in;
	type = type_in;
	backend_queue_ptr = backend_comp_md = NULL;
	queue_ETA = 0;
	CHLSelectDevice(dev_id);
	// Create the queue backend stream
	backend_queue_ptr = malloc(sizeof(cudaStream_t));
	cudaError_t err = cudaStreamCreateWithFlags((cudaStream_t*) backend_queue_ptr, cudaStreamNonBlocking);
#ifndef PRODUCTION
	massert(cudaSuccess == err, "CommandQueue(%d) - %s\n", dev_id, cudaGetErrorString(err));
#endif
	cudaStream_t stream = *((cudaStream_t*) backend_queue_ptr);
	// Also initialize computation backend if needed
	if(COMPUTATION == type || MIXED == type) CommandQueue_comp_init(this);

#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] <-----| CommandQueue(%d)::CommandQueue()\n", dev_id, id);
#endif
}

void CommandQueue_comp_delete(CQueue_p myself){
	cublasHandle_t handle = *((cublasHandle_t*) backend_comp_md);
	massert(CUBLAS_STATUS_SUCCESS == cublasDestroy(handle),
		"CommandQueue_comp_delete - cublasDestroy(handle) failed\n");
	free(backend_comp_md);
}

CommandQueue::~CommandQueue()
{
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] |-----> CommandQueue(%d)::~CommandQueue()\n", dev_id, id);
#endif
	CHLSelectDevice(dev_id);
	sync_barrier();
	if(COMPUTATION == type || MIXED == type) CommandQueue_comp_delete(this);
	cudaStream_t stream = *((cudaStream_t*) backend_queue_ptr);
	cudaError_t err = cudaStreamDestroy(stream);
#ifndef PRODUCTION
	massert(cudaSuccess == err, "CommandQueue(%d)::~CommandQueue - cudaStreamDestroy failed: %s\n", id, cudaGetErrorString(err));
#endif
	free(backend_queue_ptr);
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
	//CHLSelectDevice(dev_id);
	cudaStream_t stream = *((cudaStream_t*) backend_queue_ptr);
	cudaError_t err = cudaStreamSynchronize(stream);
#ifndef PRODUCTION
	massert(cudaSuccess == err, "CommandQueue(%d)::sync_barrier - %s\n", id, cudaGetErrorString(err));
#endif
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] <-----| CommandQueue(%d)::sync_barrier()\n", dev_id, id);
#endif
}

void CommandQueue::add_host_func(void* func, void* data){
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] |-----> CommandQueue(%d)::add_host_func()\n", dev_id, id);
#endif
	//CHLSelectDevice(dev_id);
	cudaStream_t stream = *((cudaStream_t*) backend_queue_ptr);
	cudaError_t err = cudaLaunchHostFunc(stream, (cudaHostFn_t) func, data);
#ifndef PRODUCTION
	massert(cudaSuccess == err, "CommandQueue::add_host_func - %s\n", cudaGetErrorString(err));
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
	event_status evstat = Wevent->query_status();
	if (evstat == COMPLETE) return;
	else if (evstat == UNRECORDED) {
#ifndef PRODUCTION
		warning("[dev_id=%3d] CommandQueue(%d)::wait_for_event(%d) - "
			"UNRECORDED event, looping until recorded\n",dev_id, id, Wevent->id);
#endif
		while (Wevent->query_status() == UNRECORDED);
	}
	else {
		cudaStream_t stream = *((cudaStream_t*) backend_queue_ptr);
		cudaEvent_t cuda_event = *(cudaEvent_t*) Wevent->event_backend_ptr;
		cudaError_t err = cudaStreamWaitEvent(stream, cuda_event, 0);
		massert(cudaSuccess == err, "[dev_id=%3d] CommandQueue_wait_for_event(Queue(%d), Event(%d)() - %s\n", 
		dev_id, id, Wevent->id, cudaGetErrorString(err));
	}
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
	CHLSelectDevice(dev_id);
#ifndef PRODUCTION
	event_status evstat = Wevent->query_status();
	if (evstat != UNRECORDED) error("Event(%d,dev_id = %d)::record_to_queue(%d): Recording %s event\n",
			Wevent->id, Wevent->dev_id, dev_id, print_event_status(evstat));
	else if(dev_id != Wevent->dev_id) error("[dev_id=%3d] CommandQueue(%d)::record_event(%d,dev_id = %d)"
		"Queue and event dev_id do not match\n", dev_id, id, Wevent->id, Wevent->dev_id);
#endif
	cudaEvent_t cuda_event= *(cudaEvent_t*) Wevent->event_backend_ptr;
	cudaStream_t stream = *((cudaStream_t*) backend_queue_ptr);
	cudaError_t err = cudaEventRecord(cuda_event, stream);
	Wevent->status = RECORDED;
#ifndef PRODUCTION
	massert(cudaSuccess == err, "[dev_id=%3d] CommandQueue(%d)::record_event(%d,dev_id = %d)- cudaEventRecord - %s\n",  
		id, dev_id, Wevent->id, Wevent->dev_id, cudaGetErrorString(err));
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
#ifndef PRODUCTION
	if  (type == COMMUNICATION) error("CommandQueue(%d)::run_operation(): Invoked on communication queue\n", id);
#endif
	//TODO: This seemed to be needed here for the GPU calls...
	CHLSelectDevice(target_dev_id_in);
	if(target_dev_id_in == -1 || target_dev_id_in >= CHL_WORKERS){
		if (!strcmp(opname, "Dgemm")) add_host_func( &cblas_wrap_dgemm, backend_data);
		else if (!strcmp(opname, "Daxpy")) add_host_func( &cblas_wrap_daxpy, backend_data);
		else if (!strcmp(opname, "Daxpby")) add_host_func( &cblas_wrap_daxpby, backend_data);
		else if (!strcmp(opname, "Dslaxpby")) add_host_func( &custom_cpu_wrap_dslaxpby, backend_data);
		else error("CommandQueue(%d)::run_operation(): unsupported opname = %s\n", id, opname);
	}
	else{
		if (!strcmp(opname, "Dgemm")) cublas_wrap_dgemm(backend_data, this);
		else if (!strcmp(opname, "Daxpy")) cublas_wrap_daxpy(backend_data, this);
		else if (!strcmp(opname, "Daxpby")) cublas_wrap_daxpby(backend_data, this);
		else if (!strcmp(opname, "Dslaxpby")) custom_gpu_wrap_dslaxpby(backend_data, this);
		else error("CommandQueue(%d)::run_operation(): unsupported opname = %s\n", id, opname)

	}
#ifdef UDEBUG
	fprintf(stderr, "[dev_id=%3d] <-----| CommandQueue(%d)::run_operation()\n", dev_id, id);
#endif
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

void CommandQueue::memcpy2DAsync(void* dest, long int ldest, void* src, long int ldsrc, long int rows, 
	long int cols, int elemSize, int loc_dest, int loc_src, int log_flag){
#ifdef DDEBUG
	fprintf(stderr, "CommandQueue(%d):memcpy2DAsync(dest=%p, ldest =%zu, src=%p, ldsrc = %zu, rows = %zu, cols = %zu, elemsize = %d, loc_dest = %d, loc_src = %d)\n",
		id, dest, ldest, src, ldsrc, rows, cols, elemSize, loc_dest, loc_src);
#endif	
#ifdef TTEST
	if(log_flag){
		while(__sync_lock_test_and_set(&my_log_lock, 1));
		if (n_trans_ctr > 100000) 
			error("CommandQueue::memcpy2DAsync(dest = %d, src = %d) exeeded 100000 transfers in TTEST Mode\n", 
				loc_dest, loc_src);
		if(!n_trans_ctr) n_reseTTEST();
		n_bytes[n_trans_ctr] = rows*cols*elemSize;
		n_locs[n_trans_ctr][0] = loc_src;
		n_locs[n_trans_ctr][1] = loc_dest;
	}
#endif
#ifndef PRODUCTION 
	massert(loc_dest >= -1 && loc_dest < CHL_MEMLOCS, 
		"CommandQueue(%d)::memcpy2DAsync: Invalid destination device: %d\n", id, loc_dest);
	massert(loc_src >= -1 && loc_src < CHL_MEMLOCS, 
		"CommandQueue(%d)::memcpy2DAsync: Invalid source device: %d\n", id, loc_src);
	if (loc_src == loc_dest) warning("CommandQueue(%d)::memcpy2DAsync(dest=%p, ldest =%zu, src=%p,"
	"ldsrc = %zu, rows=%zu, cols=%zu, elemSize =%d, loc_dest=%d, loc_src=%d): Source location matches destination\n",
	id, dest, ldest, src, ldsrc, rows, 
	cols, elemSize, loc_dest, loc_src);
#endif
#ifdef TTEST
	if(log_flag){
		CHLSetTimerAsync(&(n_timers[n_trans_ctr][0]));
		add_host_func((void*)&CHLSetTimerAsync, 
			(void*) &(n_timers[n_trans_ctr][1]));
	}
#endif	
	cudaStream_t stream = *((cudaStream_t*) backend_queue_ptr);
#ifndef PRODUCTION 
	//if (loc_src == -1 && loc_dest >=0) massert(CUBLAS_STATUS_SUCCESS == cublasSetMatrixAsync(rows, cols, elemSize, src, ldsrc, dest, ldest, stream), "CHLMemcpy2DAsync: cublasSetMatrixAsync failed\n");
	//else if (loc_src >=0 && loc_dest == -1) massert(CUBLAS_STATUS_SUCCESS == cublasGetMatrixAsync(rows, cols, elemSize, src, ldsrc, dest, ldest, stream),  "CHLMemcpy2DAsync: cublasGetMatrixAsync failed");
	//else 
	massert(cudaSuccess == 
#endif 
	cudaMemcpy2DAsync(dest, ldest*elemSize, 
	src, ldsrc*elemSize,
		rows*elemSize, cols, cudaMemcpyDefault, stream)
#ifndef PRODUCTION
	, "CHLMemcpy2DAsync(dest=%p, ldest =%zu, src=%p, ldsrc = %zu\nrows = %zu, cols = %zu, elemsize = %d,"
		"loc_dest = %d, loc_src = %d): cudaMemcpy2DAsync failed\n",
			dest, ldest, src, ldsrc, rows, 
			cols, elemSize, loc_dest, loc_src);
#else 
	;
#endif
#ifdef TTEST
	if(log_flag){
		add_host_func((void*)&CHLSetTimerAsync, 
			(void*) &(n_timers[n_trans_ctr][2]));
		n_trans_ctr++;
		__sync_lock_release(&my_log_lock);
	}
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
	dev_id = dev_id_in;
	cudaError_t err = cudaEventCreate(( cudaEvent_t*) event_backend_ptr);
#ifndef PRODUCTION
	massert(cudaSuccess == err, "Event::Event() - %s\n", cudaGetErrorString(err));
#endif
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
	Event_num_device[dev_id]--;
	cudaError_t err = cudaEventDestroy(*(( cudaEvent_t*) event_backend_ptr));
#ifndef PRODUCTION
	massert(cudaSuccess == err, "Event(%d)::~Event() - %s\n", id, cudaGetErrorString(err));
#endif
	free(event_backend_ptr);
#ifdef UDDEBUG
	fprintf(stderr, "[dev_id=%3d] <-----| Event(%d)::~Event()\n", dev_id, id);
#endif
}

// Block the calling thread until event is COMPLETE (using cuda backend == no busy wait)
void Event::sync_barrier()
{
#ifdef UDDEBUG
	fprintf(stderr, "[dev_id=%3d] |-----> Event(%d)::sync_barrier()\n", dev_id, id);
#endif
	event_status Evstat = query_status();
	if (status == COMPLETE) return;
	else (status == RECORDED){
		cudaEvent_t cuda_event= *(cudaEvent_t*) event_backend_ptr;
		cudaError_t err = cudaEventSynchronize(cuda_event);
		status = COMPLETE;
#ifndef PRODUCTION
		massert(cudaSuccess == err, "Event::sync_barrier() - %s\n", cudaGetErrorString(err));
	}
	else if (status == UNRECORDED) { error("[dev_id=%3d] |-----> Event(%d)::sync_barrier()"
		"Tried to sync unrecorded event\n", dev_id, id);
#endif
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
		status = COMPLETE;
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
	if (status == COMPLETE || status == UNRECORDED) return status;  // Retain status
	cudaEvent_t cuda_event= *(cudaEvent_t*) event_backend_ptr;
	cudaError_t err = cudaEventQuery(cuda_event);
	if (err == cudaSuccess) status == COMPLETE; // Update status
	else if (err == cudaErrorNotReady); // Retain status
	else error("[dev_id=%3d] |-----> Event(%d)::query_status() - %s, status = %s\n", dev_id, id,
	cudaGetErrorString(err), print_event_status(status));
#ifdef UDDEBUG
	fprintf(stderr, "[dev_id=%3d] <-----| Event(%d)::query_status() = %s\n", dev_id, id, print_event_status(status));
#endif
	return status;
}

void Event::soft_reset(){
#ifdef UDDEBUG
	fprintf(stderr, "[dev_id=%3d] |-----> Event(%d)::soft_reset()\n", dev_id, id);
#endif
	// sync_barrier();
	error("Not updated for non-lazy events\n")
	status = UNRECORDED;
#ifdef UDDEBUG
	fprintf(stderr, "[dev_id=%3d] <-----| Event(%d)::soft_reset()\n", dev_id, id);
#endif
}

void Event::reset(){
#ifdef UDDEBUG
	fprintf(stderr, "[dev_id=%3d] Event(%d)::reset()\n", dev_id, id);
#endif
	sync_barrier();
	cudaError_t err = cudaEventDestroy(*(( cudaEvent_t*) event_backend_ptr));
#ifndef PRODUCTION
	massert(cudaSuccess == err, "[dev_id=%3d] Event(%d)::reset - %s\n", dev_id, id, cudaGetErrorString(err));
#endif
	err = cudaEventCreate(( cudaEvent_t*) event_backend_ptr);
#ifndef PRODUCTION
	massert(cudaSuccess == err, "Event::Event() - %s\n", cudaGetErrorString(err));
#endif
	status = UNRECORDED;
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
