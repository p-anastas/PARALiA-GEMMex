///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief The headers for functions for general use throught CHL
///

#ifndef SMART_WRAPPERS_H
#define SMART_WRAPPERS_H

#include<iostream>
//#include <cstdlib>

#include <string>
#include <cstring>

#include <stdio.h>
#include <stdarg.h>

#include <atomic>
//#define _GNU_SOURCE 
#include <pthread.h>

extern int CHL_MEMLOCS;

#include <execinfo.h>
#include <signal.h>
#include <unistd.h>
void handler(int sig);

/*****************************************************/
/// Generalised helpers

// Return current device used
int CHLGetDevice();

// Select device for current running pthread
void CHLSelectDevice(short dev_id);

// Return free memory and total memory for current device
void CHLDevGetMemInfo(long long* free_dev_mem, long long* max_dev_mem);

// Sync all devices and search for enchountered errors.
void CHLSyncCheckErr();

// Search for enchountered errors without synchronization.
void CHLASyncCheckErr();

// Enable available links for target device with all other devices
void CHLEnableLinks(short target_dev_i, short num_devices);

void CHLSetTimerAsync(void* wrapped_timer_Ptr);

/*****************************************************/
/// Print functions
#if !defined(PRINTFLIKE)
#if defined(__GNUC__)
#define PRINTFLIKE(n,m) __attribute__((format(printf,n,m)))
#else
#define PRINTFLIKE(n,m) /* If only */
#endif /* __GNUC__ */
#endif /* PRINTFLIKE */

template<typename VALUETYPE>
extern const char *printlist(VALUETYPE *list, int length);

void lprintf(short lvl, const char *fmt, ...)PRINTFLIKE(2,3);
void massert(bool condi, const char *fmt, ...)PRINTFLIKE(2,3);
void error(const char *fmt, ...)PRINTFLIKE(1,2);
void warning(const char *fmt, ...)PRINTFLIKE(1,2);

/*****************************************************/
/// Enum(and internal symbol) functions
// Memory layout struct for matrices
enum mem_layout { ROW_MAJOR = 0, COL_MAJOR = 1 };
const char *print_mem(mem_layout mem);

/*****************************************************/
/// General benchmark functions

inline double Drandom() { return (double)rand() / (double)RAND_MAX;}
short Dtest_equality(double* C_comp, double* C, long long size);
short Stest_equality(float* C_comp, float* C, long long size);

double Gval_per_s(long long value, double time);

long long gemm_ops(long int M, long int N, long int K);
long long gemm_mem_ops(long int M, long int N, long int K);
long long gemm_mem_sz(long int M, long int N, long int K, int elemSize);

long int count_lines(FILE* fp); // TODO: Where is this used?
int check_benchmark(char *filename);

/*****************************************************/
int gcd (int a, int b, int c);
long factorial(const int n);

void DECOM_2D(int val, int* valgrid_1_ptr, int* valgrid_2_ptr);

int is_in_list(int elem, int* elem_list, int list_len);
		
void translate_binary_to_unit_list(int case_id, int* active_unit_num_p, int* active_unit_id_list);
int translate_unit_list_to_binary(int* active_unit_id_list, int active_unit_num);
int is_subset(int case_id, int case_id_set);
int binary_case_id_split(int case_id);

/// double/float arethmetic comparrison that compares based on NORMALIZE_NEAR_SPLIT_LIMIT minimum difference.
int normal_equal(double n1, double n2);
int normal_less(double n1, double n2);
int normal_larger(double n1, double n2);
int normal_lessequal(double n1, double n2);
int normal_largerequal(double n1, double n2);

void bubble_sort_3_vals(int* case_ids, double vals[][3], double weight[3], int len);
void bubble_sort_2_idx_3x3_vals(int* case_ids, int* case_ids_2, double vals[][3], double vals_2[][3], double vals_3[][3], double weight[3], int len);

/*****************************************************/
/// Generalised "Command queue" and "Event" definition (e.g. CUDA streams and Events for CUDA backend)

typedef class Event* Event_p;

enum CQueue_type{
	COMMUNICATION = 0,
	COMPUTATION = 1,
	MIXED = 2
};

typedef class CommandQueue
{
	private:
	public:
		// The (incremental) id of the queue
		int id;  
		// The Worker that this queue's backend is close to. 
		int dev_id;

		void* backend_queue_ptr, *backend_comp_md;

		/// Queues can support communication, computation or both. This is controlled by 'type'
		// type = 0: communication-only queue, 
		// type = 1: Execution-only queue backend (backend_exec_md defined)
		// type = 2: Both comm and comp queue.
		CQueue_type type; 

		//Constructor (task_id: 1)
		CommandQueue(int dev_id_in, CQueue_type type_in);

		//Destructor (task_id: -42)
		~CommandQueue();
		//Synchronize the queue backend (task_id: 2)
		void sync_barrier();
		//Add a host function to the queue - pthreadception (task_id: 3)
		void add_host_func(void* func, void* data);
		// Force the queue to wait for an event to conclude (task_id: 4)
		void wait_for_event(Event_p Wevent);

		// Record an event to the queue (task_id: 5)
		void record_event(Event_p Wevent);

		// Run a certain backend function at (mode 1 or 2) queue (task_id: 6)
		void run_operation(void* backend_data, const char* opname, int target_dev_id);

		void memcpy2DAsync(void* dest, long int ldest, void* src, long int ldsrc, 
		long int rows, long int cols, int elemSize, int loc_dest, int loc_src, int TTs_log_flag);
		std::string name;
		void print() { std::cout << "Command Queue : " << name; }

}* CQueue_p;

enum event_status{
	UNRECORDED = 0, /// Event has not been recorded yet.
	RECORDED = 1, /// Recorded but not guaranteed to be complete.
	COMPLETE = 2,  /// Complete but not yet ran 'check' function (for cache updates etc)
};

/// Returns a string representation for event_status
const char* print_event_status(event_status in_status);

//extern long int Event_num_loc[64];
//extern Event_p event_pools[64];

class Event
{
	private:
	public:
		event_status status;
		
		void* event_backend_ptr;
		int id, dev_id;

		/// Constructors
		Event();
		/// Destructors
		~Event();
		/// Functions
		void sync_barrier();
		void record_to_queue(CQueue_p Rr);
		event_status query_status();
		void checked();
		void reset();
		void soft_reset();

};

/*****************************************************/
/// Generalised data management functions

//extern long PAGE_sz;

char* mem_name(int idx);

// Checks if *ptr points to interleaved memory or not (anything else in host)
int get_hostmem_idx(void* ptr);

// Supposedly returns an estimation of the free and max host memory available in the system
void CHLGetMaxCPUmem(long long int* free_mem, long long int* max_mem);

// Pin *ptr
void pin_mem_wrap(void** ptr, long long bytes);

// Malloc in pinned host memory
void* pin_malloc(long long N_bytes, short W_flag);

// Free alloc in pinned host memory
void pin_free(void *gpuptr);

// Malloc in numa-interleaved memory and pin
void* numa_inter_pin_malloc(long long N_bytes, short W_flag);

// Malloc in 'node_num' numa node memory and pin
void* numa_bind_pin_malloc(long long count, int node_num, short W_flag);

// Unpin and free alloc in numa-interleaved memory
void numa_pin_free(void *gpuptr, long long bytes);

// Malloc in loc with error-checking
void* CHLMalloc(long long N_bytes, int loc, short W_flag);

// Malloc in Host and perform openmp default first touch. 
void* CHLMallocHostTouchDefault(long long bytes);

void* CHLMallocHostTouchSerial(long long bytes);

void* CHLMallocHostTouchSmart(long long dim1, long long dim2, int elemSize, char transpose);

// Free in loc with error-checking
void CHLFree(void * ptr, long long bytes, int loc);

short CHLGetPtrLoc(void * in_ptr);
//short CHLGetPtrAdvLoc(void * in_ptr, long long dim1, long long dim2, int elemSize);

// Memcpy between two locations with errorchecking
void CHLMemcpy(void* dest, void* src, long long N_bytes, int loc_dest, int loc_src);

// Memcpy between two locations with errorchecking
void CHLMemcpy2D(void* dest, long int ldest, void* src, long int lsrc, long int rows, long int cols, short elemSize, int loc_dest, int loc_src);

/* Deprecated...replaced by COMMUNICATION queues
// Asunchronous Memcpy between two locations WITHOUT synchronous errorchecking. Use with caution.
void CHLMemcpyAsync(void* dest, void* src, long long N_bytes, int loc_dest, int loc_src, CQueue_p transfer_medium);

// Asunchronous Memcpy between two locations WITHOUT synchronous errorchecking.
void CHLMemcpy2DAsync(void* dest, long int ldest, void* src, long int lsrc, long int rows, long int cols,
	short elemSize, int loc_dest, int loc_src, CQueue_p transfer_medium);
// Asunchronous Memcpy between two locations WITHOUT synchronous errorchecking and with TTEST logging disabled (for avoiding double logs for hop tranfers)
void CHLMemcpy2DAsync_noTTs(void* dest, long int ldest, void* src, long int lsrc, long int rows, long int cols,
	short elemSize, int loc_dest, int loc_src, CQueue_p transfer_medium);
*/

// Print and log bandwidths and links used with CHLMemcpy2DAsync. Unusable with TTEST flag
void n_HopMemcpyPrint();

// Initalize vector in loc with error-checking
template<typename VALUETYPE>
extern void CHLVecInit(VALUETYPE *vec, long long length, int seed, int loc);
// Helper for Parallel OpenMP vector initialization
template<typename VALUETYPE>
extern void CHLParallelVecInitHost(VALUETYPE *vec, long long length, int seed);

// Return the max dim size (which is a multiple of 'step') for 'Asset2DNum' square assets on 'loc'
long int CHLGetMaxDimSqAsset2D(short Asset2DNum, short dsize, long int step, int loc);

// Return the max dim size (which is a multiple of 'step') for 'Asset1DNum' assets on 'loc'
long int CHLGetMaxDimAsset1D(short Asset1DNum, short dsize, long int step, int loc);

/*****************************************************/
/// Timers for benchmarks
// CPU accurate timer
double csecond();

// Event timer for background Event timing (~usually ms accuracy)
typedef class Event_timer
{
	private:
		Event_p Event_start;
		Event_p Event_stop;
		double time_ms;
	public:
		int dev_id;

		Event_timer(int dev_id);
		void start_point(CQueue_p start_queue);
		void stop_point(CQueue_p stop_queue);
		double sync_get_time();
}* Event_timer_p;


/*****************************************************/

#endif