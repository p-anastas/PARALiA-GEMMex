#include "PARALiA.hpp"
#include "Subkernel.hpp"

//typedef class ProblemMetadata* PMD_p;
PMD_p PMD_cache[PROBLEM_MD_CACHE] = {NULL}; 
int PMD_cache_entries = 0;

CQueue_p recv_queues[64][64];
CQueue_p wb_queues[64][64];

CQueue_p exec_queue[32][MAX_BACKEND_L];
int exec_queue_ctr[32]; 

CQueue_p reduce_queue[64][REDUCE_WORKERS_PERDEV];
int reduce_queue_ctr[64]; 
int reduce_loc; 