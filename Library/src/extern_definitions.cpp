#include "PARALiA.hpp"
#include "Resource_manager.hpp"

//typedef class ProblemMetadata* PMD_p;
PMD_p PMD_cache[PROBLEM_MD_CACHE] = {NULL}; 
int PMD_cache_entries = 0;
int conserve_memory_curr = 0; 

CQueue_p recv_queues[64][64];
CQueue_p wb_queues[64][64];

CQueue_p exec_queue[32][MAX_BACKEND_L_IN];
int exec_queue_ctr[32]; 

CQueue_p reduce_queue[64][REDUCE_WORKERS_PERDEV];
int reduce_queue_ctr[64]; 
int reduce_loc; 

void ProblemMetadata::print(){
    fprintf(stderr,"ProblemMetadata::print():\n------------------------------------\nAutotuner->\n"); 
    for(int itter = 0; itter < REP_TILE; itter++){
        fprintf(stderr,"%d: ", itter); 
        autotuner[itter]->print();
    }
    fprintf(stderr,"problem_name: %s\n", problem_name); 
    fprintf(stderr,"problem_wrap: %p\n", problem_wrap); 
    fprintf(stderr,"Decomposers (%d)\n", decom_num); 
    for(int idx = 0; idx < decom_num; idx++)
        fprintf(stderr,"Decom %d -> adrs = %p\n", idx, decom[idx]->adrs); 
    fprintf(stderr,"SAB: %p\n", SAB);
    fprintf(stderr,"\n------------------------------------\n"); 
}