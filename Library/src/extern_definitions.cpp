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

void ProblemMetadata::print(){
    fprintf(stderr,"ProblemMetadata::print():\n------------------------------------\nAutotuner->\n"); 
    autotuner->print();
    fprintf(stderr,"problem_name: %s\n", problem_name); 
    fprintf(stderr,"problem_wrap: %p\n", problem_wrap); 
    fprintf(stderr,"Decomposers (%d)\n", decom_num); 
    for(int idx = 0; idx < decom_num; idx++)
        fprintf(stderr,"Decom %d -> adrs = %p\n", idx, decom[idx]->adrs); 
    fprintf(stderr,"sk_num: %d\n", sk_num);
    fprintf(stderr,"subkernel_list: %p\n", subkernel_list); 
    fprintf(stderr,"sk_dev_num: %s\n", printlist(sk_dev_num, CHL_WORKERS)); 
    fprintf(stderr,"subkernel_dev_list: %p\n", subkernel_dev_list);
    fprintf(stderr,"SAB: %p\n", SAB);
    fprintf(stderr,"\n------------------------------------\n"); 
}