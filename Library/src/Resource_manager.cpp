///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief 
///

#include "PARALiA.hpp"
#include "Resource_manager.hpp"
#include "Autotuner.hpp"
#include "DataCaching.hpp"
#include "backend_wrappers.hpp"

#include <cfloat>
#include <cmath>

int init = RMConfigResources();

int RMConfigResources(){
	for(int dev_id_idx = 0 ; dev_id_idx < 64; dev_id_idx++){
		for(int dev_id_idy = 0 ; dev_id_idy < 64; dev_id_idy++) 
			recv_queues[dev_id_idx][dev_id_idy] = wb_queues[dev_id_idx][dev_id_idy] = NULL; 
		for(int idx = 0 ; idx < REDUCE_WORKERS_PERDEV; idx++) reduce_queue[dev_id_idx][idx] = NULL;
	}
	for(int dev_id_idx = 0 ; dev_id_idx < 32; dev_id_idx++){
		for (int i = 0; i < MAX_BACKEND_L; i++) exec_queue[dev_id_idx][i] = NULL;
		exec_queue_ctr[dev_id_idx] = -1; 
	}
	return 1;
}

void RMInitResources(int* dev_list, int dev_num){
#ifdef DEBUG
	fprintf(stderr, "|-----> RMInitResources(%s)\n", printlist(dev_list, dev_num));
#endif
	for(int dev_id_idx = 0 ; dev_id_idx < CHL_WORKERS + 1; dev_id_idx++){
		for(int dev_id_idy = 0 ; dev_id_idy < CHL_WORKERS + 1; dev_id_idy++)
		if(dev_id_idy!=dev_id_idx){
			if (!recv_queues[dev_id_idx][dev_id_idy]){
				int cont_flag = 0; 
				if(best_grid_edge_active[dev_id_idx][dev_id_idy] != -1){
					cont_flag = 1;
					int queue_id = (dev_id_idy >= CHL_WORKERS)? (dev_id_idx) : (dev_id_idy);
					recv_queues[dev_id_idx][dev_id_idy] = new CommandQueue(queue_id, COMMUNICATION);	
#ifdef ENABLE_SEND_RECV_OVERLAP
					wb_queues[dev_id_idx][dev_id_idy] = new CommandQueue(queue_id, COMMUNICATION);
#else 
					wb_queues[dev_id_idx][dev_id_idy] = recv_queues[dev_id_idx][dev_id_idy];
#endif
				}
				else if(best_grid_edge_replaced[dev_id_idx][dev_id_idy][0]!= -1){
					cont_flag = 1;
					recv_queues[dev_id_idx][dev_id_idy] = recv_queues
						[best_grid_edge_replaced[dev_id_idx][dev_id_idy][0]]
						[best_grid_edge_replaced[dev_id_idx][dev_id_idy][1]];
					wb_queues[dev_id_idx][dev_id_idy] = wb_queues
						[best_grid_edge_replaced[dev_id_idx][dev_id_idy][0]]
						[best_grid_edge_replaced[dev_id_idx][dev_id_idy][1]];
				}
				if(cont_flag){
					if(dev_id_idy == CHL_WORKERS){ // The smallest index shared link allocates the queue
						for(int idx1 = CHL_WORKERS + 1 ; idx1 < CHL_MEMLOCS; idx1++){
							recv_queues[dev_id_idx][idx1] = recv_queues[dev_id_idx][dev_id_idy];
							wb_queues[dev_id_idx][idx1] = wb_queues[dev_id_idx][dev_id_idy];
						}
					}
					if(dev_id_idx == CHL_WORKERS){ // The smallest index shared link allocates the queue
						for(int idx1 = CHL_WORKERS + 1 ; idx1 < CHL_MEMLOCS; idx1++){
							recv_queues[idx1][dev_id_idy] = recv_queues[dev_id_idx][dev_id_idy];
							wb_queues[idx1][dev_id_idy] = wb_queues[dev_id_idx][dev_id_idy];
						}
					}
				}
			}
		}
		if (!exec_queue[dev_id_idx][0]) {
			int flag_is_worker = 0; 
			for (int i = 0; i < dev_num; i++) if(dev_list[i] == dev_id_idx){
				flag_is_worker = 1; 
				break;
			}
			if(flag_is_worker){
				for (int i = 0; i < MAX_BACKEND_L; i++){
					exec_queue[dev_id_idx][i] = new CommandQueue(dev_id_idx, COMPUTATION);
					exec_queue_ctr[dev_id_idx] = -1; 
				}
			}
		}
	}
	for(int dev_id_idx = 0 ; dev_id_idx < CHL_MEMLOCS; dev_id_idx++)	
		for (int i = 0; i < REDUCE_WORKERS_PERDEV; i++){
		int queue_id = (dev_id_idx >= CHL_WORKERS)? (CHL_MEMLOC_CLOSE_TO_WORKER[dev_id_idx]) : (dev_id_idx);
			if (!reduce_queue[dev_id_idx][i]) 
			//if(best_grid_edge_replaced[i][dev_id_idx][0] == -1) 
			reduce_queue[dev_id_idx][i] = new CommandQueue(queue_id, COMPUTATION);
			//else reduce_queue[dev_id_idx][i] = reduce_queue
			//			[best_grid_edge_replaced[dev_id_idx][i][0]]
			//			[best_grid_edge_replaced[dev_id_idx][i][1]];
			reduce_queue_ctr[dev_id_idx] = -1; 
		}
#ifdef DEBUG
	fprintf(stderr, "<-----|\n");
#endif
}

void RMInitWS(int* dev_list, int dev_num){
	error("CHLInitWS: INCOMPLETE\n");
	for(int dev_id_idx = 0 ; dev_id_idx < CHL_WORKERS; dev_id_idx++){
		if (!exec_queue[dev_id_idx][0]) {
			int flag_is_worker = 0; 
			for (int i = 0; i < dev_num; i++) if(dev_list[i] == (dev_id_idx)){
				flag_is_worker = 1; 
				break;
			}
			if(flag_is_worker && (dev_id_idx)!= -1){	
				void* local_ws = CHLMalloc(2048, (dev_id_idx), 1); 
				massert(CUBLAS_STATUS_SUCCESS == cublasSetWorkspace(*((cublasHandle_t*) 
					exec_queue[dev_id_idx][0]->backend_comp_md), local_ws, 2048), 
					"CommandQueue::CommandQueue(%d): cublasSetWorkspace failed\n", (dev_id_idx));	
			}
		}
	}
}

void RMFreeResources(){
	for(int dev_id_idx = 0 ; dev_id_idx < CHL_MEMLOCS; dev_id_idx++){
		for(int dev_id_idy = 0 ; dev_id_idy < CHL_MEMLOCS; dev_id_idy++)
		if(dev_id_idx!=dev_id_idy){
				if(recv_queues[dev_id_idx][dev_id_idy]) delete recv_queues[dev_id_idx][dev_id_idy];
				recv_queues[dev_id_idx][dev_id_idy] = NULL;
				if(wb_queues[dev_id_idx][dev_id_idy]) delete wb_queues[dev_id_idx][dev_id_idy];
				wb_queues[dev_id_idx][dev_id_idy] = NULL;
		}
		for (int i = 0; i < MAX_BACKEND_L; i++){
			if(exec_queue[dev_id_idx] && exec_queue[dev_id_idx][i]) delete exec_queue[dev_id_idx][i];
			exec_queue[dev_id_idx][i] = NULL;
		}
	}
}

void RMCleanResources(){
	for(int dev_id_idx = 0 ; dev_id_idx < CHL_MEMLOCS; dev_id_idx++){
		for(int dev_id_idy = 0 ; dev_id_idy < CHL_MEMLOCS; dev_id_idy++)
		for (int i = 0; i < MAX_BACKEND_L; i++)
			if(exec_queue[dev_id_idx]  && exec_queue[dev_id_idx][i]){
				exec_queue_ctr[dev_id_idx] = -1;
			}
	}
}

void RMSyncRecvQueues(){
	for(int dev_id_idx = 0 ; dev_id_idx < CHL_MEMLOCS; dev_id_idx++)
		for(int dev_id_idy = 0 ; dev_id_idy < CHL_MEMLOCS; dev_id_idy++)
		if(dev_id_idx!=dev_id_idy) recv_queues[dev_id_idx][dev_id_idy]->sync_barrier();
}

ProblemMetadata::~ProblemMetadata(){
	delete autotuner;
	//free((void*)problem_name);
	free(problem_wrap);
	for (int idx = 0; idx < decom_num; idx++) delete decom[idx];
	for (int idx = 0; idx < 64; idx++) if (SAB[idx]){
		delete SAB[idx];
		SAB[idx] = NULL;
	}
}

void PARALiADevCacheFree(int dev_id){
#ifdef DEBUG
	fprintf(stderr, "|-----> PARALiADevCacheFree(%d)\n", dev_id);
#endif
	for(int i = 0; i < PMD_cache_entries; i++) 
		if (PMD_cache[i]->SAB[dev_id] == current_SAB[dev_id]){
			PMD_cache[i]->SAB[dev_id] = NULL;
	}
	if(current_SAB[dev_id]) delete current_SAB[dev_id];
	current_SAB[dev_id] = NULL;
	/*for(int i = 0; i < PMD_cache_entries; i++) 
		if (PMD_cache[i]->SAB[dev_id]){
 			delete PMD_cache[i]->SAB[dev_id];
			PMD_cache[i]->SAB[dev_id] = NULL;
	}*/
#ifdef DEBUG
	fprintf(stderr, "<-----|\n");
#endif
}