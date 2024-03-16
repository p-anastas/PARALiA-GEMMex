
#include <iostream>
#include <cfloat>
#include <climits>
#include <algorithm>
#include <iostream>
#include <list>

#include "chl_smart_wrappers.hpp"
#include "Autotuner.hpp"

double get_edge_bw(int dest_loc, int src_loc){
	if(best_grid_edge_active[dest_loc][src_loc]!= -1){
		if(best_grid_edge_bws[dest_loc][src_loc] != -1) return best_grid_edge_bws[dest_loc][src_loc];
		else error("get_edge_bw(dest=%d,src=%d) best_grid_edge_active[%d][%d] = %d"
			"but best_grid_edge_bws[%d][%d] = %lf\n", dest_loc, src_loc, dest_loc, src_loc, 
			best_grid_edge_active[dest_loc][src_loc], dest_loc, src_loc, best_grid_edge_bws[dest_loc][src_loc]);
  	}
  	else if(best_grid_edge_active[dest_loc][src_loc] ==-1 && best_grid_edge_replaced[dest_loc][src_loc][0] == -1){
		warning("get_edge_bw(dest=%d,src=%d) called but edge is both inactive AND not replaced\n", dest_loc, src_loc);
		return 1e-9;
  	}
  	else return best_grid_edge_bws[best_grid_edge_replaced[dest_loc][src_loc][0]]
								  [best_grid_edge_replaced[dest_loc][src_loc][1]];
	return -1.0; // Should never reach, just to remove warning
}

long double LinkRoute::optimize(int* loc_map, long int size, int update_flag){
	if(!strcmp(FETCH_ROUTING, "P2P_FETCH_FROM_INIT")) return optimize_p2p_init(loc_map, size);
	else if(!strcmp(FETCH_ROUTING, "P2P_FETCH_FROM_GPU_SERIAL")) return optimize_p2p_serial(loc_map, size);
	else if(!strcmp(FETCH_ROUTING, "P2P_FETCH_FROM_GPU_DISTANCE")) return optimize_p2p_distance(loc_map, size);
	else if(!strcmp(FETCH_ROUTING, "CHAIN_FETCH_SERIAL")) return optimize_chain_serial(loc_map, size);
	else if(!strcmp(FETCH_ROUTING, "CHAIN_FETCH_RANDOM")) return optimize_chain_random(loc_map, size);
	else if(!strcmp(FETCH_ROUTING, "CHAIN_FETCH_TIME")) return optimize_chain_time(loc_map, size);
	else if(!strcmp(FETCH_ROUTING, "CHAIN_FETCH_QUEUE_WORKLOAD")) return optimize_chain_ETA(loc_map, size, update_flag);
	else error("LinkRoute::optimize() -> %s not implemented", FETCH_ROUTING);
	return 0;
}

// Naive fetch from initial data loc every time. Similar to cuBLASXt
long double LinkRoute::optimize_p2p_init(int* loc_map, long int size){
#ifdef DEBUG
	fprintf(stderr, "|-----> LinkRoute::optimize_p2p_init()\n");
#endif
	hop_num = 2;
	int start_hop = -42, end_hop = -42;
	for(int ctr = 0; ctr < CHL_MEMLOCS; ctr++)
	{
		if(loc_map[ctr] == 0) start_hop = ctr;
		if(loc_map[ctr] == 2) end_hop = ctr;
	}
	hop_uid_list[0] = start_hop;
	hop_uid_list[1] = end_hop;
	starting_hop = 0; 
	loc_map[end_hop] = 42;
#ifdef DEBUG
	fprintf(stderr, "<-----|\n");
#endif
	return 0;
}

// Outdated fetch with preference to GPU tiles but with simple serial search for src
// Similar to BLASX behaviour when its assumed topology does not fit to the interconnect
long double LinkRoute::optimize_p2p_serial(int* loc_map, long int size){
#ifdef DEBUG
	fprintf(stderr, "|-----> LinkRoute::optimize_p2p_serial()\n");
#endif
	hop_num = 2;
	int start_hop = -42, end_hop = -42;
	for(int ctr = 0; ctr < CHL_MEMLOCS; ctr++)
	{
	  if(loc_map[ctr] == 0 || loc_map[ctr] == 42) start_hop = ctr;
	  if(loc_map[ctr] == 2) end_hop = ctr;
	  if(start_hop!= -42 && end_hop != -42) break;
	}
	hop_uid_list[0] = start_hop;
	hop_uid_list[1] = end_hop;
	starting_hop = 0; 
	loc_map[end_hop] = 42;
#ifdef DEBUG
	fprintf(stderr, "<-----|\n");
#endif
	return 0;
}

// Fetch selection based on 'distance' from available sources (if multiple exist)
// Similar to XKBLAS and PARALiA 1.5
long double LinkRoute::optimize_p2p_distance(int* loc_map, long int size){
#ifdef DEBUG
	fprintf(stderr, "|-----> LinkRoute::optimize_p2p_distance()\n");
#endif
	hop_num = 2;
	int start_hop = -42, end_hop = -42;
	for(int ctr = 0; ctr < CHL_MEMLOCS; ctr++) if(loc_map[ctr] == 2) end_hop = ctr;
	int pos_max = CHL_MEMLOCS;
	double link_bw_max = 0;
	for (int pos =0; pos < CHL_MEMLOCS; pos++) if (loc_map[pos] == 0 || loc_map[pos] == 42){
		double current_link_bw = get_edge_bw(end_hop, pos);
		//fprintf(stderr, "%lf\n", current_link_bw);
		if (current_link_bw > link_bw_max){
		  link_bw_max = current_link_bw;
		  pos_max = pos;
		}
	}
	hop_uid_list[0] = pos_max;
	hop_uid_list[1] = end_hop;
	starting_hop = 0; 
	loc_map[end_hop] = 42;
#ifdef DEBUG
	fprintf(stderr, "<-----|\n");
#endif
	return 0;
}

// Naive, for comparison reasons mainly
long double LinkRoute::optimize_chain_serial(int* loc_map, long int size){
#ifdef DEBUG
	fprintf(stderr, "|-----> LinkRoute::optimize_chain_serial()\n");
#endif
	hop_num = 1;
	for(int ctr = 0; ctr < CHL_MEMLOCS; ctr++) 
		if(loc_map[ctr] == 0) hop_uid_list[0] = ctr;
		else if(loc_map[ctr] == 1 || loc_map[ctr] == 2){
			hop_uid_list[hop_num++] = ctr;
			loc_map[ctr] = 42;
		}
#ifdef DEBUG
	fprintf(stderr, "<-----|\n");
#endif
	return 0;
}

long double LinkRoute::optimize_chain_random(int* loc_map, long int size){
#ifdef DEBUG
	fprintf(stderr, "|-----> LinkRoute::optimize_chain_random()\n");
#endif
	hop_num = 0;
	int loc_list[CHL_MEMLOCS];
	for(int ctr = 0; ctr < CHL_MEMLOCS; ctr++){
	  if(loc_map[ctr] == 0) hop_uid_list[0] = ctr;
	  else if(loc_map[ctr] == 1 || loc_map[ctr] == 2)
		loc_list[hop_num++] = ctr;
	} 
	int start_idx = int(rand() % hop_num); 
	int hop_ctr = 1;
	for(int ctr = start_idx; ctr < hop_num; ctr++) hop_uid_list[hop_ctr++] = loc_list[ctr];
	for(int ctr = 0; ctr < start_idx; ctr++) hop_uid_list[hop_ctr++] = loc_list[ctr];
	hop_num++;
	for(int ctr = 1; ctr < hop_num; ctr++) loc_map[hop_uid_list[ctr]] = 42;
#ifdef DEBUG
	fprintf(stderr, "<-----|\n");
#endif
	return 0;
}

long double LinkRoute::optimize_chain_time(int* loc_map, long int size){
#ifdef DEBUG
	fprintf(stderr, "|-----> LinkRoute::optimize_chain_time()\n");
#endif
	hop_num = 0;
	std::list<int> loc_list;
	int tmp_hop = -42; 
	double fire_est = csecond();
	for(int ctr = 0; ctr < CHL_MEMLOCS; ctr++){
	  if(loc_map[ctr] == 0) hop_uid_list[0] = ctr;
	  else if(loc_map[ctr] == 1 || loc_map[ctr] == 2){
		tmp_hop = ctr;
		loc_list.push_back(tmp_hop);
		hop_num++;
	  }
	}
	double best_t = 1e9; 
	if (hop_num == 1){
	  hop_uid_list[1] = tmp_hop;
	  best_t = size/(1e9*get_edge_bw(hop_uid_list[1], hop_uid_list[0]));
	}
	else{
	  int best_list[factorial(hop_num)][hop_num]; 
	  int flag = 1, tie_list_num = 0;
	  while (flag){
		double max_t = -1, total_t = 0, temp_t;
		int temp_ctr = 0, prev = (hop_uid_list[0]), templist[hop_num]; 
		for (int x : loc_list){
		  templist[temp_ctr++] = x; 
		  temp_t = size/(1e9*get_edge_bw(x, prev));
		  if (temp_t > max_t) max_t = temp_t;
		  total_t += temp_t;
		  prev = (x);
		}
		temp_t = max_t;// + (total_t - max_t)/STREAMING_BUFFER_OVERLAP;
		//fprintf(stderr,"Checking location list[%s]: temp_t = %lf\n", printlist(templist,hop_num), temp_t);
	  
		//fprintf(stderr,"Checking location list[%s]: temp_ETA = %lf\n", printlist(templist,hop_num), temp_ETA);
		if(temp_t < best_t){
		  best_t = temp_t;
		  for (int ctr = 0; ctr < hop_num; ctr++) best_list[0][ctr] = templist[ctr];
		  tie_list_num = 1;
		}
		else if (temp_t == best_t){
		  for (int ctr = 0; ctr < hop_num; ctr++) best_list[tie_list_num][ctr] = templist[ctr];
		  tie_list_num++;
		}
		flag = std::next_permutation(loc_list.begin(), loc_list.end());
	  }
	  //fprintf(stderr,"Selecting location list[%s]: best_t = %lf\n", printlist(best_list,hop_num), best_t);
	  int rand_tie_list = int(rand() % tie_list_num); 
	  for(int ctr = 0; ctr < hop_num; ctr++)
		hop_uid_list[ctr+1] = best_list[rand_tie_list][ctr];
	}
	hop_num++; 
	for(int ctr = 1; ctr < hop_num; ctr++) loc_map[hop_uid_list[ctr]] = 42;
	return best_t;
}

/// PARALia 3.0 - simple timed queues without slowdowns
// An estimation of when the queue will be free of tasks.
typedef class P2P_queue_load{
	long double queue_ETA[64][64];
public:
	P2P_queue_load();
	void ETA_add_task(int dest, int src, long double task_fire_t, long double task_duration);
	void ETA_set(int dest, int src, long double new_ETA);
	long double ETA_get(int dest, int src);
}* Queue_load_p;
Queue_load_p queue_load_grid = NULL;


/*****************************************************/
/// PARALia 2.0 - timed queues

P2P_queue_load::P2P_queue_load(){
	for(int idx = 0; idx < 64; idx++) for(int idy = 0; idy < 64; idy++) queue_ETA[idx][idy] = 0.0;
}

void P2P_queue_load::ETA_add_task(int dest, int src, long double task_fire_t, long double task_duration){
	queue_ETA[dest][src] = std::max(queue_ETA[dest][src], task_fire_t) + task_duration;
}

void P2P_queue_load::ETA_set(int dest, int src,long double new_ETA){
	queue_ETA[dest][src] = new_ETA; 
}

long double P2P_queue_load::ETA_get(int dest, int src){
	return queue_ETA[dest][src];
}

/*****************************************************/

long double LinkRoute::optimize_chain_ETA(int* loc_map, long int size, int update_flag){
	if(!queue_load_grid) queue_load_grid = new P2P_queue_load();
	long double min_ETA = DBL_MAX, tile_t = DBL_MAX/100, fire_t = 0;//csecond();
	hop_num = 0;
	std::list<int> loc_list;
	int tmp_hop = -42; 
	for(int ctr = 0; ctr < CHL_MEMLOCS; ctr++){
		if(loc_map[ctr] == 0) hop_uid_list[0] = ctr;
		else if(loc_map[ctr] == 1 || loc_map[ctr] == 2){
			tmp_hop = ctr;
			loc_list.push_back(tmp_hop);
			hop_num++;
		}
	}
	if (hop_num == 1){
		hop_uid_list[1] = tmp_hop;
//#ifdef ENABLE_TRANSFER_HOPS
//	  min_ETA = optimize_hop_route(transfer_tile_wrapped, update_ETA_flag, hop_uid_list[1], hop_uid_list[0]);
//#else
		long double temp_t = size/(1e9*get_edge_bw(hop_uid_list[1], hop_uid_list[0]));
		//fprintf(stderr, "%lf\n", temp_t);
		if (update_flag) queue_load_grid->ETA_add_task(hop_uid_list[1], hop_uid_list[0], fire_t, temp_t);
		min_ETA = std::max(queue_load_grid->ETA_get(hop_uid_list[1], hop_uid_list[0]), fire_t) + temp_t;
		hop_num++;
//#endif
	}
	else{
		int best_list[factorial(hop_num)][hop_num]; 
		int flag = 1, tie_list_num = 0;
		while (flag){
			long double max_t = -1, total_t = 0, temp_t;
			int temp_ctr = 0, prev = hop_uid_list[0], templist[hop_num]; 
			for (int x : loc_list){
				templist[temp_ctr++] = x; 
				temp_t = 1.0 * size/(1e9*get_edge_bw(x, prev));
				if (temp_t > max_t) max_t = temp_t;
				total_t += temp_t;
				prev = x;
			}
			temp_t = max_t + (total_t - max_t)/STREAMING_BUFFER_OVERLAP;
			//fprintf(stderr, "%llf\n", temp_t);
			//fprintf(stderr,"Checking location list[%s]: temp_t = %lf\n", printlist(templist,hop_num), temp_t);
			prev = hop_uid_list[0];
			long double temp_ETA = 0, queue_ETA; 
			for(int ctr = 0; ctr < hop_num; ctr++){
				queue_ETA = std::max(queue_load_grid->ETA_get(templist[ctr], prev), fire_t) + temp_t;
				//fprintf(stderr,"queue_ETA [%d -> %d]= %lf (temp_t = %lf)\n", prev), templist[ctr], queue_ETA);
				prev = templist[ctr];
				if(temp_ETA < queue_ETA) temp_ETA = queue_ETA;
				//fprintf(stderr,"Checking location list[%s]: queue_ETA = %lf\n", printlist(templist,hop_num), queue_ETA);
			}
			//fprintf(stderr,"Checking location list[%s]: temp_ETA = %lf\n", printlist(templist,hop_num), temp_ETA);
			if(temp_ETA < min_ETA){// && BANDWIDTH_DIFFERENCE_CUTTOF_RATIO*tile_t >= temp_t){
				//if(abs(temp_ETA - min_ETA)/temp_t > NORMALIZE_NEAR_SPLIT_LIMIT && temp_ETA < min_ETA){
				min_ETA = temp_ETA;
				tile_t = temp_t;
				for (int ctr = 0; ctr < hop_num; ctr++) best_list[0][ctr] = templist[ctr];
				tie_list_num = 1;
#ifdef DPDEBUG
				fprintf(stderr,"DataTile[%d:%d,%d]: New min_ETA(%llf) for route = %s\n", 
				transfer_tile->id, transfer_tile->GridId1, transfer_tile->GridId2, min_ETA, 
				printlist(best_list[tie_list_num-1],hop_num));
#endif
			}
			else if (temp_ETA == min_ETA){
			//else if(abs(temp_ETA - min_ETA)/temp_t <= NORMALIZE_NEAR_SPLIT_LIMIT){
			for (int ctr = 0; ctr < hop_num; ctr++) best_list[tie_list_num][ctr] = templist[ctr];
			tie_list_num++;
	#ifdef DPDEBUG
			fprintf(stderr,"DataTile[%d:%d,%d]: same min_ETA(%llf) for candidate(%d) route = %s\n", 
				transfer_tile->id, transfer_tile->GridId1, transfer_tile->GridId2, temp_ETA, 
				tie_list_num, printlist(best_list[tie_list_num-1],hop_num));
	#endif
			}
			flag = std::next_permutation(loc_list.begin(), loc_list.end());
		}
		
		int rand_tie_list = int(rand() % tie_list_num); 
		for(int ctr = 0; ctr < hop_num; ctr++){
				hop_uid_list[ctr+1] = best_list[rand_tie_list][ctr];
				if (update_flag) queue_load_grid->ETA_set(hop_uid_list[ctr+1], hop_uid_list[ctr], min_ETA);
		}
		hop_num++;
#ifdef PDEBUG
		fprintf(stderr,"Selected route = %s from %d candidates with ETA = %llf\n", 
				printlist(hop_uid_list, hop_num), tie_list_num, min_ETA);
#endif
	}
	if (update_flag) for(int ctr = 1; ctr < hop_num; ctr++) loc_map[hop_uid_list[ctr]] = 42;
	return min_ETA; 
}

long double LinkRoute::optimize_reverse(int* loc_map, long int size){
	if(!strcmp(WB_ROUTING, "P2P_TO_INIT")) return optimize_reverse_p2p_init(loc_map, size);
	else error("LinkRoute::optimize_reverse() -> %s not implemented", WB_ROUTING);
	return 0;
}

long double LinkRoute::optimize_reverse_p2p_init(int* loc_map, long int size){
	hop_num = 2;
	int start_hop = -42, end_hop = -42;
	for(int ctr = 0; ctr < CHL_MEMLOCS; ctr++)
	{
		if(loc_map[ctr] == 42) start_hop = ctr;
		if(loc_map[ctr] == 0) end_hop = ctr;
	}
	hop_uid_list[0] = start_hop;
	hop_uid_list[1] = end_hop;
	starting_hop = 0; 
	return 0;
}

/*
Subkernel* SubkernelSelect(int dev_id, Subkernel** Subkernel_list, long Subkernel_list_len){
	Subkernel* curr_sk = NULL;
	int task_idx, min_ETA_tasks[Subkernel_list_len], tie_list_num = 0, doubletie_list_num = 0; 
	long double min_ETA = DBL_MAX;
	if(!Subkernel_list_len) error("SubkernelSelect: Gave 0 subkernel len with list = %p\n", Subkernel_list);
	for (sk_idx = 0; task_idx < Subkernel_list_len; task_idx++)if(!Subkernel_list[sk_idx]->launched){
		curr_sk = Subkernel_list[sk_idx];
		long double tmp_ETA = 0; 
		for (int j = 0; j < curr_sk->TileNum; j++){
			long double block_ETA = 0; 
			if ( RONLY == curr_sk->TileList[j]->WRP || WR == curr_sk->TileList[j]->WRP){
				block_ETA = curr_sk->TileList[j]->ETA_get(dev_id);
				if(-42 == block_ETA) block_ETA = curr_sk->TileList[j]->ETA_fetch_estimate(dev_id);
			}
			tmp_ETA = std::max(block_ETA, tmp_ETA);
		}
		if(tmp_ETA < min_ETA){
		//if(abs(tmp_ETA - min_ETA)/abs(tmp_ETA-csecond()) > NORMALIZE_NEAR_SPLIT_LIMIT && tmp_ETA < min_ETA){
			min_ETA = tmp_ETA;
			min_ETA_tasks[0] = task_idx;
			tie_list_num = 1; 
		}
		else if(tmp_ETA == min_ETA){
		//else if(abs(tmp_ETA - min_ETA)/abs(tmp_ETA-csecond()) <= NORMALIZE_NEAR_SPLIT_LIMIT){
			min_ETA_tasks[tie_list_num++] = task_idx;
		}
	}
	int most_fired_sks = -1, potential_tied_sks[tie_list_num];
	if (tie_list_num){
		potential_tied_sks[0] = min_ETA_tasks[0];
		doubletie_list_num = 1; 
	}
	else error("SubkernelSelect\n No sk matched search condition\n");
	for (int ctr = 0; ctr < tie_list_num; ctr++){
		curr_sk = Subkernel_list[min_ETA_tasks[ctr]];
		int tmp_fired_sks = 0; 
		for (int j = 0; j < curr_sk->TileNum; j++){
			if ( WR_LAZY == curr_sk->TileList[j]->WRP || WR == curr_sk->TileList[j]->WRP
				|| W_REDUCE == curr_sk->TileList[j]->WRP || WONLY == curr_sk->TileList[j]->WRP){
				tmp_fired_sks = curr_sk->TileList[j]->fired_times; 
			}
		}
		if(tmp_fired_sks > most_fired_sks){
			most_fired_sks = tmp_fired_sks;
			potential_tied_sks[0] = min_ETA_tasks[ctr];
			doubletie_list_num = 1; 
		}
		//else if(tmp_fired_sks == most_fired_sks){
		//else if(abs(tmp_ETA - min_ETA)/abs(tmp_ETA-csecond()) <= NORMALIZE_NEAR_SPLIT_LIMIT){
		//	potential_tied_sks[doubletie_list_num++] = min_ETA_tasks[ctr];
		//}
	}
	int selected_sk_idx = (doubletie_list_num)? 
		potential_tied_sks[int(rand() % doubletie_list_num)] : doubletie_list_num; 
	Subkernel_list[selected_sk_idx]->prepare_launch(dev_id);
	return Subkernel_list[selected_sk_idx];
}
*/

void ATC::optimize_tasks_ETA(){
	if(strcmp(FETCH_ROUTING, "CHAIN_FETCH_QUEUE_WORKLOAD")) return optimize_tasks_MinFetchNum();
	long int comp_task_ctr = 0, comp_task_perdev[active_unit_num] = {0};
	long int size = T*T*elemSize;
	task_num = 0;
	int comp_task_fired[comp_task_num] = {0}, 
		comp_task_order[active_unit_num][comp_task_num] = {0};
	while (comp_task_ctr < comp_task_num){
		for(int dev_idx = 0; dev_idx < active_unit_num; dev_idx++){
			if(comp_task_perdev[dev_idx] == comp_task_per_unit_num[dev_idx]) continue;
			int dev_id = active_unit_id_list[dev_idx];
			long double min_ETA = DBL_MAX;
			int min_ETA_tasks[comp_task_per_unit_num[dev_idx]], tie_list_num = 0;
			for(int comp_dev_idx = 0; comp_dev_idx < comp_task_per_unit_num[dev_idx]; comp_dev_idx++){
				long comp_task_cand = comp_task_per_unit_list[dev_idx][comp_dev_idx];
				if(comp_task_fired[comp_task_cand]) continue;
				long comp_task_Cidx = comp_task_cand/Grid_K;
				int im = comp_task_Cidx/Grid_N, in = comp_task_Cidx%Grid_N, ik = comp_task_cand%Grid_K;
				if(ik != 0 && !comp_task_fired[comp_task_Cidx*Grid_K]) continue;
				int k_last_rd_flag = 1;
				if(ik == Grid_K - 1) for (int ctr = 0; ctr < Grid_K - 1; ctr++) 
					if(!comp_task_fired[comp_task_Cidx*Grid_K + ctr]) k_last_rd_flag = 0; 
				if(!k_last_rd_flag) continue;
				LinkRoute_p A_tile_route = new LinkRoute(), 
					B_tile_route = new LinkRoute(), C_tile_route = new LinkRoute();
				long double temp_ETA = 0, temp_t;
				if(A_tile_loc_map[im][ik][dev_id] && A_tile_loc_map[im][ik][dev_id]!= 42){
					A_tile_loc_map[im][ik][dev_id] = 2; 
					temp_t = A_tile_route->optimize(A_tile_loc_map[im][ik], size, 0);
				}
				else temp_t = A_tile_ETA[im][ik][dev_id]; 
				temp_ETA = std::max(temp_t, temp_ETA);
				if(B_tile_loc_map[ik][in][dev_id] && B_tile_loc_map[ik][in][dev_id]!= 42){
					B_tile_loc_map[ik][in][dev_id] = 2; 
					temp_t = B_tile_route->optimize(B_tile_loc_map[ik][in], size, 0);
				}
				else temp_t = B_tile_ETA[ik][in][dev_id]; 
				temp_ETA = std::max(temp_t, temp_ETA);
				if(!strcmp(OUTPUT_ALGO_MODE, "ALGO_WR") && 
					C_tile_loc_map[im][in][dev_id] && C_tile_loc_map[im][in][dev_id]!= 42){
					C_tile_loc_map[im][in][dev_id] = 2; 
					temp_t = C_tile_route->optimize(C_tile_loc_map[im][in], size, 0);
				}
				else temp_t = C_tile_ETA[im][in][dev_id]; 
				temp_ETA = std::max(temp_t, temp_ETA);
				if(temp_ETA < min_ETA){
					min_ETA = temp_ETA;
					min_ETA_tasks[0] = comp_task_cand;
					tie_list_num = 1; 
				}
				else if (temp_ETA == min_ETA)
					min_ETA_tasks[tie_list_num++] = comp_task_cand;
				delete A_tile_route;
				delete B_tile_route;
				delete C_tile_route;
			}
			int selected_task_idx = min_ETA_tasks[int(rand() % tie_list_num)]; 
			decompose_comp_task_fetch(selected_task_idx, dev_idx);
#ifdef SUBKERNELS_FIRE_LAZY
			;
#else
			decompose_comp_task_run(selected_task_idx, dev_idx);
#ifdef ENABLE_SEND_RECV_OVERLAP
			decompose_comp_task_wb(selected_task_idx, dev_idx);
#endif
#endif
			comp_task_fired[selected_task_idx] = 1;
			comp_task_order[dev_idx][comp_task_perdev[dev_idx]] = selected_task_idx;
			comp_task_perdev[dev_idx]++;
			comp_task_ctr++;
		}
	}
	for(int dev_idx = 0; dev_idx < active_unit_num; dev_idx++) 
		for(int idx = 0; idx < comp_task_perdev[dev_idx]; idx++){
#ifdef SUBKERNELS_FIRE_LAZY
			decompose_comp_task_run(comp_task_order[dev_idx][idx], dev_idx);
			decompose_comp_task_wb(comp_task_order[dev_idx][idx], dev_idx);
#else
#ifndef ENABLE_SEND_RECV_OVERLAP
			decompose_comp_task_wb(comp_task_order[dev_idx][idx], dev_idx);
#endif
#endif
		}
}

void ATC::optimize_tasks_ETA_plus_MinPendingOps(){
	if(strcmp(FETCH_ROUTING, "CHAIN_FETCH_QUEUE_WORKLOAD")) return optimize_tasks_MinFetchNum();
	long int comp_task_ctr = 0, comp_task_perdev[active_unit_num] = {0};
	long int size = T*T*elemSize;
	task_num = 0;
	int comp_task_fired[comp_task_num] = {0}, 
		comp_task_order[active_unit_num][comp_task_num] = {0},  
		comp_Cops_fired[comp_task_num/Grid_K] = {0};
	while (comp_task_ctr < comp_task_num){
		for(int dev_idx = 0; dev_idx < active_unit_num; dev_idx++){
			if(comp_task_perdev[dev_idx] == comp_task_per_unit_num[dev_idx]) continue;
			int dev_id = active_unit_id_list[dev_idx];
			long double min_ETA = DBL_MAX;
			int min_ETA_tasks[comp_task_per_unit_num[dev_idx]], tie_list_num = 0;
			for(int comp_dev_idx = 0; comp_dev_idx < comp_task_per_unit_num[dev_idx]; comp_dev_idx++){
				long comp_task_cand = comp_task_per_unit_list[dev_idx][comp_dev_idx];
				if(comp_task_fired[comp_task_cand]) continue;
				long comp_task_Cidx = comp_task_cand/Grid_K;
				int im = comp_task_Cidx/Grid_N, in = comp_task_Cidx%Grid_N, ik = comp_task_cand%Grid_K;
				if(ik != 0 && !comp_task_fired[comp_task_Cidx*Grid_K]) continue;
				int k_last_rd_flag = 1;
				if(ik == Grid_K - 1) for (int ctr = 0; ctr < Grid_K - 1; ctr++) 
					if(!comp_task_fired[comp_task_Cidx*Grid_K + ctr]) k_last_rd_flag = 0; 
				if(!k_last_rd_flag) continue;
				LinkRoute_p A_tile_route = new LinkRoute(), 
					B_tile_route = new LinkRoute(), C_tile_route = new LinkRoute();
				long double temp_ETA = 0, temp_t;
				if(A_tile_loc_map[im][ik][dev_id] && A_tile_loc_map[im][ik][dev_id]!= 42)
					temp_t = A_tile_route->optimize(A_tile_loc_map[im][ik], size, 0);
				else temp_t = A_tile_ETA[im][ik][dev_id]; 
				temp_ETA = std::max(temp_t, temp_ETA);
				if(B_tile_loc_map[ik][in][dev_id] && B_tile_loc_map[ik][in][dev_id]!= 42)
					temp_t = B_tile_route->optimize(B_tile_loc_map[ik][in], size, 0);
				else temp_t = B_tile_ETA[ik][in][dev_id]; 
				temp_ETA = std::max(temp_t, temp_ETA);
				if(!strcmp(OUTPUT_ALGO_MODE, "ALGO_WR") && 
					C_tile_loc_map[im][in][dev_id] && C_tile_loc_map[im][in][dev_id]!= 42){
					C_tile_loc_map[im][in][dev_id] = 2; 
					temp_t = C_tile_route->optimize(C_tile_loc_map[im][in], size, 0);
				}
				else temp_t = C_tile_ETA[im][in][dev_id]; 
				temp_ETA = std::max(temp_t, temp_ETA);
				if(temp_ETA < min_ETA){
					min_ETA = temp_ETA;
					min_ETA_tasks[0] = comp_task_cand;
					tie_list_num = 1; 
				}
				else if (temp_ETA == min_ETA)
					min_ETA_tasks[tie_list_num++] = comp_task_cand;
				delete A_tile_route;
				delete B_tile_route;
				delete C_tile_route;
			}
			int min_rem_ops = INT_MAX;
			int min_ETA_plus_minops_tasks[tie_list_num], tie_list_num_minops = 0;
			for(int comp_idx = 0; comp_idx < tie_list_num; comp_idx++){
				long comp_task_cand = min_ETA_tasks[comp_idx];
				int temp_remops = Grid_K - comp_Cops_fired[comp_task_cand/Grid_K];
				if(temp_remops < min_rem_ops){
					min_rem_ops = temp_remops;
					min_ETA_plus_minops_tasks[0] = comp_task_cand;
					tie_list_num_minops = 1; 
				}
				else if (temp_remops == min_rem_ops) 
					min_ETA_plus_minops_tasks[tie_list_num_minops++] = comp_task_cand;
			}
			int selected_task_idx = min_ETA_plus_minops_tasks[int(rand() % tie_list_num_minops)]; 
			decompose_comp_task_fetch(selected_task_idx, dev_idx);
#ifdef SUBKERNELS_FIRE_LAZY
			;
#else
			decompose_comp_task_run(selected_task_idx, dev_idx);
#ifdef ENABLE_SEND_RECV_OVERLAP
			decompose_comp_task_wb(selected_task_idx, dev_idx);
#endif
#endif
			comp_task_fired[selected_task_idx] = 1;
			comp_task_order[dev_idx][comp_task_perdev[dev_idx]] = selected_task_idx;
			comp_Cops_fired[selected_task_idx/Grid_K]++;
			comp_task_perdev[dev_idx]++;
			comp_task_ctr++;
		}
	}
	for(int dev_idx = 0; dev_idx < active_unit_num; dev_idx++) 
		for(int idx = 0; idx < comp_task_perdev[dev_idx]; idx++){
#ifdef SUBKERNELS_FIRE_LAZY
			decompose_comp_task_run(comp_task_order[dev_idx][idx], dev_idx);
			decompose_comp_task_wb(comp_task_order[dev_idx][idx], dev_idx);
#else
#ifndef ENABLE_SEND_RECV_OVERLAP
			decompose_comp_task_wb(comp_task_order[dev_idx][idx], dev_idx);
#endif
#endif
		}
}

/*

long double LinkRoute::optimize_reverse(void* transfer_tile_wrapped, int update_ETA_flag){
  DataTile_p transfer_tile = (DataTile_p) transfer_tile_wrapped;
  hop_num = 2;
  long double fire_t = csecond();
  int start_hop = -42, end_hop = -42;
  for(int ctr = 0; ctr < CHL_MEMLOCS; ctr++)
  {
	  if(transfer_tile->loc_map[ctr] == 42) start_hop = ctr;
	  if(transfer_tile->loc_map[ctr] == 0) end_hop = ctr;
  0
#ifdef ENABLE_TRANSFER_HOPS
  long double min_ETA = optimize_hop_route(transfer_tile_wrapped, 0, end_hop, start_hop);
#else
	hop_uid_list[0] = start_hop;
  hop_uid_list[1] = end_hop;
  long double temp_t = transfer_tile->size()/(1e9*get_edge_bw(hop_uid_list[1], hop_uid_list[0]));
  if (update_ETA_flag)
	recv_queues[hop_uid_list[1]][hop_uid_list[0]]->ETA_add_task(fire_t, temp_t);
  long double min_ETA = std::max(recv_queues[hop_uid_list[1]][hop_uid_list[0]]->ETA_get(), 
				fire_t) + temp_t;
#endif
  return min_ETA;

}

#ifdef HOP_FETCH_BANDWIDTH
// BW-based hop optimization
long double LinkRoute::optimize_hop_route(void* transfer_tile_wrapped, int update_ETA_flag, int dest_loc, int src_loc){
	DataTile_p transfer_tile = (DataTile_p) transfer_tile_wrapped;
	if (MAX_ALLOWED_HOPS > 1) error("LinkRoute::optimize_hop_route: Not implemented for MAX_ALLOWED_HOPS = %d\n", MAX_ALLOWED_HOPS);
	hop_uid_list[0] = src_loc;
	hop_num = 1;
	long double min_ETA = 0; 
	long double fire_t = csecond();
	int best_list[CHL_MEMLOCS], tie_list_num = 0; 
	double hop_bw_best = get_edge_bw(dest_loc,src_loc);
	for(int uidx = 0; uidx < CHL_MEMLOCS; uidx++)
	  if (best_grid_edge_active[uidx][src_loc]!= -1 && best_grid_edge_active[dest_loc][uidx]!= -1){
		double hop_est_bw = (1 - HOP_PENALTY) * std::min(get_edge_bw(uidx,src_loc), 
		  get_edge_bw(dest_loc, uidx));
		if (hop_est_bw  > hop_bw_best){
		  hop_bw_best = hop_est_bw;
		  best_list[0] = uidx;
		  tie_list_num = 1; 
		}
		else if (hop_est_bw  == hop_bw_best){
		  best_list[tie_list_num++] = uidx;
		}
	  }
	if (tie_list_num) hop_uid_list[hop_num++] = best_list[int(rand() % tie_list_num)];
	hop_uid_list[hop_num++] = dest_loc;
#ifdef SDEBUG
	if(hop_num > 2) fprintf(stderr, "Optimizing transfer %d -> %d : Route = %s\n", 
	  src_loc, dest_loc, printlist<int>(hop_uid_list, hop_num));
#endif
	/// TODO: NOTE - always adding ETA to recv_queue instead of wb. 
	if (hop_num == 2){
	  long double temp_t = transfer_tile->size()/(1e9*get_edge_bw(hop_uid_list[1], hop_uid_list[0]));
	  if (update_ETA_flag)
		recv_queues[(hop_uid_list[1])][(hop_uid_list[0])]->ETA_add_task(fire_t, temp_t);
	  min_ETA = std::max(recv_queues[(hop_uid_list[1])][(hop_uid_list[0])]->ETA_get(), 
				fire_t) + temp_t;
	}
	else{
	   long double max_t = -1, total_t = 0, temp_t;
		for (int ctr = 0; ctr < hop_num - 1; ctr++){
		  temp_t = transfer_tile->size()/(1e9*get_edge_bw(hop_uid_list[ctr+1], hop_uid_list[ctr]));
		  if (temp_t > max_t) max_t = temp_t;
		  total_t += temp_t;
		}
		temp_t = max_t + (total_t - max_t)/STREAMING_BUFFER_OVERLAP;
		//fprintf(stderr,"Checking location list[%s]: temp_t = %lf\n", printlist(templist,hop_num), temp_t);
		long double queue_ETA; 
		for(int ctr = 0; ctr < hop_num - 1; ctr++){
			queue_ETA = std::max(recv_queues[(hop_uid_list[ctr+1])][(hop_uid_list[ctr])]->ETA_get(), fire_t) + temp_t;
			//fprintf(stderr,"queue_ETA [%d -> %d]= %lf (temp_t = %lf)\n", prev), templist[ctr], queue_ETA);
			if(min_ETA < queue_ETA) min_ETA = queue_ETA;
			//fprintf(stderr,"Checking location list[%s]: queue_ETA = %lf\n", printlist(templist,hop_num), queue_ETA);
		}
		for(int ctr = 0; ctr < hop_num - 1; ctr++){
		  if (update_ETA_flag) recv_queues[(hop_uid_list[ctr+1])][(hop_uid_list[ctr])]->ETA_set(min_ETA);
	  }
	}
	return min_ETA;
}
#endif

#ifdef HOP_FETCH_QUEUE_WORKLOAD
// ETA-based hop optimization
long double LinkRoute::optimize_hop_route(void* transfer_tile_wrapped, int update_ETA_flag, int dest_loc, int src_loc){
	DataTile_p transfer_tile = (DataTile_p) transfer_tile_wrapped;
	if (MAX_ALLOWED_HOPS > 1) error("LinkRoute::optimize_hop_route: Not implemented for MAX_ALLOWED_HOPS = %d\n", MAX_ALLOWED_HOPS);
	hop_uid_list[0] = src_loc;
	hop_num = 1;
	int best_list[CHL_MEMLOCS], tie_list_num = 0; 
	long double fire_t = csecond();
	double tile_t = transfer_tile->size()/(1e9*get_edge_bw(dest_loc, src_loc));
	long double min_ETA = std::max(recv_queues[(dest_loc)][(src_loc)]->ETA_get(), fire_t) + tile_t;
	for(int uidx = 0; uidx < CHL_MEMLOCS; uidx++)
	  if (best_grid_edge_active[uidx][src_loc]!= -1 && best_grid_edge_active[dest_loc][uidx]!= -1){
		long double temp_t = (1 + HOP_PENALTY) *std::max(transfer_tile->size()/(1e9*get_edge_bw(uidx, src_loc)),
		  transfer_tile->size()/(1e9*get_edge_bw(dest_loc, uidx)));
		long double total_t = (1 + HOP_PENALTY) * transfer_tile->size()/(1e9*get_edge_bw(uidx, src_loc)) +
		  transfer_tile->size()/(1e9*get_edge_bw(dest_loc, uidx)); 
		temp_t = temp_t + (total_t - temp_t)/STREAMING_BUFFER_OVERLAP;
		//fprintf(stderr,"Checking location list[%s]: temp_t = %lf\n", printlist(templist,hop_num), temp_t);
		long double queue_ETA = std::max(std::max(recv_queues[(dest_loc)][uidx]->ETA_get(), fire_t),
		  std::max(recv_queues[uidx][(src_loc)]->ETA_get(), fire_t)) + temp_t; 
		if(queue_ETA < min_ETA && BANDWIDTH_DIFFERENCE_CUTTOF_RATIO*tile_t >= temp_t){
		  min_ETA = queue_ETA;
		  best_list[0] = uidx;
		  tie_list_num = 1; 
		}
		else if (queue_ETA == min_ETA && BANDWIDTH_DIFFERENCE_CUTTOF_RATIO*tile_t >= temp_t){
		  best_list[tie_list_num++] = uidx;
		}
	  }
	if (tie_list_num) hop_uid_list[hop_num++] = best_list[int(rand() % tie_list_num)];
	hop_uid_list[hop_num++] = dest_loc;
#ifdef SDEBUG
	if(hop_num > 2) fprintf(stderr, "Optimizing transfer %d -> %d : Route = %s\n", 
	  src_loc, dest_loc, printlist<int>(hop_uid_list, hop_num));
#endif
	/// TODO: NOTE - always adding ETA to recv_queue instead of wb. 
	for(int ctr = 0; ctr < hop_num - 1; ctr++)
	  if (update_ETA_flag) recv_queues[(hop_uid_list[ctr+1])][(hop_uid_list[ctr])]->ETA_set(min_ETA);
	return min_ETA;
}
#endif

#ifdef HOP_FETCH_BW_PLUS_ETA
// BW + ETA-based hop optimization
long double LinkRoute::optimize_hop_route(void* transfer_tile_wrapped, int update_ETA_flag, int dest_loc, int src_loc){
	DataTile_p transfer_tile = (DataTile_p) transfer_tile_wrapped;
	if (MAX_ALLOWED_HOPS > 1) error("LinkRoute::optimize_hop_route: Not implemented for MAX_ALLOWED_HOPS = %d\n", MAX_ALLOWED_HOPS);
	hop_uid_list[0] = src_loc;
	hop_num = 1;
	int best_list[CHL_MEMLOCS], tie_list_num = 0; 
	long double fire_t = csecond();
	double hop_bw_best = get_edge_bw(dest_loc,src_loc);
	double tile_t = transfer_tile->size()/(1e9*get_edge_bw(dest_loc, src_loc));
	long double min_ETA = std::max(recv_queues[(dest_loc)][(src_loc)]->ETA_get(), fire_t) + tile_t;
	for(int uidx = 0; uidx < CHL_MEMLOCS; uidx++)
	  if (best_grid_edge_active[uidx][src_loc]!= -1 && best_grid_edge_active[dest_loc][uidx]!= -1){
		long double temp_t = (1 + HOP_PENALTY) *std::max(transfer_tile->size()/(1e9*get_edge_bw(uidx, src_loc)),
		  transfer_tile->size()/(1e9*get_edge_bw(dest_loc, uidx)));
		long double total_t = (1 + HOP_PENALTY) * transfer_tile->size()/(1e9*get_edge_bw(uidx, src_loc)) +
		  transfer_tile->size()/(1e9*get_edge_bw(dest_loc, uidx)); 
		temp_t = temp_t + (total_t - temp_t)/STREAMING_BUFFER_OVERLAP;
		//fprintf(stderr,"Checking location list[%s]: temp_t = %lf\n", printlist(templist,hop_num), temp_t);
		long double queue_ETA = std::max(std::max(recv_queues[(dest_loc)][uidx]->ETA_get(), fire_t),
		  std::max(recv_queues[uidx][(src_loc)]->ETA_get(), fire_t)) + temp_t;
		double hop_est_bw = (1 - HOP_PENALTY) * std::min(get_edge_bw(uidx,src_loc), 
		  get_edge_bw(dest_loc, uidx));
		if(hop_est_bw > hop_bw_best || ( hop_est_bw == hop_bw_best && 
		  queue_ETA < min_ETA)){
		  min_ETA = queue_ETA;
		  hop_bw_best = hop_est_bw; 
		  best_list[0] = uidx;
		  tie_list_num = 1; 
		}
		else if (hop_est_bw == hop_bw_best && queue_ETA == min_ETA){
		  best_list[tie_list_num++] = uidx;
		}
	  }
	if (tie_list_num) hop_uid_list[hop_num++] = best_list[int(rand() % tie_list_num)];
	hop_uid_list[hop_num++] = dest_loc;
//#ifdef SDEBUG
	if(hop_num > 2) fprintf(stderr, "Optimizing transfer %d -> %d : Route = %s\n", 
	  src_loc, dest_loc, printlist<int>(hop_uid_list, hop_num));
//#endif
	/// TODO: NOTE - always adding ETA to recv_queue instead of wb. 
	for(int ctr = 0; ctr < hop_num - 1; ctr++)
	  if (update_ETA_flag) recv_queues[(hop_uid_list[ctr+1])][(hop_uid_list[ctr])]->ETA_set(min_ETA);
	return min_ETA;
}
#endif

*/
