
#include <iostream>
#include <cfloat>

#include "chl_smart_wrappers.hpp"
#include "Autotuner.hpp"
#include "DataTile.hpp"
#include "Resource_manager.hpp"

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
  
}

// Naive fetch from initial data loc every time
// Similar to cuBLASXt
#ifdef P2P_FETCH_FROM_INIT 
long double LinkRoute::optimize(void* transfer_tile_wrapped, int update_ETA_flag){
#ifdef DEBUG
	fprintf(stderr, "|-----> LinkRoute::optimize()\n");
#endif
  DataTile_p transfer_tile = (DataTile_p) transfer_tile_wrapped;
  hop_num = 2;
  int start_hop = -42, end_hop = -42;
  for(int ctr = 0; ctr < CHL_MEMLOCS; ctr++)
  {
      if(transfer_tile->loc_map[ctr] == 0) start_hop = ctr;
      if(transfer_tile->loc_map[ctr] == 2) end_hop = ctr;
  }
  hop_uid_list[0] = start_hop;
  hop_uid_list[1] = end_hop;
#ifdef DEBUG
	fprintf(stderr, "<-----|\n");
#endif
  return 0;
}
#endif

// Outdated fetch with preference to GPU tiles but with simple serial search for src
// Similar to BLASX behaviour when its assumed topology does not fit to the interconnect
#ifdef P2P_FETCH_FROM_GPU_SERIAL 
long double LinkRoute::optimize(void* transfer_tile_wrapped, int update_ETA_flag){
  DataTile_p transfer_tile = (DataTile_p) transfer_tile_wrapped;
  hop_num = 2;
  int start_hop = -42, end_hop = -42;
  for(int ctr = 0; ctr < LOC_NUM; ctr++)
  {
      if(transfer_tile->loc_map[ctr] == 0 || transfer_tile->loc_map[ctr] == 42) start_hop = ctr;
      if(transfer_tile->loc_map[ctr] == 2) end_hop = ctr;
      if(start_hop!= -42 && end_hop != -42) break;
    }
  hop_uid_list[0] = start_hop;
  hop_uid_list[1] = end_hop;
  return 0;
}
#endif

// Fetch selection based on 'distance' from available sources (if multiple exist)
// Similar to XKBLAS and PARALiA 1.5
#ifdef P2P_FETCH_FROM_GPU_DISTANCE
long double LinkRoute::optimize(void* transfer_tile_wrapped, int update_ETA_flag){
    DataTile_p transfer_tile = (DataTile_p) transfer_tile_wrapped;
    hop_num = 2;
    int start_hop = -42, end_hop = -42;
    for(int ctr = 0; ctr < LOC_NUM; ctr++) if(transfer_tile->loc_map[ctr] == 2) end_hop = ctr;
    hop_uid_list[1] = end_hop;

    int end_hop_idx = (end_hop);    
    int pos_max = LOC_NUM;
    double link_bw_max = 0;
    for (int pos =0; pos < LOC_NUM; pos++) if (transfer_tile->loc_map[pos] == 0 || transfer_tile->loc_map[pos] == 42){
        double current_link_bw = get_edge_bw(end_hop_idx, pos);
        if (current_link_bw > link_bw_max){
          link_bw_max = current_link_bw;
          pos_max = pos;
        }
    }
    hop_uid_list[0] = pos_max;
    return 0;
}
#endif

// Fetch selection based on 'distance' from available sources (if multiple exist) 
// Exactly like PARALiA 1.5 
#ifdef P2P_FETCH_FROM_GPU_DISTANCE_PLUS
long double LinkRoute::optimize(void* transfer_tile_wrapped, int update_ETA_flag){
    DataTile_p transfer_tile = (DataTile_p) transfer_tile_wrapped;
    hop_num = 2;
    int start_hop = -42, end_hop = -42;
    for(int ctr = 0; ctr < LOC_NUM; ctr++) if(transfer_tile->loc_map[ctr] == 2) end_hop = ctr;
    hop_uid_list[1] = end_hop;

    int end_hop_idx = (end_hop);    
    int pos_max = LOC_NUM;
    double link_bw_max = 0;
    for (int pos =0; pos < LOC_NUM; pos++){
        if (pos == end_hop_idx || transfer_tile->StoreBlock[pos] == NULL) {
      //if (StoreBlock[pos]!= NULL)
      //  error("Tile2D(%d)::getClosestReadLoc(%d): Should not be called, Tile already available in %d.\n",  id, end_hop, end_hop);
      continue;
    }
    //StoreBlock[pos]->update_state(false);
    state temp = transfer_tile->StoreBlock[pos]->State;
    if (temp == AVAILABLE || temp == SHARABLE || temp == NATIVE ){
      event_status block_status = transfer_tile->StoreBlock[pos]->Available->query_status();
      if(block_status == COMPLETE || block_status == CHECKED || block_status == RECORDED){
        double current_link_bw = get_edge_bw(end_hop_idx, pos);
        if (block_status == RECORDED) current_link_bw-=current_link_bw*0.1; //FETCH_UNAVAILABLE_PENALTY=0.1
        if (current_link_bw > link_bw_max){
          link_bw_max = current_link_bw;
          pos_max = pos;
        }
        else if (current_link_bw == link_bw_max &&
        final_estimated_linkmap->link_uses[end_hop_idx][pos] < final_estimated_linkmap->link_uses[end_hop_idx][pos_max]){
          link_bw_max = current_link_bw;
          pos_max = pos;
        }
      }
    }
  }
#ifdef DEBUG
  fprintf(stderr, "|-----> LinkRoute::optimize(DataTile[%d:%d,%d]): Selecting src = %d for dest = %d\n", transfer_tile->id, transfer_tile->dim1, transfer_tile->dim2, pos_max, end_hop);
#endif
/*
  if (pos_max >= LOC_NUM) error("LinkRoute::optimize(DataTile[%d:%d,%d]): No location found for tile - bug\n", transfer_tile->id, transfer_tile->dim1, transfer_tile->dim2);
  //Global_Cache[pos_max]->lock();
  CBlock_p temp_outblock = transfer_tile->StoreBlock[pos_max];
  if(temp_outblock != NULL){
    temp_outblock->lock();
    //temp_outblock->update_state(true);
    state temp = temp_outblock->State;
    event_status block_status = temp_outblock->Available->query_status();
    if ((temp == AVAILABLE || temp == SHARABLE || temp == NATIVE) &&
    (block_status == COMPLETE || block_status == CHECKED|| block_status == RECORDED)){
      temp_outblock->add_reader(true);
      //Global_Cache[pos_max]->unlock();
      temp_outblock->unlock();
      #ifdef DDEBUG
        lprintf(lvl-1, "<-----|\n");
      #endif
      final_estimated_linkmap->link_uses[end_hop_idx][pos_max]++;
      hop_uid_list[0] = pos_max);
      return; 
    }
    else error("Tile1D(%d)::getClosestReadLoc(%d): pos_max = %d selected,\
      but something changed after locking its cache...fixme\n", transfer_tile->id, end_hop, pos_max);
  }
  else error("Tile1D(%d)::getClosestReadLoc(%d): pos_max = %d selected,\
    but StoreBlock[pos_max] was NULL after locking its cache...fixme\n", transfer_tile->id, end_hop, pos_max);
    */
   hop_uid_list[0] = pos_max);
  return 0;
}
#endif

// Naive, for comparison reasons mainly
#ifdef CHAIN_FETCH_SERIAL
long double LinkRoute::optimize(void* transfer_tile_wrapped, int update_ETA_flag){
    DataTile_p transfer_tile = (DataTile_p) transfer_tile_wrapped;
    hop_num = 1;
    for(int ctr = 0; ctr < CHL_MEMLOCS; ctr++){
      if(transfer_tile->loc_map[ctr] == 0) hop_uid_list[0] = ctr;
      else if(transfer_tile->loc_map[ctr] == 1 || transfer_tile->loc_map[ctr] == 2)
        hop_uid_list[hop_num++] = ctr;
    }
    return 0;
}
#endif

#ifdef CHAIN_FETCH_RANDOM
long double LinkRoute::optimize(void* transfer_tile_wrapped, int update_ETA_flag){
    DataTile_p transfer_tile = (DataTile_p) transfer_tile_wrapped;
    hop_num = 0;
    int loc_list[CHL_MEMLOCS];
    for(int ctr = 0; ctr < CHL_MEMLOCS; ctr++){
      if(transfer_tile->loc_map[ctr] == 0) hop_uid_list[0] = ctr;
      else if(transfer_tile->loc_map[ctr] == 1 || transfer_tile->loc_map[ctr] == 2)
        loc_list[hop_num++] = ctr;
    } 
    int start_idx = int(rand() % hop_num); 
    int hop_ctr = 1;
    for(int ctr = start_idx; ctr < hop_num; ctr++) hop_uid_list[hop_ctr++] = loc_list[ctr];
    for(int ctr = 0; ctr < start_idx; ctr++) hop_uid_list[hop_ctr++] = loc_list[ctr];
    hop_num++;
    return 0;
}
#endif

#include <algorithm>
#include <iostream>
#include <list>

#ifdef CHAIN_FETCH_TIME
long double LinkRoute::optimize(void* transfer_tile_wrapped, int update_ETA_flag){
    DataTile_p transfer_tile = (DataTile_p) transfer_tile_wrapped;
    hop_num = 0;
    std::list<int> loc_list;
    int tmp_hop = -42; 
    double fire_est = csecond();
    for(int ctr = 0; ctr < CHL_MEMLOCS; ctr++){
      if(transfer_tile->loc_map[ctr] == 0) hop_uid_list[0] = ctr;
      else if(transfer_tile->loc_map[ctr] == 1 || transfer_tile->loc_map[ctr] == 2){
        tmp_hop = ctr;
        loc_list.push_back(tmp_hop);
        hop_num++;
      }
    }
    double best_t = 1e9; 
    if (hop_num == 1){
      hop_uid_list[1] = tmp_hop;
      best_t = transfer_tile->size()/(1e9*get_edge_bw(hop_uid_list[1], hop_uid_list[0]));
    }
    else{
      int best_list[factorial(hop_num)][hop_num]; 
      int flag = 1, tie_list_num = 0;
      while (flag){
        double max_t = -1, total_t = 0, temp_t;
        int temp_ctr = 0, prev = (hop_uid_list[0]), templist[hop_num]; 
        for (int x : loc_list){
          templist[temp_ctr++] = x; 
          temp_t = transfer_tile->size()/(1e9*get_edge_bw(x, prev));
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
    return best_t;
}
#endif

#ifdef CHAIN_FETCH_QUEUE_WORKLOAD
long double LinkRoute::optimize(void* transfer_tile_wrapped, int update_ETA_flag){
    DataTile_p transfer_tile = (DataTile_p) transfer_tile_wrapped;
    long double min_ETA = DBL_MAX, tile_t = DBL_MAX/100, fire_t = csecond();
    hop_num = 0;
    std::list<int> loc_list;
    int tmp_hop = -42; 
    for(int ctr = 0; ctr < CHL_MEMLOCS; ctr++){
      if(transfer_tile->loc_map[ctr] == 0) hop_uid_list[0] = ctr;
      else if(transfer_tile->loc_map[ctr] == 1 || transfer_tile->loc_map[ctr] == 2){
        tmp_hop = ctr;
        loc_list.push_back(tmp_hop);
        hop_num++;
      }
    }
    if (hop_num == 1){
      hop_uid_list[1] = tmp_hop;
#ifdef ENABLE_TRANSFER_HOPS
      min_ETA = optimize_hop_route(transfer_tile_wrapped, update_ETA_flag, hop_uid_list[1], hop_uid_list[0]);
#else
      long double temp_t = transfer_tile->size()/(1e9*get_edge_bw(hop_uid_list[1], hop_uid_list[0]));
      if (update_ETA_flag)
        recv_queues[(hop_uid_list[1])][(hop_uid_list[0])]->ETA_add_task(fire_t, temp_t);
      min_ETA = std::max(recv_queues[(hop_uid_list[1])][(hop_uid_list[0])]->ETA_get(), 
                fire_t) + temp_t;
      hop_num++;
#endif
    }
    else{
      int best_list[factorial(hop_num)][hop_num]; 
      int flag = 1, tie_list_num = 0;
      while (flag){
        long double max_t = -1, total_t = 0, temp_t;
        int temp_ctr = 0, prev = hop_uid_list[0], templist[hop_num]; 
        for (int x : loc_list){
          templist[temp_ctr++] = x; 
          temp_t = transfer_tile->size()/(1e9*get_edge_bw(x, prev));
          if (temp_t > max_t) max_t = temp_t;
          total_t += temp_t;
          prev = x;
        }
        temp_t = max_t + (total_t - max_t)/STREAMING_BUFFER_OVERLAP;
        //fprintf(stderr,"Checking location list[%s]: temp_t = %lf\n", printlist(templist,hop_num), temp_t);
        prev = hop_uid_list[0];
        long double temp_ETA = 0, queue_ETA; 
        for(int ctr = 0; ctr < hop_num; ctr++){
            queue_ETA = std::max(recv_queues[(templist[ctr])][(prev)]->ETA_get(), fire_t) + temp_t;
            //fprintf(stderr,"queue_ETA [%d -> %d]= %lf (temp_t = %lf)\n", prev), templist[ctr], queue_ETA);
            prev = templist[ctr];
            if(temp_ETA < queue_ETA) temp_ETA = queue_ETA;
            //fprintf(stderr,"Checking location list[%s]: queue_ETA = %lf\n", printlist(templist,hop_num), queue_ETA);
        }
        //fprintf(stderr,"Checking location list[%s]: temp_ETA = %lf\n", printlist(templist,hop_num), temp_ETA);
        if(temp_ETA < min_ETA && BANDWIDTH_DIFFERENCE_CUTTOF_RATIO*tile_t >= temp_t){
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
#ifdef SDEBUG
      fprintf(stderr,"DataTile[%d:%d,%d]: Selected route = %s from %d candidates with ETA = %llf\n", 
          transfer_tile->id, transfer_tile->GridId1, transfer_tile->GridId2,
          printlist(best_list[tie_list_num-1],hop_num), tie_list_num, min_ETA);
#endif
      for(int ctr = 0; ctr < hop_num; ctr++){
        hop_uid_list[ctr+1] = best_list[rand_tie_list][ctr];
        if (update_ETA_flag) recv_queues[(hop_uid_list[ctr+1])][(hop_uid_list[ctr])]->ETA_set(min_ETA);
      }
      hop_num++;
    }
    return min_ETA; 
}
#endif

long double LinkRoute::optimize_reverse(void* transfer_tile_wrapped, int update_ETA_flag){
  DataTile_p transfer_tile = (DataTile_p) transfer_tile_wrapped;
  hop_num = 2;
  long double fire_t = csecond();
  int start_hop = -42, end_hop = -42;
  for(int ctr = 0; ctr < CHL_MEMLOCS; ctr++)
  {
      if(transfer_tile->loc_map[ctr] == 42) start_hop = ctr;
      if(transfer_tile->loc_map[ctr] == 0) end_hop = ctr;
  }
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
