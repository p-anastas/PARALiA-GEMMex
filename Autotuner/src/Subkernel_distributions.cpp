///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief The possible subkernel distributions to different execution units.
///

#include "Autotuner.hpp"

#include <cmath>
/*
void CoCoDistributeTasksRoundRobin(ATC_p autotune_controller){
  #ifdef DEBUG
  	fprintf(stderr, "|-----> CoCoDistributeTasksRoundRobin(%p)\n", autotune_controller);
  #endif
  if (autotune_controller->subkernel_num < autotune_controller->active_unit_num){ // < or <= here?
    int pred_active_unit_num = autotune_controller->active_unit_num;
    autotune_controller->active_unit_num = autotune_controller->subkernel_num;
    warning("CoCoDistributeTasksRoundRobin: Problem with predicted active_unit_num(%d) < subkernel_num(%d) will be run with active_unit_num = %d\n",
    	pred_active_unit_num, autotune_controller->subkernel_num, autotune_controller->active_unit_num);
    for (int d = 0 ; d < autotune_controller->active_unit_num; d++){
      autotune_controller->Tasks_per_unit_num[d] = 1;
      autotune_controller->Tasks_per_unit_list[d][0] = d;
    }
  }
  else{
    int rem_dev = autotune_controller->subkernel_num;
    for (int d = 0 ; d < autotune_controller->active_unit_num; d++){
       autotune_controller->Tasks_per_unit_num[d] =
        (int) (1.0* autotune_controller->active_unit_score[d]* autotune_controller->subkernel_num);
       rem_dev-= autotune_controller->Tasks_per_unit_num[d];
    }
    while(rem_dev!= 0){
      for (int d = 0 ; d < autotune_controller->active_unit_num; d++){
         if(rem_dev!= 0){
           autotune_controller->Tasks_per_unit_num[d] += 1;
           rem_dev--;
         }
         else break;
      }
    }
    int total_sk_ctr = 0;
    short dev_sk_ctr_list[autotune_controller->active_unit_num];
    for(int devidx = 0; devidx < autotune_controller->active_unit_num; devidx++) dev_sk_ctr_list[devidx] = 0;
    while(total_sk_ctr<autotune_controller->subkernel_num){
      for(int devidx = 0; devidx < autotune_controller->active_unit_num; devidx++){
        if(total_sk_ctr == autotune_controller->subkernel_num) break;
        else if(dev_sk_ctr_list[devidx] == autotune_controller->Tasks_per_unit_num[devidx]) continue;
        else{
          autotune_controller->Tasks_per_unit_list[devidx][dev_sk_ctr_list[devidx]] = total_sk_ctr;
          dev_sk_ctr_list[devidx]++;
          total_sk_ctr++;
        }
      }
    }
  }
#ifdef PDEBUG
  fprintf(stderr, "CoCoDistributeTasksRoundRobin:\nDistributing %ld Tasks to %d devices\n",
    autotune_controller->subkernel_num, autotune_controller->active_unit_num);
  fprintf(stderr, "Device Ids : [ ");
  for (int i =0; i < autotune_controller->active_unit_num; i++) fprintf(stderr, "%d ", autotune_controller->active_unit_id_list[i]);
  lprintf(0, "]\n");
  fprintf(stderr, "Subker Num : [ ");
  for (int i =0; i < autotune_controller->active_unit_num; i++) fprintf(stderr, "%d ",
    autotune_controller->Tasks_per_unit_num[i]);
  lprintf(0, "]\n");
  for (int i =0; i < autotune_controller->active_unit_num; i++){
    fprintf(stderr, "Subker Id list for dev_id = %d: [ ", autotune_controller->active_unit_id_list[i]);
    for (int j =0; j < autotune_controller->Tasks_per_unit_num[i]; j++) fprintf(stderr, "%d ",
      autotune_controller->Tasks_per_unit_list[i][j]);
    lprintf(0, "]\n");
  }
#endif
#ifdef DEBUG
	fprintf(stderr, "<-----|\n");
#endif
}

void CoCoDistributeTasksNaive(ATC_p autotune_controller){
  #ifdef DEBUG
  	fprintf(stderr, "|-----> CoCoDistributeTasksNaive(%p)\n", autotune_controller);
  #endif
  if (autotune_controller->subkernel_num <= autotune_controller->active_unit_num){
    int pred_active_unit_num = autotune_controller->active_unit_num;
    autotune_controller->active_unit_num = autotune_controller->subkernel_num;
    warning("CoCoDistributeTasksNaive: Problem with predicted active_unit_num(%d) < subkernel_num(%d) will be run with active_unit_num = %d\n",
    	pred_active_unit_num, autotune_controller->subkernel_num, autotune_controller->active_unit_num);
    for (int d = 0 ; d < autotune_controller->active_unit_num; d++){
      autotune_controller->Tasks_per_unit_num[d] = 1;
      autotune_controller->Tasks_per_unit_list[d][0] = d;
    }
  }
  else{
    int total_sk_ctr = 0;
    int rem_dev = autotune_controller->subkernel_num;
    for (int d = 0 ; d < autotune_controller->active_unit_num; d++){
       autotune_controller->Tasks_per_unit_num[d] =
        (int) (1.0* autotune_controller->active_unit_score[d]* autotune_controller->subkernel_num);
       rem_dev-= autotune_controller->Tasks_per_unit_num[d];
    }
    while(rem_dev!= 0){
      for (int d = 0 ; d < autotune_controller->active_unit_num; d++){
         if(rem_dev!= 0){
           autotune_controller->Tasks_per_unit_num[d] += 1;
           rem_dev--;
         }
         else break;
      }
    }
    short dev_sk_ctr = 0, cur_dev_id_ctr = 0;
    while(total_sk_ctr<autotune_controller->subkernel_num && cur_dev_id_ctr < autotune_controller->active_unit_num){
      while(dev_sk_ctr == autotune_controller->Tasks_per_unit_num[cur_dev_id_ctr]){
        dev_sk_ctr = 0;
        cur_dev_id_ctr++;
      }
      autotune_controller->Tasks_per_unit_list[cur_dev_id_ctr][dev_sk_ctr] = total_sk_ctr;
      dev_sk_ctr++;
      total_sk_ctr++;
    }
  }
#ifdef PDEBUG
    fprintf(stderr, "CoCoDistributeTasksNaive:\nDistributing %ld Tasks to %d devices\n",
      autotune_controller->subkernel_num, autotune_controller->active_unit_num);
    fprintf(stderr, "Device Ids : [ ");
    for (int i =0; i < autotune_controller->active_unit_num; i++) fprintf(stderr, "%d ", autotune_controller->active_unit_id_list[i]);
    lprintf(0, "]\n");
    fprintf(stderr, "Subker Num : [ ");
    for (int i =0; i < autotune_controller->active_unit_num; i++) fprintf(stderr, "%d ",
      autotune_controller->Tasks_per_unit_num[i]);
    lprintf(0, "]\n");
    for (int i =0; i < autotune_controller->active_unit_num; i++){
      fprintf(stderr, "Subker Id list for dev_id = %d: [ ", autotune_controller->active_unit_id_list[i]);
      for (int j =0; j < autotune_controller->Tasks_per_unit_num[i]; j++) fprintf(stderr, "%d ",
        autotune_controller->Tasks_per_unit_list[i][j]);
      lprintf(0, "]\n");
    }
#endif
#ifdef DEBUG
	fprintf(stderr, "<-----|\n");
#endif
}

void CoCoDistributeTasksRoundRobinChunk(ATC_p autotune_controller,  int Chunk_size){
#ifdef DEBUG
  	fprintf(stderr, "|-----> CoCoDistributeTasksRoundRobinChunk(%p, %d)\n", autotune_controller, Chunk_size);
#endif
#ifdef PDEBUG
fprintf(stderr, "CoCoDistributeTasksRoundRobinChunk(%d): Devices = %d (scores = %s), sk_num = %d, sk_buckets = %d\n",
  Chunk_size, autotune_controller->active_unit_num, 
  printlist<double>(autotune_controller->active_unit_score, autotune_controller->active_unit_num),
  autotune_controller->subkernel_num, autotune_controller->subkernel_num/Chunk_size);
#endif
  if (autotune_controller->subkernel_num/Chunk_size <= autotune_controller->active_unit_num){
    int pred_active_unit_num = autotune_controller->active_unit_num;
    autotune_controller->active_unit_num = autotune_controller->subkernel_num/Chunk_size;
    warning("CoCoDistributeTasksRoundRobinChunk: Problem with predicted active_unit_num(%d) < subkernel_num(%d) will be run with active_unit_num = %d\n",
    	pred_active_unit_num, autotune_controller->subkernel_num, autotune_controller->active_unit_num);
    for (int d = 0 ; d < autotune_controller->active_unit_num; d++){
      autotune_controller->Tasks_per_unit_num[d] = Chunk_size;
      for (int idx3 = 0; idx3 < Chunk_size; idx3++)
        autotune_controller->Tasks_per_unit_list[d][idx3] = d*Chunk_size + idx3;
    }
  }
  else{
  for (int d = 0 ; d < autotune_controller->active_unit_num; d++){
    autotune_controller->Tasks_per_unit_num[d] = Chunk_size*
      (int) (1.0* autotune_controller->active_unit_score[d]* (autotune_controller->subkernel_num/Chunk_size));
  }
  int sks_accounted_for = 0;
  for (int d = 0 ; d < autotune_controller->active_unit_num; d++)
    sks_accounted_for += autotune_controller->Tasks_per_unit_num[d];
#ifdef PDEBUG
fprintf(stderr, "Assigned kernel num to devices kernels (first pass): %s\n",
  printlist<int>(autotune_controller->Tasks_per_unit_num, autotune_controller->active_unit_num));
#endif
  int sk_ctr = 0, dev_sk_ctr_list[autotune_controller->active_unit_num] = {0}, devidx = 0;
  for (int D1 = 0; D1 < autotune_controller->subkernel_num/Chunk_size; D1++){
    for (int D3 = 0; D3 < Chunk_size; D3++){
      int full_circle = autotune_controller->active_unit_num;
      while(dev_sk_ctr_list[devidx] == autotune_controller->Tasks_per_unit_num[devidx] && sk_ctr < sks_accounted_for){ 
        if(!full_circle) error("CoCoDistributeTasks2DBlockCyclic: would enter infinite loop due to wrong Tasks_per_unit_num, terminating\n");
        if(devidx == autotune_controller->active_unit_num - 1) devidx = 0;
        else devidx++;
        full_circle--;
      }
      autotune_controller->Tasks_per_unit_list[devidx][dev_sk_ctr_list[devidx]] = sk_ctr;
      if(sk_ctr >= sks_accounted_for) autotune_controller->Tasks_per_unit_num[devidx] ++;
      dev_sk_ctr_list[devidx]++;
#ifdef PDEBUG
      fprintf(stderr, "CoCoDistributeTasksRoundRobinChunk: sk_ctr[%d,%d] = %d, devidx = %d\n",
        D1, D3, sk_ctr, devidx);
#endif
      sk_ctr++;
    }
    if(devidx == autotune_controller->active_unit_num - 1) devidx = 0;
    else devidx++;
  }
  }
#ifdef PDEBUG
    fprintf(stderr, "CoCoDistributeTasksRoundRobinChunk:\nDistributing %ld Tasks to %d devices\n",
      autotune_controller->subkernel_num, autotune_controller->active_unit_num);
    fprintf(stderr, "Device Ids : [ ");
    for (int i =0; i < autotune_controller->active_unit_num; i++) fprintf(stderr, "%d ", autotune_controller->active_unit_id_list[i]);
    lprintf(0, "]\n");
    fprintf(stderr, "Subker Num : [ ");
    for (int i =0; i < autotune_controller->active_unit_num; i++) fprintf(stderr, "%d ",
      autotune_controller->Tasks_per_unit_num[i]);
    lprintf(0, "]\n");
    for (int i =0; i < autotune_controller->active_unit_num; i++){
      fprintf(stderr, "Subker Id list for dev_id = %d: [ ", autotune_controller->active_unit_id_list[i]);
      for (int j =0; j < autotune_controller->Tasks_per_unit_num[i]; j++) fprintf(stderr, "%d ",
        autotune_controller->Tasks_per_unit_list[i][j]);
      lprintf(0, "]\n");
    }
#endif
#ifdef DEBUG
	fprintf(stderr, "<-----|\n");
#endif
}

void CoCoDistributeTasksRoundRobinChunkReverse(ATC_p autotune_controller,  int Chunk_size){
  #ifdef DEBUG
  	fprintf(stderr, "|-----> CoCoDistributeTasksRoundRobinChunkReverse(%p, %d)\n", autotune_controller, Chunk_size);
  #endif
  if (autotune_controller->subkernel_num <= autotune_controller->active_unit_num){
    int pred_active_unit_num = autotune_controller->active_unit_num;
    autotune_controller->active_unit_num = autotune_controller->subkernel_num;
    warning("CoCoDistributeTasksRoundRobinChunkReverse: Problem with predicted active_unit_num(%d) < subkernel_num(%d) will be run with active_unit_num = %d\n",
    	pred_active_unit_num, autotune_controller->subkernel_num, autotune_controller->active_unit_num);
    for (int d = 0 ; d < autotune_controller->active_unit_num; d++){
      autotune_controller->Tasks_per_unit_num[d] = 1;
      autotune_controller->Tasks_per_unit_list[d][0] = d;
    }
  }
  else{
    int rem_dev = autotune_controller->subkernel_num;
    for (int d = 0 ; d < autotune_controller->active_unit_num; d++){
       autotune_controller->Tasks_per_unit_num[d] =
        (int) (1.0* autotune_controller->active_unit_score[d]* autotune_controller->subkernel_num);
       rem_dev-= autotune_controller->Tasks_per_unit_num[d];
    }
    while(rem_dev!= 0){
      for (int d = 0 ; d < autotune_controller->active_unit_num; d++){
         if(rem_dev!= 0){
           autotune_controller->Tasks_per_unit_num[d] += 1;
           rem_dev--;
         }
         else break;
      }
    }
    int total_sk_ctr = 0, total_sk_prev = 0;
    short dev_sk_ctr_list[autotune_controller->active_unit_num];
    for(int devidx = 0; devidx < autotune_controller->active_unit_num; devidx++) dev_sk_ctr_list[devidx] = 0;
    while(total_sk_ctr<autotune_controller->subkernel_num){

      for(int devidx = 0; devidx < autotune_controller->active_unit_num; devidx++){
        total_sk_prev = total_sk_ctr;
        if(total_sk_ctr == autotune_controller->subkernel_num) break;
        else if(dev_sk_ctr_list[devidx] == autotune_controller->Tasks_per_unit_num[devidx]) continue;
        else{
          autotune_controller->Tasks_per_unit_list[devidx][dev_sk_ctr_list[devidx]] = total_sk_ctr;
          dev_sk_ctr_list[devidx]++;
          total_sk_ctr++;
        }
        while(total_sk_ctr%Chunk_size!=0){
          if(total_sk_ctr == autotune_controller->subkernel_num) break;
          else if(dev_sk_ctr_list[devidx] == autotune_controller->Tasks_per_unit_num[devidx]) break;
          else{
            autotune_controller->Tasks_per_unit_list[devidx][dev_sk_ctr_list[devidx]] = total_sk_ctr;
            dev_sk_ctr_list[devidx]++;
            total_sk_ctr++;
          }
        }
        if (devidx%2 == 0){
          for(int local_ctr = dev_sk_ctr_list[devidx] - total_sk_ctr + total_sk_prev; local_ctr < dev_sk_ctr_list[devidx]; local_ctr++){
            if (local_ctr < dev_sk_ctr_list[devidx] - local_ctr - 1){
              int temp_sk_id = autotune_controller->Tasks_per_unit_list[devidx][local_ctr];
              autotune_controller->Tasks_per_unit_list[devidx][local_ctr] =
                autotune_controller->Tasks_per_unit_list[devidx]
                  [dev_sk_ctr_list[devidx] - local_ctr - 1];
              autotune_controller->Tasks_per_unit_list[devidx]
                [dev_sk_ctr_list[devidx] - local_ctr - 1] = temp_sk_id;
            }
            else break;
          }
        }
        if(total_sk_ctr == autotune_controller->subkernel_num) break;
      }
    }
  }
#ifdef PDEBUG
  fprintf(stderr, "CoCoDistributeTasksRoundRobinChunkReverse:\nDistributing %ld Tasks to %d devices\n",
    autotune_controller->subkernel_num, autotune_controller->active_unit_num);
  fprintf(stderr, "Device Ids : [ ");
  for (int i =0; i < autotune_controller->active_unit_num; i++) fprintf(stderr, "%d ", autotune_controller->active_unit_id_list[i]);
  lprintf(0, "]\n");
  fprintf(stderr, "Subker Num : [ ");
  for (int i =0; i < autotune_controller->active_unit_num; i++) fprintf(stderr, "%d ",
    autotune_controller->Tasks_per_unit_num[i]);
  lprintf(0, "]\n");
  for (int i =0; i < autotune_controller->active_unit_num; i++){
    fprintf(stderr, "Subker Id list for dev_id = %d: [ ", autotune_controller->active_unit_id_list[i]);
    for (int j =0; j < autotune_controller->Tasks_per_unit_num[i]; j++) fprintf(stderr, "%d ",
      autotune_controller->Tasks_per_unit_list[i][j]);
    lprintf(0, "]\n");
  }
#endif
#ifdef DEBUG
	fprintf(stderr, "<-----|\n");
#endif
}

void CoCoDistributeTasks2DBlockCyclic(ATC_p autotune_controller, int D1GridSz, int D2GridSz, int D3GridSz){
#ifdef DEBUG
  	fprintf(stderr, "|-----> CoCoDistributeTasks2DBlockCyclic(%p, %d, %d, %d)\n", autotune_controller, D1GridSz, D2GridSz, D3GridSz);
#endif
  if ((D2GridSz == D3GridSz) &&  (D2GridSz == 1)){
    warning("CoCoDistributeTasks2DBlockCyclic: D2GridSz==D3GridSz==1 -> using CoCoDistributeTasksRoundRobin\n");
    return CoCoDistributeTasksRoundRobin(autotune_controller);
  }

  // 2D Block cyclic
  int D1_parts = std::sqrt(autotune_controller->active_unit_num);
  int D2_parts = D1_parts;
  if (D1_parts ==0) { D2_parts = autotune_controller->active_unit_num; D1_parts = 1; }
  else {
    // find the most square decomposition of autotune_controller->active_unit_num in D1_parts x D2_parts
    int g;
    for (g = D1_parts+1; g>0; --g)
       if (autotune_controller->active_unit_num % g == 0) break;
    if (g==0) { D1_parts = autotune_controller->active_unit_num; D2_parts = 1; }
    //if (g==0) { D1_parts = 1; D2_parts = autotune_controller->active_unit_num; }
    else { D1_parts = g; D2_parts = autotune_controller->active_unit_num/g; }
  }
  //TODO: reverse layout
  //int tmp = D1_parts;
  //D1_parts = D2_parts;
  //D2_parts = tmp;
  if(D1GridSz < D1_parts || D2GridSz < D2_parts){
    warning("CoCoDistributeTasks2DBlockCyclic:\nGrid(%d,%d) smaller than {D1,D2}_parts = (%d,%d)\
    using CoCoDistributeTasksRoundRobinChunk instead\n", D1GridSz, D2GridSz, D1_parts, D2_parts);
    CoCoDistributeTasksRoundRobinChunk(autotune_controller, D3GridSz);
    return;
  }
  int D1GridSz_div = D1GridSz/D1_parts*D1_parts, D2GridSz_div = D2GridSz/D2_parts*D2_parts,
      D1GridSz_mod = D1GridSz%D1_parts, D2GridSz_mod = D2GridSz%D2_parts;
#ifdef PDEBUG
fprintf(stderr, "CoCoDistributeTasks2DBlockCyclic(%d, %d, %d): Devices = %d (scores = %s), D1_parts = %d, D2_parts = %d\n",
  D1GridSz, D2GridSz, D3GridSz, autotune_controller->active_unit_num, 
  printlist<double>(autotune_controller->active_unit_score, autotune_controller->active_unit_num),D1_parts, D2_parts);
#endif

  if ((D1GridSz*D2GridSz) < autotune_controller->active_unit_num){
    warning("CoCoDistributeTasks2DBlockCyclic: D1GridSz*D2GridSz(%d) < autotune_controller->active_unit_num(%d)\n", 
      D1GridSz*D2GridSz, autotune_controller->active_unit_num);
    int pred_active_unit_num = D1GridSz*D2GridSz;
    autotune_controller->active_unit_num = autotune_controller->subkernel_num;
    warning("CoCoDistributeTasks2DBlockCyclic: Problem with predicted active_unit_num(%d) < subkernel_num(%d) will be run with active_unit_num = %d\n",
    	pred_active_unit_num, autotune_controller->subkernel_num, autotune_controller->active_unit_num);
    for (int d = 0 ; d < autotune_controller->active_unit_num; d++){
      autotune_controller->Tasks_per_unit_num[d] = D3GridSz;
      for (int idx3 = 0; idx3 < D3GridSz; idx3++)
        autotune_controller->Tasks_per_unit_list[d][idx3] = d*D3GridSz + idx3;
    }
  }
  else{
  int sks_accounted_for = 0;
  for (int d = 0 ; d < autotune_controller->active_unit_num; d++){
     autotune_controller->Tasks_per_unit_num[d] = D3GridSz * (
      (int) (autotune_controller->active_unit_score[d]* D1GridSz_div*D2GridSz_div));
      /// TODO: this is a fix because int (some_double) does not work for all doubles as intended 
      /// Will disrupt non-homogeneous splits!
      sks_accounted_for+= autotune_controller->Tasks_per_unit_num[d]; 
  }
  if(!D1GridSz_mod && !D2GridSz_mod && sks_accounted_for < autotune_controller->subkernel_num){
    warning("CoCoDistributeTasks2DBlockCyclic: Questionable remainder from first pass %d / %d sub-kernels\n",
      autotune_controller->subkernel_num - sks_accounted_for, autotune_controller->subkernel_num);
    int buckets =  D1GridSz_div*D2GridSz_div, 
        buckets_rem = (autotune_controller->subkernel_num - sks_accounted_for)/D3GridSz,
        buckets_intended = buckets/autotune_controller->active_unit_num; 
    for (int d = 0 ; d < autotune_controller->active_unit_num; d++) 
      if(autotune_controller->Tasks_per_unit_num[d]/D3GridSz < buckets_intended && buckets_rem){
        autotune_controller->Tasks_per_unit_num[d]+= D3GridSz; 
        buckets_rem--;
      }
  }
#ifdef PDEBUG
fprintf(stderr, "Assigned kernel num to devices kernels (first pass): %s\n",
  printlist<int>(autotune_controller->Tasks_per_unit_num, autotune_controller->active_unit_num));
#endif
  int sk_ctr, dev_sk_ctr_list[autotune_controller->active_unit_num] = {0}, devidx = 0;
  for (int D1 = 0; D1 < D1GridSz_div; D1++)
    for (int D2 = 0; D2 < D2GridSz_div; D2++)
        for (int D3 = 0; D3 < D3GridSz; D3++){
          sk_ctr = D1*D2GridSz*D3GridSz + D2*D3GridSz+D3;
          devidx = D1/(D1GridSz/D1_parts)*D2_parts + D2/(D2GridSz/D2_parts);
#ifdef PDEBUG
          fprintf(stderr, "CoCoDistributeTasks2DBlockCyclic: sk_ctr[%d,%d,%d] = %d, devidx = %d\n",
            D1,D2,D3, sk_ctr, devidx);
#endif
          int full_circle = autotune_controller->active_unit_num; 
          while(dev_sk_ctr_list[devidx] == autotune_controller->Tasks_per_unit_num[devidx]){
            if(!full_circle) error("CoCoDistributeTasks2DBlockCyclic: would enter infinite loop due to wrong Tasks_per_unit_num, terminating\n");
            if(devidx == autotune_controller->active_unit_num - 1) devidx = 0;
            else devidx++;
            full_circle--; 
          }
          autotune_controller->Tasks_per_unit_list[devidx][dev_sk_ctr_list[devidx]] = sk_ctr;
          dev_sk_ctr_list[devidx]++;
        }

  devidx = 0;
  for (int D1 = 0; D1 < D1GridSz_div; D1++)
    for (int D2 = D2GridSz_div; D2 < D2GridSz_div + D2GridSz_mod; D2++){
      for (int D3 = 0; D3 < D3GridSz; D3++){
        sk_ctr = D1*D2GridSz*D3GridSz + D2*D3GridSz+D3;
#ifdef PDEBUG
        fprintf(stderr, "CoCoDistributeTasks2DBlockCyclic: D1-mod part\n sk_ctr[%d,%d,%d] = %d, devidx = %d\n",
          D1,D2,D3, sk_ctr, devidx);
#endif
        autotune_controller->Tasks_per_unit_list[devidx][dev_sk_ctr_list[devidx]] = sk_ctr;
        dev_sk_ctr_list[devidx]++;
      }
      autotune_controller->Tasks_per_unit_num[devidx] += D3GridSz;
      if(devidx == autotune_controller->active_unit_num - 1) devidx = 0;
      else devidx++;
    }

  for (int D1 = D1GridSz_div; D1 < D1GridSz_div + D1GridSz_mod; D1++)
    for (int D2 = 0; D2 < D2GridSz_div; D2++){
      for (int D3 = 0; D3 < D3GridSz; D3++){
        sk_ctr = D1*D2GridSz*D3GridSz + D2*D3GridSz+D3;
#ifdef PDEBUG
        fprintf(stderr, "CoCoDistributeTasks2DBlockCyclic: D2-mod part\nsk_ctr[%d,%d,%d] = %d, devidx = %d\n",
          D1,D2,D3, sk_ctr, devidx);
#endif
        autotune_controller->Tasks_per_unit_list[devidx][dev_sk_ctr_list[devidx]] = sk_ctr;
        dev_sk_ctr_list[devidx]++;
      }
      autotune_controller->Tasks_per_unit_num[devidx] += D3GridSz;
      if(devidx == autotune_controller->active_unit_num - 1) devidx = 0;
      else devidx++;
    }

  for (int D1 = D1GridSz_div; D1 < D1GridSz_div + D1GridSz_mod; D1++)
    for (int D2 = D2GridSz_div; D2 < D2GridSz_div + D2GridSz_mod; D2++){
        for (int D3 = 0; D3 < D3GridSz; D3++){
          sk_ctr = D1*D2GridSz*D3GridSz + D2*D3GridSz+D3;
#ifdef PDEBUG
          fprintf(stderr, "CoCoDistributeTasks2DBlockCyclic: D1 & D2 mod part\nsk_ctr[%d,%d,%d] = %d, devidx = %d\n",
            D1,D2,D3, sk_ctr, devidx);
#endif
        autotune_controller->Tasks_per_unit_list[devidx][dev_sk_ctr_list[devidx]] = sk_ctr;
        dev_sk_ctr_list[devidx]++;
      }
      autotune_controller->Tasks_per_unit_num[devidx] += D3GridSz;
      if(devidx == autotune_controller->active_unit_num - 1) devidx = 0;
      else devidx++;
      }
    }
#ifdef PDEBUG
  fprintf(stderr, "CoCoDistributeTasks2DBlockCyclic:\nDistributing %ld Tasks to %d devices\n",
    autotune_controller->subkernel_num, autotune_controller->active_unit_num);
  fprintf(stderr, "Device Ids : [ ");
  for (int i =0; i < autotune_controller->active_unit_num; i++) fprintf(stderr, "%d ", autotune_controller->active_unit_id_list[i]);
  lprintf(0, "]\n");
  fprintf(stderr, "Subker Num : [ ");
  for (int i =0; i < autotune_controller->active_unit_num; i++) fprintf(stderr, "%d ",
    autotune_controller->Tasks_per_unit_num[i]);
  lprintf(0, "]\n");
  for (int i =0; i < autotune_controller->active_unit_num; i++){
    fprintf(stderr, "Subker Id list for dev_id = %d: [ ", autotune_controller->active_unit_id_list[i]);
    for (int j =0; j < autotune_controller->Tasks_per_unit_num[i]; j++) fprintf(stderr, "%d ",
      autotune_controller->Tasks_per_unit_list[i][j]);
    lprintf(0, "]\n");
  }
#endif
#ifdef DEBUG
	fprintf(stderr, "<-----|\n");
#endif

}
*/
/*****************************************************/
/// PARALia 2.0 - timed queues and blocks
/*long double Subkernel::run_op_estimate(MD_p modeler){
	run_op_est_t = modeler->getGPUexecFull(); 
#ifdef PDEBUG
	fprintf(stderr, "|-----> Subkernel(dev=%d,id=%d):run_op_estimate() -> run_op_est_t = %lf\n", 
		run_dev_id, id, run_op_est_t);
#endif
	return run_op_est_t; 
}*/
/*
#ifdef SUBKERNEL_SELECT_FETCH_ETA_PLUS_MIN_PENDING
Subkernel* SubkernelSelect(int dev_id, Subkernel** Subkernel_list, long Subkernel_list_len){
#ifdef SERIAL_SUBKERNEL_SELECTION
	for (int sk_idx = 0; sk_idx < Subkernel_list_len; sk_idx++)
		if(!Subkernel_list[sk_idx]->launched){
			Subkernel_list[sk_idx]->prepare_launch(dev_id);
			return Subkernel_list[sk_idx];
		}
	error("SubkernelSelect(SERIAL)\n No sk matched search condition\n");
#endif
	Subkernel* curr_sk = NULL;
	int sk_idx, potential_sks[Subkernel_list_len], tie_list_num = 0, doubletie_list_num = 0; 
	long double min_ETA = DBL_MAX;
	if(!Subkernel_list_len) error("SubkernelSelect: Gave 0 subkernel len with list = %p\n", Subkernel_list);
	for (sk_idx = 0; sk_idx < Subkernel_list_len; sk_idx++)if(!Subkernel_list[sk_idx]->launched){
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
			potential_sks[0] = sk_idx;
			tie_list_num = 1; 
		}
		else if(tmp_ETA == min_ETA){
		//else if(abs(tmp_ETA - min_ETA)/abs(tmp_ETA-csecond()) <= NORMALIZE_NEAR_SPLIT_LIMIT){
			potential_sks[tie_list_num++] = sk_idx;
		}
	}
	int most_fired_sks = -1, potential_tied_sks[tie_list_num];
	if (tie_list_num){
		potential_tied_sks[0] = potential_sks[0];
		doubletie_list_num = 1; 
	}
	else error("SubkernelSelect\n No sk matched search condition\n");
	for (int ctr = 0; ctr < tie_list_num; ctr++){
		curr_sk = Subkernel_list[potential_sks[ctr]];
		int tmp_fired_sks = 0; 
		for (int j = 0; j < curr_sk->TileNum; j++){
			if ( WR_LAZY == curr_sk->TileList[j]->WRP || WR == curr_sk->TileList[j]->WRP
				|| W_REDUCE == curr_sk->TileList[j]->WRP || WONLY == curr_sk->TileList[j]->WRP){
				tmp_fired_sks = curr_sk->TileList[j]->fired_times; 
			}
		}
		if(tmp_fired_sks > most_fired_sks){
			most_fired_sks = tmp_fired_sks;
			potential_tied_sks[0] = potential_sks[ctr];
			doubletie_list_num = 1; 
		}
		//else if(tmp_fired_sks == most_fired_sks){
		//else if(abs(tmp_ETA - min_ETA)/abs(tmp_ETA-csecond()) <= NORMALIZE_NEAR_SPLIT_LIMIT){
		//	potential_tied_sks[doubletie_list_num++] = potential_sks[ctr];
		//}
	}
	int selected_sk_idx = (doubletie_list_num)? 
		potential_tied_sks[int(rand() % doubletie_list_num)] : doubletie_list_num; 
	Subkernel_list[selected_sk_idx]->prepare_launch(dev_id);
	return Subkernel_list[selected_sk_idx];
}
#endif

#ifdef SUBKERNEL_SELECT_MIN_RONLY_ETA
Subkernel* SubkernelSelect(int dev_id, Subkernel** Subkernel_list, long Subkernel_list_len){
#ifdef SERIAL_SUBKERNEL_SELECTION
	Subkernel_list[0]->prepare_launch(dev_id);
	return Subkernel_list[0];
#endif
	Subkernel* curr_sk = NULL;
	int sk_idx;
	int potential_sks[Subkernel_list_len], tie_list_num = 0; 
	int max_fetches = -1;
	if(!Subkernel_list_len) error("SubkernelSelect: Gave 0 subkernel len with list = %p\n", Subkernel_list);
	for (sk_idx = 0; sk_idx < Subkernel_list_len; sk_idx++){
		curr_sk = Subkernel_list[sk_idx];
		int tmp_fetches = 0; 
		for (int j = 0; j < curr_sk->TileNum; j++){
			if(RONLY == curr_sk->TileList[j]->WRP)
				if(curr_sk->TileList[j]->loc_map[(dev_id)] == 0 || 
					curr_sk->TileList[j]->loc_map[(dev_id)] == 42) tmp_fetches++;
		}
		if(tmp_fetches > max_fetches){
			max_fetches = tmp_fetches;
			potential_sks[0] = sk_idx;
			tie_list_num = 1; 
		}
		else if(tmp_fetches == max_fetches){
			potential_sks[tie_list_num++] = sk_idx;
		}
	}
	int selected_sk_idx = (tie_list_num)? potential_sks[int(rand() % tie_list_num)] : tie_list_num; 
	swap_sk(&(Subkernel_list[0]), &(Subkernel_list[selected_sk_idx])); 
	Subkernel_list[0]->prepare_launch(dev_id);
	return Subkernel_list[0];
}
#endif

/*****************************************************/