///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief The possible subkernel distributions to different execution units.
///

#include "Autotuner.hpp"

#include <cmath>

void DistributeCompTasksRoundRobin(ATC_p autotune_controller){
#ifdef DEBUG
  	fprintf(stderr, "|-----> DistributeCompTasksRoundRobin(%p)\n", autotune_controller);
#endif
	error("PARALiA 3.0 does not support W-tile cimputation in multiple devices, use DistributeCompTasksRoundRobinChunk with Chunk = D3GridSz\n");
	if (autotune_controller->comp_task_num < autotune_controller->active_unit_num){
		int pred_active_unit_num = autotune_controller->active_unit_num;
		autotune_controller->active_unit_num = autotune_controller->comp_task_num;
		warning("DistributeCompTasksRoundRobin: Problem with predicted active_unit_num(%d) < comp_task_num(%ld) will be run with active_unit_num = %d\n",
			pred_active_unit_num, autotune_controller->comp_task_num, autotune_controller->active_unit_num);
		for (int d = 0 ; d < autotune_controller->active_unit_num; d++){
			autotune_controller->comp_task_per_unit_num[d] = 1;
			autotune_controller->comp_task_unit_list[d] = autotune_controller->active_unit_id_list[d];
		}
  	}
	else{
		int rem_dev = autotune_controller->comp_task_num;
		for (int d = 0 ; d < autotune_controller->active_unit_num; d++){
			autotune_controller->comp_task_per_unit_num[d] =
				(int) (1.0* autotune_controller->active_unit_score[d]* autotune_controller->comp_task_num);
			rem_dev-= autotune_controller->comp_task_per_unit_num[d];
		}
		while(rem_dev!= 0){
			for (int d = 0 ; d < autotune_controller->active_unit_num; d++){
				if(rem_dev!= 0){
					autotune_controller->comp_task_per_unit_num[d] += 1;
					rem_dev--;
				}
				else break;
			}
		}
	int total_sk_ctr = 0;
	short dev_sk_ctr_list[autotune_controller->active_unit_num];
	for(int devidx = 0; devidx < autotune_controller->active_unit_num; devidx++) dev_sk_ctr_list[devidx] = 0;
	while(total_sk_ctr<autotune_controller->comp_task_num){
		for(int devidx = 0; devidx < autotune_controller->active_unit_num; devidx++){
			if(total_sk_ctr == autotune_controller->comp_task_num) break;
			else if(dev_sk_ctr_list[devidx] == autotune_controller->comp_task_per_unit_num[devidx]) continue;
			else{
				autotune_controller->comp_task_unit_list[total_sk_ctr] = devidx;
				dev_sk_ctr_list[devidx]++;
				total_sk_ctr++;
			}
		}
	}
  }
#ifdef PDEBUG
	fprintf(stderr, "DistributeCompTasksRoundRobin:\nDistributing %ld Tasks to %d devices\n",
		autotune_controller->comp_task_num, autotune_controller->active_unit_num);
	fprintf(stderr, "Device Ids : [ ");
	for (int i =0; i < autotune_controller->active_unit_num; i++) fprintf(stderr, "%d ", 
		autotune_controller->active_unit_id_list[i]);
	fprintf(stderr, "]\n");
	fprintf(stderr, "Subker Num : [ ");
	for (int i =0; i < autotune_controller->active_unit_num; i++) fprintf(stderr, "%ld ", 
		autotune_controller->comp_task_per_unit_num[i]);
	fprintf(stderr, "]\n");
	fprintf(stderr, "Subker Id list: [ ");
	for (long int i =0; i < autotune_controller->comp_task_num; i++)
 		fprintf(stderr, "%d ", autotune_controller->comp_task_unit_list[i]);
	fprintf(stderr, "]\n");
#endif
#ifdef DEBUG
	fprintf(stderr, "<-----|\n");
#endif
}


void DistributeCompTasksRoundRobinChunk(ATC_p autotune_controller,  int Chunk_size){
#ifdef DEBUG
  	fprintf(stderr, "|-----> DistributeCompTasksRoundRobinChunk(%p, %d)\n", autotune_controller, Chunk_size);
#endif
#ifdef PDEBUG
	fprintf(stderr, "DistributeCompTasksRoundRobinChunk(%d): Devices = %d (scores = %s), task_num = %ld, task_buckets = %ld\n",
		Chunk_size, autotune_controller->active_unit_num, 
		printlist<double>(autotune_controller->active_unit_score, autotune_controller->active_unit_num),
		autotune_controller->comp_task_num, autotune_controller->comp_task_num/Chunk_size);
#endif
	autotune_controller->disable_caching = 1;
	autotune_controller->D1_parts = autotune_controller->D2_parts = 1; 
	if (autotune_controller->comp_task_num/Chunk_size + autotune_controller->comp_task_num%Chunk_size/1 <= autotune_controller->active_unit_num){
		int pred_active_unit_num = autotune_controller->active_unit_num;
		autotune_controller->active_unit_num = autotune_controller->comp_task_num/Chunk_size;
		warning("DistributeCompTasksRoundRobinChunk: Problem with predicted active_unit_num(%d) <= comp_task_num(=%ld)/Chunk_size"
		"will be run with active_unit_num = %d\n", pred_active_unit_num, autotune_controller->comp_task_num, autotune_controller->active_unit_num);
	for (int d = 0 ; d < autotune_controller->active_unit_num; d++){
		autotune_controller->comp_task_per_unit_num[d] = Chunk_size;
		for (int idx3 = 0; idx3 < Chunk_size; idx3++)
			autotune_controller->comp_task_unit_list[d*Chunk_size + idx3] = autotune_controller->active_unit_id_list[d];
	}
	}
  	else{
 		for (int d = 0 ; d < autotune_controller->active_unit_num; d++){
			autotune_controller->comp_task_per_unit_num[d] = Chunk_size*
			(int) (1.0* autotune_controller->active_unit_score[d]* (autotune_controller->comp_task_num/Chunk_size));
  		}
  		int sks_accounted_for = 0;
  		for (int d = 0 ; d < autotune_controller->active_unit_num; d++)
			sks_accounted_for += autotune_controller->comp_task_per_unit_num[d];
#ifdef PDEBUG
		fprintf(stderr, "Assigned kernel num to devices kernels (first pass): %s\n",
			printlist<long int>(autotune_controller->comp_task_per_unit_num, autotune_controller->active_unit_num));
#endif
  		int task_ctr = 0, dev_sk_ctr_list[autotune_controller->active_unit_num] = {0}, devidx = 0;
		for (int D1 = 0; D1 < autotune_controller->comp_task_num/Chunk_size; D1++){
			for (int D3 = 0; D3 < Chunk_size; D3++){
				int full_circle = autotune_controller->active_unit_num;
				while(dev_sk_ctr_list[devidx] == autotune_controller->comp_task_per_unit_num[devidx] && task_ctr < sks_accounted_for){ 
					if(!full_circle) error("DistributeCompTasks2DBlockCyclic: would enter infinite loop due to wrong"
					"comp_task_per_unit_num, terminating\n");
					if(devidx == autotune_controller->active_unit_num - 1) devidx = 0;
					else devidx++;
					full_circle--;
				}
				autotune_controller->comp_task_unit_list[task_ctr] = autotune_controller->active_unit_id_list[devidx];
				if(task_ctr >= sks_accounted_for) autotune_controller->comp_task_per_unit_num[devidx] ++;
				dev_sk_ctr_list[devidx]++;
#ifdef PDEBUG
			fprintf(stderr, "DistributeCompTasksRoundRobinChunk: task_ctr[%d,%d] = %d, devidx = %d\n",
				D1, D3, task_ctr, devidx);
#endif
			task_ctr++;
			}
			if(devidx == autotune_controller->active_unit_num - 1) devidx = 0;
			else devidx++;
		}
	}
	int dev_task_ctr[autotune_controller->active_unit_num] = {0};
	for (long int cidx = 0 ; cidx < autotune_controller->comp_task_num; cidx++){
		int dev_tmp = autotune_controller->comp_task_unit_list[cidx], dev_tmp_idx = -1;
		for (int dev_idx = 0; dev_idx < autotune_controller->active_unit_num; dev_idx++)
			if(dev_tmp == autotune_controller->active_unit_id_list[dev_idx]){
				dev_tmp_idx = dev_idx;
				break;
			} 
		autotune_controller->comp_task_per_unit_list[dev_tmp_idx][dev_task_ctr[dev_tmp_idx]++] = cidx;
	}
#ifdef PDEBUG
	fprintf(stderr, "DistributeCompTasksRoundRobinChunk(Chunk_size=%d):\nDistributing %ld Tasks to %d devices\n",
		Chunk_size, autotune_controller->comp_task_num, autotune_controller->active_unit_num);
	fprintf(stderr, "Device Ids : [ ");
	for (int i =0; i < autotune_controller->active_unit_num; i++) fprintf(stderr, "%d ", 
		autotune_controller->active_unit_id_list[i]);
	fprintf(stderr, "]\n");
	fprintf(stderr, "Subker Num : [ ");
	for (int i =0; i < autotune_controller->active_unit_num; i++) fprintf(stderr, "%ld ", 
		autotune_controller->comp_task_per_unit_num[i]);
	fprintf(stderr, "]\n");
	fprintf(stderr, "Subker Id list: [ ");
	for (long int i =0; i < autotune_controller->comp_task_num; i++)
 		fprintf(stderr, "%d ", autotune_controller->comp_task_unit_list[i]);
	fprintf(stderr, "]\n");
#endif
#ifdef DEBUG
	fprintf(stderr, "<-----|\n");
#endif
}

void DistributeCompTasks2DBlockCyclic(ATC_p autotune_controller, int D1GridSz, int D2GridSz, int D3GridSz){
#ifdef DEBUG
  	fprintf(stderr, "|-----> DistributeCompTasks2DBlockCyclic(%p, %d, %d, %d)\n", autotune_controller, D1GridSz, D2GridSz, D3GridSz);
#endif
 	if ((D2GridSz == D3GridSz) &&  (D2GridSz == 1)){
		warning("DistributeCompTasks2DBlockCyclic: D2GridSz==D3GridSz==1 -> using DistributeCompTasksRoundRobinChunk\n");
		return DistributeCompTasksRoundRobinChunk(autotune_controller, D3GridSz);
  	}
	if ((D1GridSz*D2GridSz) < autotune_controller->active_unit_num){
		warning("DistributeCompTasks2DBlockCyclic: D1GridSz*D2GridSz(%d) < autotune_controller->active_unit_num(%d)"
		", using DistributeCompTasksRoundRobinChunk instead\n", 
	  		D1GridSz*D2GridSz, autotune_controller->active_unit_num);
		return DistributeCompTasksRoundRobinChunk(autotune_controller, D3GridSz);
  	}

	// 2D Block cyclic device decomposition
	autotune_controller->D1_parts = std::sqrt(autotune_controller->active_unit_num);
	autotune_controller->D2_parts = autotune_controller->D1_parts;
	if (autotune_controller->D1_parts ==0) { autotune_controller->D2_parts = autotune_controller->active_unit_num; autotune_controller->D1_parts = 1; }
	else {
		// find the most square decomposition of autotune_controller->active_unit_num in autotune_controller->D1_parts x autotune_controller->D2_parts
		int g;
		for (g = autotune_controller->D1_parts+1; g>0; --g)
		if (autotune_controller->active_unit_num % g == 0) break;
		if (g==0) { autotune_controller->D1_parts = autotune_controller->active_unit_num; autotune_controller->D2_parts = 1; }
		//if (g==0) { autotune_controller->D1_parts = 1; autotune_controller->D2_parts = autotune_controller->active_unit_num; }
		else { autotune_controller->D1_parts = g; autotune_controller->D2_parts = autotune_controller->active_unit_num/g; }
	}
	//if(!autotune_controller->disable_caching){
	/// If ORDER_2DBC="D1_lesseq_D2", reverse layout. 
	if (!strcmp(ORDER_2DBC, "D1_lesseq_D2")){
		int tmp = autotune_controller->D1_parts;
		autotune_controller->D1_parts = autotune_controller->D2_parts;
		autotune_controller->D2_parts = tmp;
	}
	if(D1GridSz < autotune_controller->D1_parts || D2GridSz < autotune_controller->D2_parts){
		warning("DistributeCompTasks2DBlockCyclic:\nGrid(%d,%d) smaller than {D1,D2}_parts = (%d,%d)\
			using DistributeCompTasksRoundRobinChunk instead\n", D1GridSz, D2GridSz, autotune_controller->D1_parts, autotune_controller->D2_parts);
		return DistributeCompTasksRoundRobinChunk(autotune_controller, D3GridSz);
	}

	int D1GridSz_div = D1GridSz/autotune_controller->D1_parts, 
		D2GridSz_div = D2GridSz/autotune_controller->D2_parts,
		D1GridSz_mod = D1GridSz%autotune_controller->D1_parts, 
		D2GridSz_mod = D2GridSz%autotune_controller->D2_parts;
#ifdef PDEBUG
	fprintf(stderr, "DistributeCompTasks2DBlockCyclic(%d, %d, %d): Devices = %d (scores = %s), autotune_controller->D1_parts = %d, autotune_controller->D2_parts = %d\n",
		D1GridSz, D2GridSz, D3GridSz, autotune_controller->active_unit_num, 
		printlist<double>(autotune_controller->active_unit_score, autotune_controller->active_unit_num),autotune_controller->D1_parts, autotune_controller->D2_parts);
#endif

	/// Actual 2D block-cyclic distribution.
	int D1GridSz_mod_rem[64], D2GridSz_mod_rem[64];
	for(int idx = 0; idx < 64; idx++){
		D1GridSz_mod_rem[idx] = D1GridSz_mod;
		D2GridSz_mod_rem[idx] = D2GridSz_mod;
		for(int idy = 0; idy < 64; idy++)
			autotune_controller->C_Decom_grid[idx][idy][0] = autotune_controller->C_Decom_grid[idx][idy][1] = -42;
	}
	for(int idx = 0; idx < autotune_controller->D1_parts; idx++) 
	for(int idy = 0; idy < autotune_controller->D2_parts; idy++){
			if(!idx) autotune_controller->C_Decom_grid[idx][idy][0] = D1GridSz_div;
			else  autotune_controller->C_Decom_grid[idx][idy][0] = autotune_controller->C_Decom_grid[idx-1][idy][0] + D1GridSz_div;
			if(D1GridSz_mod_rem[idy]){
				autotune_controller->C_Decom_grid[idx][idy][0]++;
				D1GridSz_mod_rem[idy]--;
			}
			if (!idy) autotune_controller->C_Decom_grid[idx][idy][1] = D2GridSz_div;
			else  autotune_controller->C_Decom_grid[idx][idy][1] = autotune_controller->C_Decom_grid[idx][idy-1][1] + D2GridSz_div;
			if(D2GridSz_mod_rem[idx]){
				autotune_controller->C_Decom_grid[idx][idy][1]++;
				D2GridSz_mod_rem[idx]--;
			}
#ifdef PDEBUG
		fprintf(stderr, "DistributeCompTasks2DBlockCyclic: autotune_controller->C_Decom_grid[%d,%d] = (%d, %d)\n",
			idx, idy, autotune_controller->C_Decom_grid[idx][idy][0], autotune_controller->C_Decom_grid[idx][idy][1]);
#endif
	}
#ifndef PRODUCTION
	for(int idx = 0; idx < autotune_controller->D1_parts; idx++) 
		if(D2GridSz_mod_rem[idx]) error("DistributeCompTasks2DBlockCyclic: Remainder dimensions"
		" D2GridSz_mod_rem[%d] = %d not empty\n", idx, D2GridSz_mod_rem[idx]);
	
	for(int idy = 0; idy < autotune_controller->D2_parts; idy++)
		if(D1GridSz_mod_rem[idy]) error("DistributeCompTasks2DBlockCyclic: Remainder dimensions"
		" D1GridSz_mod_rem[%d] = %d not empty\n", idy, D1GridSz_mod_rem[idy]);
#endif

	int D1GridIdx = -1, D2GridIdx = -1, D3GridIdx = -1; 
	for(long int task_ctr = 0; task_ctr < autotune_controller->comp_task_num; task_ctr++){
		int comp_task_idx = task_ctr/D3GridSz;
		D1GridIdx = comp_task_idx/D2GridSz;
		D2GridIdx = comp_task_idx%D2GridSz;
		D3GridIdx = task_ctr%(D3GridSz);
		int devidx = -1;
		for(int idx = 0; idx < autotune_controller->D1_parts; idx++) 
			if(devidx == -1) for(int idy = 0; idy < autotune_controller->D2_parts; idy++)
				if(autotune_controller->C_Decom_grid[idx][idy][0] > D1GridIdx && autotune_controller->C_Decom_grid[idx][idy][1] > D2GridIdx){
					devidx = idx*autotune_controller->D2_parts + idy;
					break;
				}
#ifdef PDEBUG
		fprintf(stderr, "DistributeCompTasks2DBlockCyclic: task_ctr[%d,%d,%d] = %ld, devidx = %d\n",
			D1GridIdx,D2GridIdx,D3GridIdx, task_ctr, devidx);
#endif
		autotune_controller->comp_task_unit_list[task_ctr] = devidx;
		autotune_controller->comp_task_per_unit_list[devidx][autotune_controller->comp_task_per_unit_num[devidx]] = task_ctr;
		autotune_controller->comp_task_per_unit_num[devidx]++;
	}
#ifdef PDEBUG
	fprintf(stderr, "DistributeCompTasks2DBlockCyclic:\nDistributing %ld Tasks to %d devices\n",
		autotune_controller->comp_task_num, autotune_controller->active_unit_num);
	fprintf(stderr, "Device Ids : [ ");
	for (int i =0; i < autotune_controller->active_unit_num; i++) fprintf(stderr, "%d ", 
		autotune_controller->active_unit_id_list[i]);
	fprintf(stderr, "]\n");
	fprintf(stderr, "Subker Num : [ ");
	for (int i =0; i < autotune_controller->active_unit_num; i++) fprintf(stderr, "%ld ", 
		autotune_controller->comp_task_per_unit_num[i]);
	fprintf(stderr, "]\n");
	fprintf(stderr, "Subker Id list: [ ");
	for (long int i =0; i < autotune_controller->comp_task_num; i++)
 		fprintf(stderr, "%d ", autotune_controller->comp_task_unit_list[i]);
	fprintf(stderr, "]\n");
#endif
#ifdef DEBUG
	fprintf(stderr, "<-----|\n");
#endif

}
/* 
/// Previous version, created for D1_parts >= D2_parts and row-wise split priority. Div/mods don't work well for reverse layout.
void DistributeCompTasks2DBlockCyclic(ATC_p autotune_controller, int D1GridSz, int D2GridSz, int D3GridSz){
#ifdef DEBUG
  	fprintf(stderr, "|-----> DistributeCompTasks2DBlockCyclic(%p, %d, %d, %d)\n", autotune_controller, D1GridSz, D2GridSz, D3GridSz);
#endif
 	if ((D2GridSz == D3GridSz) &&  (D2GridSz == 1)){
		warning("DistributeCompTasks2DBlockCyclic: D2GridSz==D3GridSz==1 -> using DistributeCompTasksRoundRobinChunk\n");
		return DistributeCompTasksRoundRobinChunk(autotune_controller, D3GridSz);
  	}

	// 2D Block cyclic
	autotune_controller->D1_parts = std::sqrt(autotune_controller->active_unit_num);
	autotune_controller->D2_parts = autotune_controller->D1_parts;
	if (autotune_controller->D1_parts ==0) { autotune_controller->D2_parts = autotune_controller->active_unit_num; autotune_controller->D1_parts = 1; }
	else {
		// find the most square decomposition of autotune_controller->active_unit_num in autotune_controller->D1_parts x autotune_controller->D2_parts
		int g;
		for (g = autotune_controller->D1_parts+1; g>0; --g)
		if (autotune_controller->active_unit_num % g == 0) break;
		if (g==0) { autotune_controller->D1_parts = autotune_controller->active_unit_num; autotune_controller->D2_parts = 1; }
		//if (g==0) { autotune_controller->D1_parts = 1; autotune_controller->D2_parts = autotune_controller->active_unit_num; }
		else { autotune_controller->D1_parts = g; autotune_controller->D2_parts = autotune_controller->active_unit_num/g; }
	}
	//TODO: reverse layout
	//int tmp = autotune_controller->D1_parts;
	//autotune_controller->D1_parts = autotune_controller->D2_parts;
	//autotune_controller->D2_parts = tmp;
	if(D1GridSz < autotune_controller->D1_parts || D2GridSz < autotune_controller->D2_parts){
		warning("DistributeCompTasks2DBlockCyclic:\nGrid(%d,%d) smaller than {D1,D2}_parts = (%d,%d)\
			using DistributeCompTasksRoundRobinChunk instead\n", D1GridSz, D2GridSz, autotune_controller->D1_parts, autotune_controller->D2_parts);
			DistributeCompTasksRoundRobinChunk(autotune_controller, D3GridSz);
		return;
	}
	int D1GridSz_div = D1GridSz/autotune_controller->D1_parts*autotune_controller->D1_parts, D2GridSz_div = D2GridSz/autotune_controller->D2_parts*autotune_controller->D2_parts,
		D1GridSz_mod = D1GridSz%autotune_controller->D1_parts, D2GridSz_mod = D2GridSz%autotune_controller->D2_parts;
#ifdef PDEBUG
	fprintf(stderr, "DistributeCompTasks2DBlockCyclic(%d, %d, %d): Devices = %d (scores = %s), autotune_controller->D1_parts = %d, autotune_controller->D2_parts = %d\n",
		D1GridSz, D2GridSz, D3GridSz, autotune_controller->active_unit_num, 
		printlist<double>(autotune_controller->active_unit_score, autotune_controller->active_unit_num),autotune_controller->D1_parts, autotune_controller->D2_parts);
#endif

	if ((D1GridSz*D2GridSz) < autotune_controller->active_unit_num){
		warning("DistributeCompTasks2DBlockCyclic: D1GridSz*D2GridSz(%d) < autotune_controller->active_unit_num(%d)\n", 
	  		D1GridSz*D2GridSz, autotune_controller->active_unit_num);
		int pred_active_unit_num = D1GridSz*D2GridSz;
		autotune_controller->active_unit_num = autotune_controller->comp_task_num;
		warning("DistributeCompTasks2DBlockCyclic: Problem with predicted active_unit_num(%d) < comp_task_num(%ld) will be run with active_unit_num = %d\n",
			pred_active_unit_num, autotune_controller->comp_task_num, autotune_controller->active_unit_num);
		for (int d = 0 ; d < autotune_controller->active_unit_num; d++){
	  		autotune_controller->comp_task_per_unit_num[d] = D3GridSz;
	  		for (int idx3 = 0; idx3 < D3GridSz; idx3++)
				autotune_controller->comp_task_unit_list[d*D3GridSz + idx3] = autotune_controller->active_unit_id_list[d];
		}
  	}
  	else{
		long sks_accounted_for = 0;
		for (int d = 0 ; d < autotune_controller->active_unit_num; d++){
			autotune_controller->comp_task_per_unit_num[d] = D3GridSz * (
				(int) (autotune_controller->active_unit_score[d]* D1GridSz_div*D2GridSz_div));
	  		/// TODO: this is a fix because int (some_double) does not work for all doubles as intended 
	  		/// Will disrupt non-homogeneous splits!
			sks_accounted_for+= autotune_controller->comp_task_per_unit_num[d]; 
  		}
		if(!D1GridSz_mod && !D2GridSz_mod && sks_accounted_for < autotune_controller->comp_task_num){
			warning("DistributeCompTasks2DBlockCyclic: Questionable remainder from first pass %ld / %ld sub-kernels\n",
				autotune_controller->comp_task_num - sks_accounted_for, autotune_controller->comp_task_num);
			int buckets =  D1GridSz_div*D2GridSz_div; 
			int buckets_rem = (autotune_controller->comp_task_num - sks_accounted_for)/D3GridSz;
			int buckets_intended = buckets/autotune_controller->active_unit_num; 
			for (int d = 0 ; d < autotune_controller->active_unit_num; d++) 
			if(autotune_controller->comp_task_per_unit_num[d]/D3GridSz < buckets_intended && buckets_rem){
				autotune_controller->comp_task_per_unit_num[d]+= D3GridSz; 
				buckets_rem--;
			}
		}
#ifdef PDEBUG
		fprintf(stderr, "Assigned kernel num to devices kernels (first pass): %s\n",
			printlist<long int>(autotune_controller->comp_task_per_unit_num, autotune_controller->active_unit_num));
#endif
		int task_ctr, dev_sk_ctr_list[autotune_controller->active_unit_num] = {0}, devidx = 0;
		for (int D1 = 0; D1 < D1GridSz_div; D1++)
		for (int D2 = 0; D2 < D2GridSz_div; D2++)
		for (int D3 = 0; D3 < D3GridSz; D3++){
			task_ctr = D1*D2GridSz*D3GridSz + D2*D3GridSz+D3;
			devidx = D1/(D1GridSz/autotune_controller->D1_parts)*autotune_controller->D2_parts + D2/(D2GridSz/autotune_controller->D2_parts);
#ifdef PDEBUG
			fprintf(stderr, "DistributeCompTasks2DBlockCyclic: task_ctr[%d,%d,%d] = %d, devidx = %d\n",
				D1,D2,D3, task_ctr, devidx);
#endif
			int full_circle = autotune_controller->active_unit_num; 
			while(dev_sk_ctr_list[devidx] == autotune_controller->comp_task_per_unit_num[devidx]){
				if(!full_circle) error("DistributeCompTasks2DBlockCyclic: would enter infinite loop due to wrong comp_task_per_unit_num, terminating\n");
				if(devidx == autotune_controller->active_unit_num - 1) devidx = 0;
				else devidx++;
				full_circle--; 
			}
			autotune_controller->comp_task_unit_list[task_ctr] = autotune_controller->active_unit_id_list[devidx];
			dev_sk_ctr_list[devidx]++;
		}
		devidx = 0;
		for (int D1 = 0; D1 < D1GridSz_div; D1++)
		for (int D2 = D2GridSz_div; D2 < D2GridSz_div + D2GridSz_mod; D2++){
			for (int D3 = 0; D3 < D3GridSz; D3++){
				task_ctr = D1*D2GridSz*D3GridSz + D2*D3GridSz+D3;
#ifdef PDEBUG
				fprintf(stderr, "DistributeCompTasks2DBlockCyclic: D1-mod part\n task_ctr[%d,%d,%d] = %d, devidx = %d\n",
					D1,D2,D3, task_ctr, devidx);
#endif
				autotune_controller->comp_task_unit_list[task_ctr] = autotune_controller->active_unit_id_list[devidx];
				dev_sk_ctr_list[devidx]++;
			}
			autotune_controller->comp_task_per_unit_num[devidx] += D3GridSz;
			if(devidx == autotune_controller->active_unit_num - 1) devidx = 0;
			else devidx++;
		}
		for (int D1 = D1GridSz_div; D1 < D1GridSz_div + D1GridSz_mod; D1++)
		for (int D2 = 0; D2 < D2GridSz_div; D2++){
			for (int D3 = 0; D3 < D3GridSz; D3++){
				task_ctr = D1*D2GridSz*D3GridSz + D2*D3GridSz+D3;
#ifdef PDEBUG
				fprintf(stderr, "DistributeCompTasks2DBlockCyclic: D2-mod part\nsk_ctr[%d,%d,%d] = %d, devidx = %d\n",
					D1,D2,D3, task_ctr, devidx);
#endif
				autotune_controller->comp_task_unit_list[task_ctr] = autotune_controller->active_unit_id_list[devidx];
				dev_sk_ctr_list[devidx]++;
			}
			autotune_controller->comp_task_per_unit_num[devidx] += D3GridSz;
			if(devidx == autotune_controller->active_unit_num - 1) devidx = 0;
			else devidx++;
		}
		for (int D1 = D1GridSz_div; D1 < D1GridSz_div + D1GridSz_mod; D1++)
		for (int D2 = D2GridSz_div; D2 < D2GridSz_div + D2GridSz_mod; D2++){
			for (int D3 = 0; D3 < D3GridSz; D3++){
				task_ctr = D1*D2GridSz*D3GridSz + D2*D3GridSz+D3;
#ifdef PDEBUG
				fprintf(stderr, "DistributeCompTasks2DBlockCyclic: D1 & D2 mod part\nsk_ctr[%d,%d,%d] = %d, devidx = %d\n",
					D1,D2,D3, task_ctr, devidx);
#endif
				autotune_controller->comp_task_unit_list[task_ctr] = autotune_controller->active_unit_id_list[devidx];
				dev_sk_ctr_list[devidx]++;
			}
			autotune_controller->comp_task_per_unit_num[devidx] += D3GridSz;
			if(devidx == autotune_controller->active_unit_num - 1) devidx = 0;
			else devidx++;
		}
	}
#ifdef PDEBUG
	fprintf(stderr, "DistributeCompTasks2DBlockCyclic:\nDistributing %ld Tasks to %d devices\n",
		autotune_controller->comp_task_num, autotune_controller->active_unit_num);
	fprintf(stderr, "Device Ids : [ ");
	for (int i =0; i < autotune_controller->active_unit_num; i++) fprintf(stderr, "%d ", 
		autotune_controller->active_unit_id_list[i]);
	fprintf(stderr, "]\n");
	fprintf(stderr, "Subker Num : [ ");
	for (int i =0; i < autotune_controller->active_unit_num; i++) fprintf(stderr, "%ld ", 
		autotune_controller->comp_task_per_unit_num[i]);
	fprintf(stderr, "]\n");
	fprintf(stderr, "Subker Id list: [ ");
	for (long int i =0; i < autotune_controller->comp_task_num; i++)
 		fprintf(stderr, "%d ", autotune_controller->comp_task_unit_list[i]);
	fprintf(stderr, "]\n");
#endif
#ifdef DEBUG
	fprintf(stderr, "<-----|\n");
#endif

}*/

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
	for (int task_idx = 0; task_idx < Subkernel_list_len; task_idx++)
		if(!Subkernel_list[sk_idx]->launched){
			Subkernel_list[sk_idx]->prepare_launch(dev_id);
			return Subkernel_list[sk_idx];
		}
	error("SubkernelSelect(SERIAL)\n No sk matched search condition\n");
#endif
	Subkernel* curr_sk = NULL;
	int task_idx, potential_sks[Subkernel_list_len], tie_list_num = 0, doubletie_list_num = 0; 
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
			potential_sks[0] = task_idx;
			tie_list_num = 1; 
		}
		else if(tmp_ETA == min_ETA){
		//else if(abs(tmp_ETA - min_ETA)/abs(tmp_ETA-csecond()) <= NORMALIZE_NEAR_SPLIT_LIMIT){
			potential_sks[tie_list_num++] = task_idx;
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
	int task_idx;
	int potential_sks[Subkernel_list_len], tie_list_num = 0; 
	int max_fetches = -1;
	if(!Subkernel_list_len) error("SubkernelSelect: Gave 0 subkernel len with list = %p\n", Subkernel_list);
	for (sk_idx = 0; task_idx < Subkernel_list_len; task_idx++){
		curr_sk = Subkernel_list[sk_idx];
		int tmp_fetches = 0; 
		for (int j = 0; j < curr_sk->TileNum; j++){
			if(RONLY == curr_sk->TileList[j]->WRP)
				if(curr_sk->TileList[j]->loc_map[(dev_id)] == 0 || 
					curr_sk->TileList[j]->loc_map[(dev_id)] == 42) tmp_fetches++;
		}
		if(tmp_fetches > max_fetches){
			max_fetches = tmp_fetches;
			potential_sks[0] = task_idx;
			tie_list_num = 1; 
		}
		else if(tmp_fetches == max_fetches){
			potential_sks[tie_list_num++] = task_idx;
		}
	}
	int selected_sk_idx = (tie_list_num)? potential_sks[int(rand() % tie_list_num)] : tie_list_num; 
	swap_sk(&(Subkernel_list[0]), &(Subkernel_list[selected_sk_idx])); 
	Subkernel_list[0]->prepare_launch(dev_id);
	return Subkernel_list[0];
}
#endif

/*****************************************************/