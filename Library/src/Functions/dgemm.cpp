///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief The DGEMM implementation using the new mission-agent-asset C++ classes.
///

#include "backend_wrappers.hpp"
#include "chl_smart_wrappers.hpp"
#include "Autotuner.hpp"
#include "Decomposer.hpp"
#include "PARALiA.hpp"
#include "Resource_manager.hpp"
#include "DataCaching.hpp"

//#include <pthread.h>

ATC_p predef_controller_dgemm = NULL;

#ifdef STEST
double gemm_entry_ts;
#endif

void ManageCachesDgemm(PMD_p local_PMD, int aut_idx){
	for(int loc = 0; loc < CHL_MEMLOCS; loc++) if(local_PMD->autotuner[aut_idx]->Block_num[loc]){
#ifndef PRODUCTION
		if(local_PMD->autotuner[aut_idx]->Block_num[loc] == -42) 
			error("PARALiADgemm: local_PMD->autotuner[aut_idx]->Block_num[%d] is -42\n", loc);
#endif
		//long long prev_DevCache_sz = 0;
		//if (current_SAB[loc] != NULL) prev_DevCache_sz = (long long)
		//	current_SAB[loc]->BlockSize * current_SAB[loc]->BlockNum;
#ifdef BUFFER_REUSE_ENABLE
		if(current_SAB[loc] == NULL) current_SAB[loc] = 
			new Buffer(loc, local_PMD->autotuner[aut_idx]->Block_num[loc], local_PMD->autotuner[aut_idx]->Block_sz);
		else if (current_SAB[loc]->BlockSize != local_PMD->autotuner[aut_idx]->Block_sz 
			|| current_SAB[loc]->BlockNum < local_PMD->autotuner[aut_idx]->Block_num[loc]){
			error("PARALiADgemm:ManageCachesDgemm -> PARALiA 3.0 should not enter this\n");
#ifdef DEBUG
			fprintf(stderr, "PARALiADgemm: Previous Cache smaller than requested:\
			current_SAB[%d]->BlockSize=%lld vs local_PMD->autotuner[aut_idx]->Block_sz = %lld,\
			current_SAB[%d]->BlockNum=%d vs local_PMD->autotuner[aut_idx]->Block_num[loc] = %d\n",
			loc, current_SAB[loc]->BlockSize, local_PMD->autotuner[aut_idx]->Block_sz,
			loc, current_SAB[loc]->BlockNum, local_PMD->autotuner[aut_idx]->Block_num[loc]);
#endif
			delete current_SAB[loc];
			current_SAB[loc] = new Buffer(loc, local_PMD->autotuner[aut_idx]->Block_num[loc], local_PMD->autotuner[aut_idx]->Block_sz);
		}
#else
		if(current_SAB[loc]!= NULL) 
			error("PARALiADgemm: current_SAB[%d] was not NULL with reuse disabled\n", loc);
		current_SAB[loc] = new Buffer(loc, local_PMD->autotuner[aut_idx]->Block_num[loc], local_PMD->autotuner[aut_idx]->Block_sz);
#endif
	}
	for (int i = 0; i < CHL_MEMLOCS; i++) local_PMD->SAB[i] = current_SAB[i];
	return;
}

void CreateTasksDgemm(PMD_p local_PMD, int aut_idx){

#ifdef DEBUG
	fprintf(stderr, "|-----> CreateTasksDgemm(%p,%ld,%ld)\n",
		local_PMD, local_PMD->autotuner[aut_idx]->T, local_PMD->autotuner[aut_idx]->comp_task_num);
	fprintf(stderr,"MgridSz = %d, NgridSz = %d, KgridSz = %d\n",
		local_PMD->decom[0]->GridSz1, local_PMD->decom[1]->GridSz2, local_PMD->decom[0]->GridSz2);
	fprintf(stderr,"Mlast = %d, Nlast = %d, Klast = %d\n",
	local_PMD->decom[0]->Tile_map[local_PMD->decom[0]->GridSz1*local_PMD->decom[0]->GridSz2-1]->dim1,
	local_PMD->decom[1]->Tile_map[local_PMD->decom[1]->GridSz1*local_PMD->decom[1]->GridSz2-1]->dim2,
	local_PMD->decom[0]->Tile_map[local_PMD->decom[0]->GridSz1*local_PMD->decom[0]->GridSz2-1]->dim2);
#endif
	gemm_backend_in* initial_dgemm = (gemm_backend_in*) local_PMD->problem_wrap;
	//int current_ctr = 0;
	for (int mi = 0; mi < local_PMD->decom[0]->GridSz1; mi++){
		for (int ni = 0; ni < local_PMD->decom[1]->GridSz2; ni++){
			Tile2D_p C_tile = local_PMD->decom[2]->getTile(mi,ni);
			C_tile->W_op_num = local_PMD->decom[0]->GridSz2;
			C_tile->reduce_mult = initial_dgemm->beta; 
			C_tile->W_op_name = "MM_FP64";
			C_tile->W_op_params = (void**) malloc(C_tile->W_op_num*sizeof(void*));
			long int comp_task_idx = mi*local_PMD->decom[1]->GridSz2*local_PMD->decom[0]->GridSz2 + ni*local_PMD->decom[0]->GridSz2;
			C_tile->W_op_dev_id = local_PMD->autotuner[aut_idx]->comp_task_unit_list[comp_task_idx];
			C_tile->Block_reuses[C_tile->W_op_dev_id] = C_tile->W_op_num;
			if ((C_tile->WRP == WR_LAZY || C_tile->WRP == W_REDUCE) && C_tile->W_init_loc == C_tile->W_op_dev_id) C_tile->set_WRP(WR);
			C_tile->W_op_complete = new Event();
			C_tile->W_wb_complete = new Event();
			C_tile->W_ready = new Event();
			for (int ki = 0; ki < C_tile->W_op_num; ki++){
				Tile2D_p A_tile = local_PMD->decom[0]->getTile(mi,ki);
				if(A_tile->Block_reuses[C_tile->W_op_dev_id] == -42) A_tile->Block_reuses[C_tile->W_op_dev_id] = 1;
				else A_tile->Block_reuses[C_tile->W_op_dev_id]++;
				Tile2D_p B_tile = local_PMD->decom[1]->getTile(ki,ni);
				if(B_tile->Block_reuses[C_tile->W_op_dev_id] == -42) B_tile->Block_reuses[C_tile->W_op_dev_id] = 1;
				else B_tile->Block_reuses[C_tile->W_op_dev_id]++; 
				C_tile->W_op_params[ki] = (void*) malloc(sizeof(gemm_backend_in));
				gemm_backend_in*  ptr_ker_translate = 
					(gemm_backend_in*) C_tile->W_op_params[ki];
				ptr_ker_translate->TransA = initial_dgemm->TransA;
				ptr_ker_translate->TransB = initial_dgemm->TransB;
				ptr_ker_translate->M = C_tile->dim1;
				ptr_ker_translate->N = C_tile->dim2;
				if (ptr_ker_translate->TransA == 'N') ptr_ker_translate->K = A_tile->dim2;
				else if (ptr_ker_translate->TransA == 'T') ptr_ker_translate->K = A_tile->dim1;
				else error("CreateTasksDgemm: Unknown transpose type\n");
				ptr_ker_translate->A = NULL;
				ptr_ker_translate->B = NULL;
				ptr_ker_translate->C = NULL;
				ptr_ker_translate->alpha = initial_dgemm->alpha;
				if (!ki){
					if(WR_LAZY == C_tile->WRP || W_REDUCE == C_tile->WRP)
						set_val(C_tile->dtype, &ptr_ker_translate->beta, 0);
					else ptr_ker_translate->beta = initial_dgemm->beta;
				}
				else set_val(C_tile->dtype, &ptr_ker_translate->beta, 1);
				ptr_ker_translate->ldA = A_tile->ldim[C_tile->W_op_dev_id];
				ptr_ker_translate->ldB = B_tile->ldim[C_tile->W_op_dev_id];
				ptr_ker_translate->ldC = C_tile->ldim[C_tile->W_op_dev_id];
				//current_ctr = mi*local_PMD->decom[1]->GridSz2*local_PMD->decom[0]->GridSz2 
				//			+ ni*local_PMD->decom[0]->GridSz2 + ki;
				ptr_ker_translate->A_tile_v = (void*) A_tile; 
				ptr_ker_translate->B_tile_v = (void*) B_tile; 
				ptr_ker_translate->C_tile_v = (void*) C_tile; 

			}
		}
	}
#ifdef DEBUG
	fprintf(stderr, "<-----|\n");
#endif
}

void UpdateTasksDgemm(PMD_p local_PMD){
	gemm_backend_in* initial_dgemm = (gemm_backend_in*) local_PMD->problem_wrap;
	//int current_ctr = 0;
	for (int mi = 0; mi < local_PMD->decom[0]->GridSz1; mi++){
		for (int ni = 0; ni < local_PMD->decom[1]->GridSz2; ni++){
			Tile2D_p C_tile = local_PMD->decom[2]->getTile(mi,ni);
			C_tile->W_op_fired = 0;
			C_tile->reduce_mult = initial_dgemm->beta; 
			C_tile->Block_reuses[C_tile->W_op_dev_id] = C_tile->W_op_num;
			// TODO: These were done at the Tile2D reset already.
			//C_tile->W_op_complete->reset();
			//C_tile->W_wb_complete->reset();
			//C_tile->W_ready->reset();
			for (int ki = 0; ki < C_tile->W_op_num; ki++){
				Tile2D_p A_tile = local_PMD->decom[0]->getTile(mi,ki);
				if(A_tile->Block_reuses[C_tile->W_op_dev_id] == -42) A_tile->Block_reuses[C_tile->W_op_dev_id] = 1;
				else A_tile->Block_reuses[C_tile->W_op_dev_id]++;
				Tile2D_p B_tile = local_PMD->decom[1]->getTile(ki,ni);
				if(B_tile->Block_reuses[C_tile->W_op_dev_id] == -42) B_tile->Block_reuses[C_tile->W_op_dev_id] = 1;
				else B_tile->Block_reuses[C_tile->W_op_dev_id]++; 
				gemm_backend_in*  ptr_ker_translate = (gemm_backend_in*) C_tile->W_op_params[ki];
				ptr_ker_translate->A = NULL;
				ptr_ker_translate->B = NULL;
				ptr_ker_translate->C = NULL;
				ptr_ker_translate->alpha = initial_dgemm->alpha;
				if (!ki){
					if(WR_LAZY == C_tile->WRP || W_REDUCE == C_tile->WRP)
						set_val(C_tile->dtype, &ptr_ker_translate->beta, 0);
					else ptr_ker_translate->beta = initial_dgemm->beta;
				}
				else set_val(C_tile->dtype, &ptr_ker_translate->beta, 1.0);
				ptr_ker_translate->ldA = A_tile->ldim[C_tile->W_op_dev_id];
				ptr_ker_translate->ldB = B_tile->ldim[C_tile->W_op_dev_id];
				ptr_ker_translate->ldC = C_tile->ldim[C_tile->W_op_dev_id];
				//current_ctr = mi*local_PMD->decom[1]->GridSz2*local_PMD->decom[0]->GridSz2 
				//			+ ni*local_PMD->decom[0]->GridSz2 + ki;
			}
		}
	}
}

/// A dgemm wrapper including auto-tuning of T and cache_size, as well as device management
ATC_p PARALiADgemm(char TransA,  char TransB, long int M, long int N, long int K, double alpha, double* A, long int ldA,
		double* B, long int ldB, double beta, double* C, long int ldC)
{
	short lvl = 1;
#ifdef DEBUG
	fprintf(stderr, "|-----> PARALiADgemm(%c,%c,%zu,%zu,%zu,%lf,A=%p(%d),%zu,B=%p(%d),%zu,%lf,C=%p(%d),%zu)\n",
		TransA, TransB, M, N, K, alpha, A, CHLGetPtrLoc(A), ldA,
		B, CHLGetPtrLoc(B), ldB, beta, C, CHLGetPtrLoc(C), ldC);
#endif
#ifdef STEST
	gemm_entry_ts = csecond();
#endif
#ifdef TEST
	fprintf(stderr, "|-----> PARALiADgemm\n");
	double cpu_timer = csecond();
#endif
  	signal(SIGSEGV, handler);   // install segfault handler

	int prev_dev_id = CHLGetDevice();

	int reuse_problem_flag = 0;
	PMD_p local_PMD = NULL; 
	gemm_backend_in* initial_dgemm = NULL;

	for(int cache_entries = 0; cache_entries < PMD_cache_entries; cache_entries++)
	if(PMD_cache[cache_entries] && !strcmp(PMD_cache[cache_entries]->problem_name, "MM_FP64")){
			initial_dgemm = (gemm_backend_in*) PMD_cache[cache_entries]->problem_wrap;
#ifdef DEBUG
			PMD_cache[cache_entries]->print();
#endif
			reuse_problem_flag = 1; 
			if(initial_dgemm->TransA != TransA)
			reuse_problem_flag = 0;
			if(initial_dgemm->TransB != TransB)
				reuse_problem_flag = 0;
			if(PMD_cache[cache_entries]->autotuner[PMD_cache[cache_entries]->autotuner_best_idx]->M != M)
				reuse_problem_flag = 0;
			if(PMD_cache[cache_entries]->autotuner[PMD_cache[cache_entries]->autotuner_best_idx]->N != N)
				reuse_problem_flag = 0;
			if(PMD_cache[cache_entries]->autotuner[PMD_cache[cache_entries]->autotuner_best_idx]->K != K)
				reuse_problem_flag = 0;
			if(PMD_cache[cache_entries]->decom[0] && PMD_cache[cache_entries]->decom[0]->loc != CHLGetPtrLoc(A))
				reuse_problem_flag = 0;
			if(PMD_cache[cache_entries]->decom[1] && PMD_cache[cache_entries]->decom[1]->loc != CHLGetPtrLoc(B))
				reuse_problem_flag = 0;
			if(PMD_cache[cache_entries]->decom[2] && PMD_cache[cache_entries]->decom[2]->loc != CHLGetPtrLoc(C))
				reuse_problem_flag = 0;
			if(reuse_problem_flag){
				local_PMD = PMD_cache[cache_entries]; 
				break; 
			}
	}
	if (!local_PMD){
#ifdef DEBUG
			fprintf(stderr, "No previous problem metadata cache entry found, starting...\n");
#endif
			if (PMD_cache_entries == PROBLEM_MD_CACHE){
#ifndef PRODUCTION
				warning("PARALiADgemm - problem cache full, removing first entry\n");
#endif
				delete PMD_cache[0];
				local_PMD = PMD_cache[0] = new ProblemMetadata();
			}
			else{
				PMD_cache[PMD_cache_entries] = new ProblemMetadata();
				local_PMD = PMD_cache[PMD_cache_entries++];
			}
			local_PMD->problem_wrap = malloc(sizeof(gemm_backend_in));
	}
	else{;
#ifdef DEBUG
			fprintf(stderr, "Reusing local_PMD = %p with similar characteristics\n", local_PMD);
#endif
	}
	initial_dgemm = (gemm_backend_in*) local_PMD->problem_wrap;
	initial_dgemm->TransA = TransA;
	initial_dgemm->TransB = TransB;
	initial_dgemm->M = M;
	initial_dgemm->N = N;
	initial_dgemm->K = K;
	initial_dgemm->A = (void**) &A;
	initial_dgemm->B = (void**) &B;
	initial_dgemm->C = (void**) &C;
	initial_dgemm->alpha = &alpha;
	initial_dgemm->beta = &beta;
	initial_dgemm->ldA = ldA;
	initial_dgemm->ldB = ldB;
	initial_dgemm->ldC = ldC;
	initial_dgemm->dev_id = -1;

	double autotune_timer = 0;
	long int T; 
	int remaining_tasks = 0;

	int curr_autotuner_ctr = 0;
	if(!reuse_problem_flag){
		local_PMD->problem_name = "MM_FP64";
		local_PMD->decom_num = 3;
		local_PMD->decom[0] = new Decom2D( (void*) A, M, K, ldA, TransA, DOUBLE);
		local_PMD->decom[1] = new Decom2D( (void*) B, K, N, ldB, TransB, DOUBLE);
		local_PMD->decom[2] = new Decom2D( (void*) C, M, N, ldC, 'N', DOUBLE);

		local_PMD->decom[0]->prepareAsync();
		local_PMD->decom[1]->prepareAsync();
		local_PMD->decom[2]->prepareAsync();

		local_PMD->autotuner_ctr = 1;
		curr_autotuner_ctr = local_PMD->autotuner_best_idx = 0; 
		local_PMD->autotuner[curr_autotuner_ctr] = new ATC();
		if (predef_controller_dgemm && local_PMD->autotuner[curr_autotuner_ctr]->diff_intialized_params_ATC(predef_controller_dgemm))
			local_PMD->autotuner[curr_autotuner_ctr]->mimic_ATC(predef_controller_dgemm);
		autotune_timer = local_PMD->autotuner[curr_autotuner_ctr]->autotune_problem(local_PMD->problem_name, CHLGetPtrLoc(A), 
		CHLGetPtrLoc(B), CHLGetPtrLoc(C), 
		CHLGetPtrLoc(C), M, N, K, sizeof(double));
		if(predef_controller_dgemm && predef_controller_dgemm->T > 0)
			/// Cannot tune tile repetitively if its user-defined
			local_PMD->autotuner_ctr = REP_TILE*2 + 2;
		for(int idx = local_PMD->autotuner_ctr; idx < REP_TILE; idx++){
			local_PMD->autotuner[idx] = new ATC();
			if (predef_controller_dgemm && local_PMD->autotuner[idx]->diff_intialized_params_ATC(predef_controller_dgemm)) 
				local_PMD->autotuner[idx]->mimic_ATC(predef_controller_dgemm);
			local_PMD->autotuner[idx]->T = local_PMD->autotuner[curr_autotuner_ctr]->T_best_candidates[idx];
			local_PMD->autotuner[idx]->autotune_problem(local_PMD->problem_name, CHLGetPtrLoc(A), 
			CHLGetPtrLoc(B), CHLGetPtrLoc(C), 
			CHLGetPtrLoc(C), M, N, K, sizeof(double));
		}

		for(int d=0; d < CHL_MEMLOCS; d++) CHLEnableLinks(d, CHL_MEMLOCS);

		if (!strcmp(OUTPUT_ALGO_MODE, "ALGO_AUTO")) local_PMD->autotuner[curr_autotuner_ctr]->select_algo();
#ifdef TEST
		cpu_timer = csecond() - cpu_timer;
		fprintf(stderr, "Preparing decomposers -> t_prep = %lf ms\n", cpu_timer*1000);
		cpu_timer = csecond();
#endif
		for (int i = 0; i < CHL_MEMLOCS; i++) current_SAB[i] = NULL;
		ManageCachesDgemm(local_PMD, 0);
		T = local_PMD->autotuner[curr_autotuner_ctr]->T;
		local_PMD->decom[0]->InitTileMap(T, T, local_PMD->SAB, RONLY);
		local_PMD->decom[1]->InitTileMap(T, T, local_PMD->SAB, RONLY);
		WR_properties C_tile_prop; 
		if (!strcmp(OUTPUT_ALGO_MODE, "ALGO_WR")) C_tile_prop = WR;
		else if (!strcmp(OUTPUT_ALGO_MODE, "ALGO_WR_LAZY")) C_tile_prop = WR_LAZY;
		//else if (!strcmp(OUTPUT_ALGO_MODE, "ALGO_WREDUCE")) C_tile_prop = W_REDUCE;
		else error("PARALiADgemm: Unsupported OUTPUT_ALGO_MODE = %s\n", OUTPUT_ALGO_MODE);
		local_PMD->decom[2]->InitTileMap(T, T, local_PMD->SAB, C_tile_prop);
#ifdef TEST
		cpu_timer = csecond() - cpu_timer;
		fprintf(stderr, "Decomposing data to tiles -> t_tile = %lf ms\n", cpu_timer*1000);
		cpu_timer = csecond();
#endif
		CreateTasksDgemm(local_PMD, 0);
		remaining_tasks = local_PMD->autotuner[curr_autotuner_ctr]->task_num;
	}
	else if(REP_TILE!= 1 && local_PMD->autotuner_ctr <= REP_TILE*2 && local_PMD->autotuner_ctr%2 == 0){
		curr_autotuner_ctr = local_PMD->autotuner_ctr/2;
		if(curr_autotuner_ctr == REP_TILE) curr_autotuner_ctr =  local_PMD->autotuner_best_idx;
		local_PMD->decom[0]->MatrixReset((void*) A, ldA);
		local_PMD->decom[1]->MatrixReset((void*) B, ldB);
		local_PMD->decom[2]->MatrixReset((void*) C, ldC);
		for (int i = 0; i < CHL_MEMLOCS; i++) current_SAB[i] = NULL;
		ManageCachesDgemm(local_PMD, curr_autotuner_ctr);
		T = local_PMD->autotuner[curr_autotuner_ctr]->T;

		local_PMD->decom[0]->InitTileMap(T, T, local_PMD->SAB, RONLY);
		local_PMD->decom[1]->InitTileMap(T, T, local_PMD->SAB, RONLY);
		WR_properties C_tile_prop;
		if (!strcmp(OUTPUT_ALGO_MODE, "ALGO_WR")) C_tile_prop = WR;
		else if (!strcmp(OUTPUT_ALGO_MODE, "ALGO_WR_LAZY")) C_tile_prop = WR_LAZY;
		//else if (!strcmp(OUTPUT_ALGO_MODE, "ALGO_WREDUCE")) C_tile_prop = W_REDUCE;
		else error("PARALiADgemm: Unsupported OUTPUT_ALGO_MODE = %s\n", OUTPUT_ALGO_MODE);
		local_PMD->decom[2]->InitTileMap(T, T, local_PMD->SAB, C_tile_prop);
#ifdef TEST
		cpu_timer = csecond() - cpu_timer;
		fprintf(stderr, "Decomposing data to tiles -> t_tile = %lf ms\n", cpu_timer*1000);
		cpu_timer = csecond();
#endif
		CreateTasksDgemm(local_PMD, curr_autotuner_ctr);
		remaining_tasks = local_PMD->autotuner[curr_autotuner_ctr]->task_num;
		local_PMD->autotuner_ctr++;
	}
	else if(REP_TILE == 1 || local_PMD->autotuner_ctr > REP_TILE*2 || local_PMD->autotuner_ctr%2 == 1){
		if(REP_TILE != 1 && local_PMD->autotuner_ctr <= REP_TILE*2){
			curr_autotuner_ctr = local_PMD->autotuner_ctr/2;
			local_PMD->autotuner_ctr++;
		}
		else curr_autotuner_ctr = local_PMD->autotuner_best_idx;
		int buffer_freed = 0; 
		for (int i = 0; i < CHL_MEMLOCS; i++){
			current_SAB[i] = local_PMD->SAB[i];
			if(is_in_list (i, local_PMD->autotuner[curr_autotuner_ctr]->active_memlocs, 
			local_PMD->autotuner[curr_autotuner_ctr]->active_memloc_num)) buffer_freed = 1; 
		}
		if(buffer_freed) ManageCachesDgemm(local_PMD, curr_autotuner_ctr);
		T = local_PMD->autotuner[curr_autotuner_ctr]->T;
		local_PMD->decom[0]->Reset((void*) A, T, T, ldA, local_PMD->SAB);
		local_PMD->decom[1]->Reset((void*) B, T, T, ldB, local_PMD->SAB);
		local_PMD->decom[2]->Reset((void*) C, T, T, ldC, local_PMD->SAB);
#ifdef TEST
		cpu_timer = csecond() - cpu_timer;
		fprintf(stderr, "Re-assigning cache blocks to tiles -> t_tile = %lf ms\n", cpu_timer*1000);
		cpu_timer = csecond();
#endif
		UpdateTasksDgemm(local_PMD);
		remaining_tasks = local_PMD->autotuner[curr_autotuner_ctr]->task_num;
	}
	else error("PARALiADgemm: Should not reach here - local_PMD->autotuner_ctr = %d\n", local_PMD->autotuner_ctr);
	conserve_memory_curr = local_PMD->autotuner[curr_autotuner_ctr]->conserve_memory; 

#ifdef TEST
	cpu_timer = csecond() - cpu_timer;
	fprintf(stderr, "Updated subkernels for devices: t_update = %lf ms\n", cpu_timer*1000);
	cpu_timer = csecond();
#endif
	for(int d=0; d < CHL_MEMLOCS; d++)
		if(current_SAB[d]) current_SAB[d]->allocate();

#ifdef TEST
	cpu_timer = csecond() - cpu_timer;
	fprintf(stderr, "Memory management: t_mem = %lf ms\n", cpu_timer*1000);
	cpu_timer = csecond();
#endif

	RMInitResources(local_PMD->autotuner[curr_autotuner_ctr]->active_unit_id_list, local_PMD->autotuner[curr_autotuner_ctr]->active_unit_num);
	//RMInitWS(local_PMD->autotuner[curr_autotuner_ctr]->active_unit_id_list, local_PMD->autotuner[curr_autotuner_ctr]->active_unit_num);
	CHLSyncCheckErr();
#ifdef TEST
	cpu_timer = csecond() - cpu_timer;
	fprintf(stderr, "Queue/Handle init: t_resource = %lf ms\n", cpu_timer*1000);
	cpu_timer = csecond();
#endif
	double run_timer = csecond();
	while (remaining_tasks){
		Ttask_p curr = local_PMD->autotuner[curr_autotuner_ctr]->task_list[local_PMD->autotuner[curr_autotuner_ctr]->task_num - remaining_tasks];
		Tile2D_p target_tile = local_PMD->decom[curr->DecomIdx]->getTile(curr->TileIdx, curr->TileIdy);
		if (curr->type == FETCH) target_tile->fetch(curr->predef_route);
		else if (curr->type == COMPUTE)target_tile->run_operation(curr->op_id, curr->predef_route);
		else if (curr->type == WRITEBACK) target_tile->writeback(curr->predef_route);
		remaining_tasks--;
		//usleep(3000);
	}
#ifdef TEST
	cpu_timer = csecond() - cpu_timer;
	fprintf(stderr, "Subkernels launched: t_sk_fire = %lf ms\n", cpu_timer*1000);
	cpu_timer = csecond();
#endif
	//local_PMD->decom[2]->WBTileMap();
#ifdef TEST
	cpu_timer = csecond() - cpu_timer;
	fprintf(stderr, "Writebacks launched -> t_wb_fire = %lf ms\n", cpu_timer*1000);
	cpu_timer = csecond();
#endif
//#ifndef ENABLE_SEND_RECV_OVERLAP
//	sync_recv_queues();
//	local_PMD->decom[2]->WBTileMap();
//#endif
	local_PMD->decom[2]->SyncTileMap();
	CHLSyncCheckErr();
	run_timer = csecond() - run_timer;
	if(local_PMD->best_t > run_timer){
		local_PMD->autotuner_best_idx = curr_autotuner_ctr;
		local_PMD->best_t = run_timer;
	}
#ifdef TEST
	cpu_timer = run_timer;
	fprintf(stderr, "Synced result -> t_exec_full = %lf ms\n", cpu_timer*1000);
	fprintf(stderr, "t_predicted for T=%zu was %.2lf ms : %lf percentile error\n", T, local_PMD->autotuner[curr_autotuner_ctr]->pred_t*1000,
	(local_PMD->autotuner[curr_autotuner_ctr]->pred_t==0)? 0.0: (local_PMD->autotuner[curr_autotuner_ctr]->pred_t - cpu_timer )/local_PMD->autotuner[curr_autotuner_ctr]->pred_t*100);
	cpu_timer = csecond();
#endif
#ifdef PDEBUG
	fprintf(stderr, "PARALiADgemm(): completed for PMD_cache_entries = %d\n", PMD_cache_entries);
#endif
#ifdef STEST
	STEST_print_SK(thread_dev_data, gemm_entry_ts, local_PMD->autotuner[curr_autotuner_ctr]->active_unit_num);
#endif

#ifdef TTEST
	HopMemcpyPrint();
	//n_HopMemcpyPrint();
#endif

#ifdef DDEBUG
  local_PMD->decom[0]->DrawTileMap();
  local_PMD->decom[1]->DrawTileMap();
  local_PMD->decom[2]->DrawTileMap();
#endif

#ifdef CDEBUG
	for(int i = 0; i < CHL_MEMLOCS; i++) if(local_PMD->SAB[i]) local_PMD->SAB[i]->draw_buffer(true);
#endif

#ifdef BUFFER_REUSE_ENABLE
	for(int i = 0; i < CHL_MEMLOCS; i++) if(local_PMD->SAB[i]) local_PMD->SAB[i]->reset(true);
#else
	for(int i = 0 ; i < CHL_MEMLOCS; i++) if(local_PMD->SAB[i]){
		delete local_PMD->SAB[i];
		local_PMD->SAB[i] = current_SAB[i] = NULL;
	}
#endif

#ifdef QUEUE_REUSE_ENABLE
	RMCleanResources();
#else 
	RMFreeResources();
#endif

#ifdef TEST
	CHLSyncCheckErr();
	cpu_timer = csecond() - cpu_timer;
	fprintf(stderr, "Invalidate caches -> t_invalidate = %lf ms\n", cpu_timer*1000);
	cpu_timer = csecond();
#endif
#ifdef METADATA_REUSE_PROBLEMS
	;
#else
	local_PMD->decom[0]->DestroyTileMap();
	local_PMD->decom[1]->DestroyTileMap();
	local_PMD->decom[2]->DestroyTileMap();
#ifdef TEST
	cpu_timer = csecond() - cpu_timer;
	fprintf(stderr, "Destroyed Tilemaps -> t_invalidate = %lf ms\n", cpu_timer*1000);
	cpu_timer = csecond();
#endif

	CHLSelectDevice(prev_dev_id);
    local_PMD->decom[0]->resetProperties();
    local_PMD->decom[1]->resetProperties();
    local_PMD->decom[2]->resetProperties();
	delete local_PMD->decom[0];
	delete local_PMD->decom[1];
	delete local_PMD->decom[2];
#ifdef TEST
	cpu_timer = csecond() - cpu_timer;
	fprintf(stderr, "Unregistering decomposers -> t_unpin = %lf ms\n", cpu_timer*1000);
#endif
#endif

#ifdef DEBUG
	fprintf(stderr, "<-----|\n");
#endif
#ifdef TEST
	fprintf(stderr, "<-----|\n");
#endif

	predef_controller_dgemm = NULL;
	// Better not return our global to the user, he can accidentally do stuff to it.
	ATC_p result = new ATC();
	result->mimic_ATC(local_PMD->autotuner[curr_autotuner_ctr]);
	return result;
}

/// A modification of PARALiADgemm but with given parameters (mainly for performance/debug purposes)
ATC_p PARALiADgemmControled(char TransA,  char TransB, long int M, long int N, long int K, double alpha, double* A, long int ldA,
		double* B, long int ldB, double beta, double* C, long int ldC, ATC_p predef_controller){
	if (predef_controller == NULL){
		warning("Calling PARALiADgemmControled with empty controller -> "
		"falling back to full autotune version \'PARALiADgemm\'\n");
	}
	else predef_controller_dgemm = predef_controller;
	return PARALiADgemm(TransA, TransB,  M, N, K, alpha, A, ldA, B, ldB, beta, C, ldC);
}

/// A modification of PARALiADgemm but with given parameters (mainly for performance/debug purposes)
ATC_p PARALiADgemmLarge(char TransA,  char TransB, long int M, long int N, long int K, double alpha, double* A, long int ldA,
		double* B, long int ldB, double beta, double* C, long int ldC){
	long ex_tile = 40960;
	if(predef_controller_dgemm!=NULL){
		if(predef_controller_dgemm->active_unit_num == 1) ex_tile = 16384;
		else if(predef_controller_dgemm->active_unit_num == 2) ex_tile = 32768;
	}
	if(M > ex_tile && N > ex_tile && K > ex_tile){ // Run large decom version to fit in GPU mems
		warning("Running PARALiADgemmLarge version with ex_tile = %ld, still under development\n", ex_tile);
		ATC_p tile_force_atc = new ATC();
		if(predef_controller_dgemm) tile_force_atc->mimic_ATC(predef_controller_dgemm);
		tile_force_atc->T = ex_tile/8;
		warning("Setting tiling size to %ld\n", tile_force_atc->T);
		long Grid_M = M/ex_tile, Grid_N = N/ex_tile, Grid_K = K/ex_tile;
		int Grid_M_div = M%ex_tile, Grid_N_div = N%ex_tile, Grid_K_div = K%ex_tile;
		// TODO: Padding instead of resize so all data fit in buffer without complex mechanism.
		// Can degrade performance for small div sizes.
		if (Grid_M_div > 0) Grid_M++;
		else Grid_M_div = ex_tile;
		if (Grid_N_div > 0) Grid_N++;
		else Grid_N_div = ex_tile;
		if (Grid_K_div > 0) Grid_K++;
		else Grid_K_div = ex_tile;
		ATC_p temp;
		int curr_chunk_M, curr_chunk_N, curr_chunk_K;
		for (long int itt1 = 0; itt1 < Grid_M; itt1++){
    		if ( itt1 == Grid_M - 1) curr_chunk_M = Grid_M_div;
    		else  curr_chunk_M = ex_tile;
			for (long int itt2 = 0 ; itt2 < Grid_N; itt2++){
      			if ( itt2 == Grid_N - 1) curr_chunk_N = Grid_N_div;
				else  curr_chunk_N = ex_tile;
				for (long int itt3 = 0 ; itt3 < Grid_K; itt3++){
      				if ( itt3 == Grid_K - 1) curr_chunk_K = Grid_K_div;
					else  curr_chunk_K = ex_tile;
					double beta_ck;
					if(!itt3) beta_ck = beta;
					else beta_ck = 1.0; 
					double* addr_chunk_A, *addr_chunk_B, *addr_chunk_C;
					if (TransA == 'N') addr_chunk_A = ((double*)A) + ((long)(itt1*ex_tile + itt3*ex_tile*ldA)); 
					else addr_chunk_A = ((double*)A) + ((long)(itt1*ldA*ex_tile + itt3*ex_tile)); 
					if (TransB == 'N') addr_chunk_B = ((double*)B) + ((long)(itt3*ex_tile + itt2*ex_tile*ldB)); 
					else  addr_chunk_B = ((double*)B) + ((long)(itt3*ldB*ex_tile + itt2*ex_tile)); 
					addr_chunk_C = ((double*)C) + ((long)(itt1*ex_tile + itt2*ex_tile*ldC)); 
					temp = PARALiADgemmControled(TransA, TransB,  curr_chunk_M, curr_chunk_N, curr_chunk_K, alpha, 
						addr_chunk_A, ldA, addr_chunk_B, ldB, beta_ck, addr_chunk_C, ldC, tile_force_atc);
					CHLSyncCheckErr();
				}
			}
		}
		return temp;
	}
	else return PARALiADgemm(TransA, TransB,  M, N, K, alpha, A, ldA, B, ldB, beta, C, ldC);
}

/// A modification of PARALiADgemm but with given parameters (mainly for performance/debug purposes)
ATC_p PARALiADgemmLargeControled(char TransA,  char TransB, long int M, long int N, long int K, double alpha, double* A, long int ldA,
		double* B, long int ldB, double beta, double* C, long int ldC, ATC_p predef_controller){
	if (predef_controller == NULL){
		warning("Calling PARALiADgemmLargeControled with empty controller -> "
		"falling back to full autotune version \'PARALiADgemmLarge\'\n");
	}
	else predef_controller_dgemm = predef_controller;
	return PARALiADgemmLarge(TransA, TransB,  M, N, K, alpha, A, ldA, B, ldB, beta, C, ldC);
}
