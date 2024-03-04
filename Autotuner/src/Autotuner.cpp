///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief The Autotuner controller functions.
///

#include "Autotuner.hpp"
#include "TaskDistribution.hpp"

#include <cfloat> /// For DBL_MAX
#include <cmath>

TileTask::TileTask(TileTaskType type_in, int DecomIdx_in, int TileIdx_in, 
	int TileIdy_in, int op_id_in, LinkRoute_p predef_route_in){
	type = type_in;
	DecomIdx = DecomIdx_in;
	TileIdx = TileIdx_in; 
	TileIdy = TileIdy_in; 
	op_id = op_id_in;
	predef_route = predef_route_in;
	print();
}

void TileTask::print(){
	fprintf(stderr,  "| --------- |\n");
	if(type == FETCH) fprintf(stderr,  "| --FETCH-- |\n");
	else if(type == COMPUTE) fprintf(stderr,  "| -COMPUTE- |\n");
	else if(type == WRITEBACK) fprintf(stderr,  "| WRITEBACK |\n");
	fprintf(stderr,  "| %1d.[%2d,%2d] |\n", DecomIdx, TileIdx, TileIdy);
	fprintf(stderr,  "|  Oper=%2d  |\n", op_id);
	fprintf(stderr,  "| %9s |\n", printlist<int>(predef_route->hop_uid_list, predef_route->hop_num));
}

ATC::ATC(){
	short lvl = 2;
#ifdef DEBUG
		fprintf(stderr,  "|-----> ATC::ATC\n");
#endif
	active_memlocs = (int*) malloc (CHL_MEMLOCS*sizeof(int));
	active_unit_id_list = (int*) malloc(CHL_WORKERS*sizeof(int));
	active_unit_score = (double*) malloc(CHL_WORKERS*sizeof(double));
	comp_task_per_unit_num = (long int*) malloc(CHL_WORKERS*sizeof(long int));
	task_list = NULL;
	comp_task_unit_list = NULL;
	comp_task_per_unit_list = NULL;
	T = active_unit_num = task_num = comp_task_num = -1;
	pred_t = pred_J = power_delay = energy_delay = -1.0;
	T_aggregate_sl = T_remainder_sl = T_small_sl = T_sknum_sl = T_big_sl = 0.0;

	cache_limit = 0;

	inter_grid = NULL;
#ifdef DEBUG
	fprintf(stderr,  "<-----|\n");
#endif
}

ATC::~ATC(){
	short lvl = 2;
#ifdef DEBUG
		fprintf(stderr,  "|-----> ATC::~ATC\n");
#endif
	free(active_memlocs);
	free(active_unit_id_list);
	free(active_unit_score);
	free(comp_task_per_unit_num);
	free(comp_task_unit_list); 
	for (int idx = 0; idx < active_unit_num; idx++){
		free(comp_task_per_unit_list[idx]); 
	}
	free(comp_task_per_unit_list); 
	for (int idx = 0; idx < task_num; idx++){
		delete task_list[idx]; 
	}
	free(task_list);
	delete inter_grid; 
#ifdef DEBUG
	fprintf(stderr,  "<-----|\n");
#endif
}

void ATC::reset(){
	short lvl = 4;
#ifdef DEBUG
		fprintf(stderr,  "|-----> ATC::reset\n");
#endif
	free(comp_task_per_unit_num);
	free(comp_task_unit_list);
	for (int idx = 0; idx < active_unit_num; idx++){
		free(comp_task_per_unit_list[idx]); 
	}
	free(comp_task_per_unit_list); 
	comp_task_per_unit_list = NULL;
	for (int idx = 0; idx < task_num; idx++){
		delete task_list[idx]; 
	}
	free(task_list);
	task_list = NULL;
	delete inter_grid;
	T = active_unit_num = task_num = comp_task_num = -1;
	pred_t = -1.0;
	cache_limit = 0;
#ifdef DEBUG
	fprintf(stderr,  "<-----|\n");
#endif
}

int ATC::diff_intialized_params_ATC(ATC_p other_ATC){
	short lvl = 3;
	#ifdef DEBUG
		fprintf(stderr,  "|-----> ATC::diff_intialized_params_ATC(other_ATC = %p)\n", other_ATC);
	#endif
	int result = 0;
	if(other_ATC->T != -1 && other_ATC->T != T){
		result++;
#ifdef PDEBUG
		fprintf(stderr,  "ATC::diff_intialized_params_ATC(): other_ATC->T = %ld, T = %ld\n", other_ATC->T, T);
#endif
		}
	if(other_ATC->active_unit_num != -1 && other_ATC->active_unit_num != active_unit_num){
		result++;
#ifdef PDEBUG
		fprintf(stderr,  "ATC::diff_intialized_params_ATC(): other_ATC->active_unit_num = %d, active_unit_num = %d\n",
			other_ATC->active_unit_num, active_unit_num);
#endif
	}
	else if(other_ATC->active_unit_num != -1 && other_ATC->active_unit_num == active_unit_num){
		for (int ctr = 0; ctr < active_unit_num; ctr++) if(other_ATC->active_unit_id_list[ctr] != active_unit_id_list[ctr]){
			result++;
#ifdef PDEBUG
		fprintf(stderr,  "ATC::diff_intialized_params_ATC(): other_ATC->active_unit_id_list[%d] = %d, active_unit_id_list[%d] = %d\n",
			ctr, other_ATC->active_unit_id_list[ctr], ctr, active_unit_id_list[ctr]);
#endif
			break;
		}
	}
	if(other_ATC->cache_limit != 0 && other_ATC->cache_limit != cache_limit){
		result++;
#ifdef PDEBUG
		fprintf(stderr,  "ATC::diff_intialized_params_ATC(): other_ATC->cache_limit = %lld, cache_limit = %lld\n",
			other_ATC->cache_limit, cache_limit);
#endif
	}
#ifdef DEBUG
	fprintf(stderr,  "<-----|\n");
#endif
	return result;
}

void ATC::mimic_ATC(ATC_p other_ATC){
#ifdef DEBUG
	fprintf(stderr,  "|-----> ATC::mimic_ATC(other_ATC = %p)\n", other_ATC);
#endif
	// TODO: Is this needed?
	//reset();
	T = other_ATC->T;
	T_aggregate_sl = other_ATC->T_aggregate_sl;
	T_imbalance_sl = other_ATC->T_imbalance_sl;
	T_remainder_sl = other_ATC->T_remainder_sl;
	T_small_sl = other_ATC->T_small_sl;
	T_sknum_sl = other_ATC->T_sknum_sl;
	T_big_sl = other_ATC->T_big_sl;
	active_unit_num = other_ATC->active_unit_num;
	for (int d = 0; d < other_ATC->active_unit_num; d++) active_unit_id_list[d] = other_ATC->active_unit_id_list[d];
	for (int d = 0; d < other_ATC->active_unit_num; d++) active_unit_score[d] = other_ATC->active_unit_score[d];
	pred_t = other_ATC->pred_t;
	pred_J = other_ATC->pred_J;
	power_delay = other_ATC->power_delay;
	energy_delay = other_ATC->energy_delay;
	pred_t_pesimistic = other_ATC->pred_t_pesimistic;
	pred_J_pesimistic = other_ATC->pred_J_pesimistic;
	power_delay_pesimistic = other_ATC->power_delay_pesimistic;
	energy_delay_pesimistic = other_ATC->energy_delay_pesimistic;
	cache_limit = other_ATC->cache_limit;
	// The following is not implemented
	//inter_grid->copy(other_ATC->inter_grid);

	if (other_ATC->task_num != -1){
		task_num = other_ATC->task_num;
		update_comp_task_num(other_ATC->comp_task_num);
		for (long int taskidx = 0; taskidx < comp_task_num; taskidx++) 
			comp_task_unit_list[taskidx] = other_ATC->comp_task_unit_list[taskidx];
		for (int dev = 0; dev < CHL_WORKERS; dev++) 
			comp_task_per_unit_num[dev] = other_ATC->comp_task_per_unit_num[dev];
		for (int dev_id = 0; dev_id < active_unit_num; dev_id++)
			for (long int task_comp_idx = 0; task_comp_idx < comp_task_per_unit_num[dev_id]; task_comp_idx++)
				comp_task_per_unit_list[dev_id][task_comp_idx] = 
					other_ATC->comp_task_per_unit_list[dev_id][task_comp_idx];
			
		if(!task_list) task_list = (Ttask_p*) malloc(task_num*sizeof(Ttask_p));
		for (int sk = 0; sk < task_num; sk++) task_list[sk] = other_ATC->task_list[sk];
	}

#ifdef DEBUG
	fprintf(stderr,  "<-----|\n");
#endif
}

void ATC::update_comp_task_num(long int comp_task_num_in){
#ifdef DEBUG
	fprintf(stderr,  "|-----> ATC::update_comp_task_num\n");
#endif
	int prev_comp_task_num = comp_task_num;
	comp_task_num = comp_task_num_in;
	if (prev_comp_task_num == -1){
		comp_task_unit_list = (int*) malloc(comp_task_num*sizeof(int));
		comp_task_per_unit_list = (int**) malloc(active_unit_num*sizeof(int*));
		for(int dev_id = 0; dev_id < active_unit_num; dev_id++) 
			comp_task_per_unit_list[dev_id] = (int*) malloc(comp_task_num*sizeof(int));
	}
	else if (prev_comp_task_num < comp_task_num){
		warning("ATC::update_comp_task_num: updating predefined comp_task_unit_list is untested\n");
		free(comp_task_unit_list);
		for(int dev_id = 0; dev_id < active_unit_num; dev_id++) free(comp_task_per_unit_list[dev_id]);
		free(comp_task_per_unit_list);
		comp_task_unit_list = (int*) malloc(comp_task_num*sizeof(int));
		comp_task_per_unit_list = (int**) malloc(active_unit_num*sizeof(int*));
		for(int dev_id = 0; dev_id < active_unit_num; dev_id++)
			comp_task_per_unit_list[dev_id] = (int*) malloc(comp_task_num*sizeof(int));
	}
	else error("ATC::update_comp_task_num was not updated correctly\n");
#ifdef DEBUG
	fprintf(stderr,  "<-----|\n");
#endif
}

void ATC::distribute_comp_tasks(){
	if (!strcmp(DISTRIBUTION, "ROUND-ROBIN"))
		DistributeCompTasksRoundRobin(this);
	else if (!strcmp(DISTRIBUTION, "SPLIT-NAIVE"))
		DistributeCompTasksNaive(this);
	else if (!strcmp(DISTRIBUTION, "SPLIT-CHUNKS-ROBIN"))
		DistributeCompTasksRoundRobinChunk(this, Grid_K);
	else if (!strcmp(DISTRIBUTION, "SPLIT-CHUNKS-ROBIN-REVERSE"))
		DistributeCompTasksRoundRobinChunkReverse(this, Grid_K);
	else if (!strcmp(DISTRIBUTION, "2D-BLOCK-CYCLIC"))
		DistributeCompTasks2DBlockCyclic(this, Grid_M, Grid_N, Grid_K);
	else error("ATC::distribute_comp_tasks: Unknown Subkernel Distribution %s\n", DISTRIBUTION);
	int dev_task_ctr[active_unit_num] = {0};
	for (long int cidx = 0 ; cidx < comp_task_num; cidx++){
		int dev_tmp = comp_task_unit_list[cidx]; 
		comp_task_per_unit_list[dev_tmp][dev_task_ctr[dev_tmp]++] = cidx;
	}
#ifdef PDEBUG
    print();
#endif
}

void ATC::initialize_tasks(){
	if(comp_task_num == -1) error("ATC::initialize_tasks: Called with undefined comp_task_num\n");
	task_num = comp_task_num*5;
	task_list = (Ttask_p*) malloc(task_num*sizeof(Ttask_p));

	/// Tile maps moved to autotuner in PARALiA 3.0 (removed from DataTile)
	A_tile_loc_map = (int***) malloc(Grid_M*sizeof(int**));
	for(int im = 0; im < Grid_M; im++){
		A_tile_loc_map[im] = (int**) malloc(Grid_K*sizeof(int*));
		for(int ik = 0; ik < Grid_K; ik++){
			A_tile_loc_map[im][ik] = (int*) malloc(CHL_MEMLOCS*sizeof(int));
			for(int loc = 0; loc < CHL_MEMLOCS; loc++) 
				if (loc == A_loc) A_tile_loc_map[im][ik][loc] = 0;
				else A_tile_loc_map[im][ik][loc] = -42;
		}
	}
	B_tile_loc_map = (int***) malloc(Grid_K*sizeof(int**));
	for(int ik = 0; ik < Grid_K; ik++){
		B_tile_loc_map[ik] = (int**) malloc(Grid_N*sizeof(int*));
		for(int in = 0; in < Grid_N; in++){
			B_tile_loc_map[ik][in] = (int*) malloc(CHL_MEMLOCS*sizeof(int));
			for(int loc = 0; loc < CHL_MEMLOCS; loc++) 
				if (loc == A_loc) B_tile_loc_map[ik][in][loc] = 0;
				else B_tile_loc_map[ik][in][loc] = -42;
		}
	}
	C_tile_loc_map = (int***) malloc(Grid_M*sizeof(int**));
	for(int im = 0; im < Grid_M; im++){
		C_tile_loc_map[im] = (int**) malloc(Grid_N*sizeof(int*));
		for(int in = 0; in < Grid_N; in++){
			C_tile_loc_map[im][in] = (int*) malloc(CHL_MEMLOCS*sizeof(int));
			for(int loc = 0; loc < CHL_MEMLOCS; loc++) 
				if (loc == A_loc) C_tile_loc_map[im][in][loc] = 0;
				else C_tile_loc_map[im][in][loc] = -42;
		}
	}
}

void ATC::optimize_tasks_serial(){
	long int comp_task_ctr, task_ctr = comp_task_ctr = 0, comp_task_perdev[active_unit_num] = {0};
	while (comp_task_ctr < comp_task_num){
		for(int ik = 0; ik < Grid_K; ik++){
			for(int dev_idx = 0; dev_idx < active_unit_num; dev_idx++){
				if(comp_task_perdev[dev_idx] == comp_task_per_unit_num[dev_idx]) continue;
				int dev_id = active_unit_id_list[dev_idx]; 
				int comp_task_idx = comp_task_per_unit_list[dev_idx][comp_task_perdev[dev_idx]++]/Grid_K;
				int im = comp_task_idx/Grid_N, in = comp_task_idx%Grid_N;
				if(A_tile_loc_map[im][ik][dev_id]!= 0 && A_tile_loc_map[im][ik][dev_id]!= 42){
					A_tile_loc_map[im][ik][dev_id] = 2; 
					long int size = T*T*elemSize;
					LinkRoute_p A_tile_route = new LinkRoute();
					A_tile_route->optimize(A_tile_loc_map[im][ik], size);
					task_list[task_ctr++] = new TileTask(FETCH, 0, im, ik, 0, A_tile_route);
				}
				if(B_tile_loc_map[ik][in][dev_id] && B_tile_loc_map[ik][in][dev_id]!= 42){
					B_tile_loc_map[ik][in][dev_id] = 2; 
					long int size = T*T*elemSize;
					LinkRoute_p B_tile_route = new LinkRoute();
					B_tile_route->optimize(B_tile_loc_map[ik][in], size);
					task_list[task_ctr++] = new TileTask(FETCH, 1, ik, in, 0, B_tile_route);
				}
				LinkRoute_p C_tile_route = new LinkRoute();
				if(C_tile_loc_map[im][in][dev_id] && C_tile_loc_map[im][in][dev_id]!= 42){
					C_tile_loc_map[im][in][dev_id] = 2; 
					long int size = T*T*elemSize;
					C_tile_route->optimize(C_tile_loc_map[im][in], size);
					if(!strcmp(OUTPUT_ALGO_MODE, "ALGO_WR"))
						task_list[task_ctr++] = new TileTask(FETCH, 2, im, in, 0, C_tile_route);
				}
				task_list[task_ctr++] = new TileTask(COMPUTE, 2, im, in, ik, C_tile_route);
				if(ik == Grid_K - 1 && C_tile_loc_map[im][in][dev_id]){
					long int size = T*T*elemSize;
					LinkRoute_p C_tile_out_route = new LinkRoute();
					C_tile_out_route->optimize_reverse(C_tile_loc_map[im][in], size);
					task_list[task_ctr++] = new TileTask(WRITEBACK, 2, im, in, 0, C_tile_out_route);
				}
				comp_task_ctr++;
			}
		}
	}
	task_num = task_ctr;
}

void ATC::optimize_tasks(){
	if(!strcmp(TASK_ORDER, "SERIAL")) optimize_tasks_serial();
}

double ATC::autotune_problem(int A_loc_in, int B_loc_in, int C_loc_in, int D_loc_in, 
    int M_in, int N_in, int K_in, int elemSize){
	short lvl = 3;
	double cpu_timer = csecond();
#ifdef DEBUG
	fprintf(stderr,  "|-----> ATC::autotune_problem(%d, %d, %d, %d, %d, %d, %d, %d)\n", 
		A_loc_in, B_loc_in, C_loc_in, D_loc_in, M_in, N_in, K_in, elemSize);
	print();
#endif
	A_loc = A_loc_in;
	B_loc = B_loc_in;
	C_loc = C_loc_in;
	D_loc = D_loc_in;
	M = M_in;
	N = N_in;
	K = K_in;
	int autotune_eval_devices = 0;
	if (active_unit_num > 0){
		if (active_unit_id_list){
//#ifdef PDEBUG
		fprintf(stderr,  "Running on %d devices with dev_ids=[ ", active_unit_num);
		for (int i =0; i < active_unit_num; i++) fprintf(stderr, "%d ", active_unit_id_list[i]);
		fprintf(stderr, "]\n");
//#endif
		}
		else{
			autotune_eval_devices = 1;
//#ifdef PDEBUG
		fprintf(stderr,  "Running on %d devices with tunable dev_ids\n", active_unit_num);
//#endif
		}
	}
	else{
		active_unit_num = CHL_WORKERS;
		autotune_eval_devices = 1;
		for (int i =0; i < active_unit_num; i++)
			active_unit_id_list[i] = i;
	}
	active_memloc_num = 0;
	for( int idx = 0; idx < CHL_MEMLOCS; idx++){
		if (is_in_list(idx, active_unit_id_list, active_unit_num) 
		|| idx == A_loc || idx == B_loc || idx == C_loc || idx == D_loc) 
			active_memlocs[active_memloc_num++] = idx;
	}
	if (autotune_eval_devices){
		error("ATC::autotune_problem: Not implemented (yet) for autotune_eval_devices\n");
	}
	else{
		int initial_T = T;
		double tile_selection_t = 0, split_selection_t = 0, best_t = DBL_MAX;
		Gamalg_p test_grid;
		for (int idx = 0; idx < MAX_WORKER_CONFIG; idx++){
			test_grid = new Grid_amalgamation(CHL_INPUT_QUEUES_CASE_IDS[active_unit_num-1][idx]);
			if(!test_grid->load_edges(CHL_INPUT_QUEUES_CASE_IDS[active_unit_num-1][idx], 
				CHL_OUTPUT_QUEUES_CASE_IDS[active_unit_num-1][idx])) continue;
			if(test_grid->active_nodes_id != translate_unit_list_to_binary(active_unit_id_list,active_unit_num)) continue;
			for (int idx = 0; idx < active_unit_num; idx++) active_unit_score[idx] = 1.0/active_unit_num;
			long long edge_load[64][64];
			gemm_translate_problem_comm(edge_load, A_loc, B_loc, C_loc, D_loc, M, N, K, elemSize, active_unit_num, active_unit_id_list, active_unit_score);
			test_grid->set_edge_load(edge_load);
			test_grid->update_problem_edges();
			test_grid->load_nodes();
			long long node_ops[CHL_WORKERS], node_mem_ops[CHL_WORKERS];
			gemm_translate_problem_ops(node_ops, node_mem_ops, M, N, K, active_unit_num, active_unit_id_list, active_unit_score);
			test_grid->set_node_load(node_ops, node_mem_ops);
			double temp_t = test_grid->get_problem_perf_estimation();
			if(temp_t <= best_t){
				inter_grid = test_grid;
				best_t = temp_t; 
			}
			else{
				delete test_grid; 
			}
		}
		if (best_t == DBL_MAX) error("ATC::autotune_problem: No device configuration found matching devices %s\n",
		printlist(active_unit_id_list,active_unit_num));
		if(initial_T <= 0) tile_selection_t += optimize_tile();
		else{
			double* c_T_sl = (double*) calloc(6,sizeof(double));
			get_T_slowdowns(c_T_sl, initial_T);
			set_T_slowdowns(c_T_sl);
			free(c_T_sl);
		}		// TODO: Must decide if workload ratio should be tuned when there is a predefined number of devices... Currently == off for paper
		split_homogeneously = 1;
	}
	
	for(int i = 0; i< CHL_MEMLOCS; i++)	for(int j = 0; j< CHL_MEMLOCS; j++){
			best_grid_edge_bws[i][j] = inter_grid->problem_edge_bw[i][j];
			best_grid_edge_active[i][j] = inter_grid->edge_active[i][j];
			best_grid_edge_replaced[i][j][0] = inter_grid->edge_replaced[i][j][0];
			best_grid_edge_replaced[i][j][1] = inter_grid->edge_replaced[i][j][1];
	}
#ifdef SDEBUG
	inter_grid->print_edge_active();
	inter_grid->print_edge_replaced();
	inter_grid->print_edge_bws();
	inter_grid->print_simu_edge_bws();
	inter_grid->print_edge_load();
	inter_grid->print_problem_edge_bws();
	inter_grid->print_nodes();
	inter_grid->print_node_load();
#endif
	Grid_M = M/T + ((M%T) ? 1 : 0);
	Grid_N = N/T + ((N%T) ? 1 : 0);
	Grid_K = K/T + ((K%T) ? 1 : 0);

	// Calculate compute tasks and allocate comp_task_unit_list
	update_comp_task_num(Grid_M*Grid_N*Grid_K);
	// Distribute compute tasks to devices
	distribute_comp_tasks(); 
	// Calculate (max) total tasks and allocate task_list
	initialize_tasks();
	// Translate compute tasks to sub-tasks based on 1) an ordering algorithm and 2) a routing algorithm.
	optimize_tasks();

	cpu_timer = csecond() - cpu_timer;
	if(T!=-1){
		if (T_imbalance_sl > 0) 
			warning("ATC::optimize_tile -> T = %ld: C1 (NO-imbalance) was not satisfied, estimated sl = %lf\n", 
				T, T_imbalance_sl);
		if (T_remainder_sl > 0) 
			warning("ATC::optimize_tile -> T = %ld: C2 (NO-remainder) was not satisfied, estimated sl = %lf\n", 
				T, T_remainder_sl);
		if (T_small_sl > 0) 
			warning("ATC::optimize_tile -> T = %ld: C3 (T >= %d) was not satisfied, estimated sl = %lf\n", 
				T, TILE_MIN, T_small_sl);
		double sl_too_many_tasks = 0;
		if (comp_task_num/active_unit_num > MAX_DESIRED_SK_DEV){
			sl_too_many_tasks = (1.0*comp_task_num/active_unit_num/MAX_DESIRED_SK_DEV)*MAX_DESIRED_SK_DEV_SLOWDOWN;
			warning("ATC::optimize_tile -> T = %ld: C4 (SK_DEV <= %d) was not satisfied, estimated sl = %lf\n", 
				T, MIN_DESIRED_SK_DEV, sl_too_many_tasks);
		}
		fprintf(stderr, "ATC::optimize_tile -> T = %ld: estimated sl from overlap = %lf\n", 
			T, T_sknum_sl - sl_too_many_tasks);
		if (T_big_sl > 0 ) 
			warning("ATC::optimize_tile -> T = %ld: C5 (T <= %d) was not satisfied, estimated sl = %lf\n", 
				T, TILE_MAX, T_big_sl);
	}
	fprintf(stderr, "====================================\n");
	fprintf(stderr, "ATC::autotune_problem: Autotuning complete-> t_autotune = %lf ms\n", cpu_timer*1000);
	fprintf(stderr, "autotune_controller: T=%ld, T_aggregate_sl=%lf active_unit_num=%d, Problem split = %s -> %s :\n"
		"\t -> pred_t = %lf ms, pred_J = %lf kJ, pred_PDP = %lf Gflops/J, pred_EDP = %lf Gflops^2/J\n"
		"\t -> pred_t_pesimistic = %lf ms, pred_J_pesimistic = %lf kJ, PDP_pesimistic = %lf Gflops/J, EDP_pesimistic = %lf Gflops^2/J\n",
		T, T_aggregate_sl, active_unit_num, printlist<int>(active_unit_id_list, active_unit_num),
		printlist<long int>(comp_task_per_unit_num, active_unit_num), pred_t*1000, pred_J/1000, power_delay, energy_delay,
		pred_t_pesimistic*1000, pred_J_pesimistic/1000, power_delay_pesimistic, energy_delay_pesimistic);
	fprintf(stderr, "====================================\n");
#ifdef DEBUG
	fprintf(stderr,  "<-----|\n");
#endif
	return cpu_timer;
}

void ATC::set_T_slowdowns(double* slowdowns){
	T_aggregate_sl = slowdowns[0];
	T_imbalance_sl = slowdowns[1];
	T_remainder_sl = slowdowns[2];
	T_small_sl = slowdowns[3];
	T_sknum_sl = slowdowns[4];
	T_big_sl = slowdowns[5];
}

void ATC::get_T_slowdowns(double* slowdown, int candidate_T){
	for(int idx = 0; idx < 6; idx++) slowdown[idx] = 0.0;
	long long sk_num = (M/candidate_T + ((M%candidate_T) ? 1 : 0))*
		(N/candidate_T + ((N%candidate_T) ? 1 : 0))*(K/candidate_T + ((K%candidate_T) ? 1 : 0));
	// Condition 1
	int bucket_num = sk_num/(K/candidate_T + ((K%candidate_T) ? 1:0));
	slowdown[1] = (bucket_num / active_unit_num ) ? 
		1.0*(bucket_num % active_unit_num)/(bucket_num / active_unit_num) * active_unit_num 
		: 1.0* bucket_num / active_unit_num;
	// Condition 2
	if(M%candidate_T) slowdown[2] += 1.0/(M/candidate_T);
	if(N%candidate_T) slowdown[2] += 1.0/(N/candidate_T);
	if(K%candidate_T) slowdown[2] += 1.0/(K/candidate_T);
	// Condition 3
	if(candidate_T < TILE_MIN) slowdown[3]+= 1.0*TILE_MIN/candidate_T*TILE_MIN_SLOWDOWN;
	// Condition 4.1
	long int dev_sks = (1.0*sk_num)/active_unit_num; 
	slowdown[4]+= 2.0/(dev_sks); // This slowdown like this removes the need for MIN_DESIRED_SK_DEV
	//if(dev_sks < MIN_DESIRED_SK_DEV) slowdown+= 1/dev_sks;
	// Condition 4.2
	if(dev_sks > MAX_DESIRED_SK_DEV) slowdown[4]+= (1.0*dev_sks/MAX_DESIRED_SK_DEV)*MAX_DESIRED_SK_DEV_SLOWDOWN;
	// Condition 5
	if(candidate_T > TILE_MAX) slowdown[5]+=candidate_T/TILE_MAX*TILE_MAX_SLOWDOWN;
	slowdown[0] = slowdown[1] + slowdown[2] + slowdown[3] + slowdown[4]  + slowdown[5];
#ifdef DPDEBUG
	fprintf(stderr,  "====================================\n");
	fprintf(stderr,  "ATC::get_T_slowdowns(D1=%d, D2 = %d, D3 = %d) T=%d with T_aggregate_sl = %lf, T_imbalance_sl= %lf, T_remainder_sl= %lf, T_small_sl= %lf, "
	"T_sknum_sl= %lf, T_big_sl = %lf\n", M, N, K, candidate_T, slowdown[0], slowdown[1], slowdown[2], slowdown[3], slowdown[4], slowdown[5]);
#endif
	return;
}

double ATC::optimize_tile(){
	double timer = csecond();
#ifdef DEBUG
fprintf(stderr,  "|-----> ATC::optimize_tile( autotune_controller{ T=%ld, active_unit_num=%d, Problem split = %s -> %s : t_pred = %lf ms})\n",
	T, active_unit_num, printlist<int>(active_unit_id_list, active_unit_num),
	printlist<double>(active_unit_score, active_unit_num), pred_t*1000);
#endif
	int best_idx = -1;
	double temp_score = 0;

	if (active_unit_num <= 0)
	error("ATC::optimize_tile: Called with active_unit_num = %d\n", active_unit_num);
	int max_allowed_T = std::min(M, std::min(N, K));
	int best_T = -1;
	double* best_T_sl = (double*) calloc(6,sizeof(double));
	for(int idx = 0; idx < 6; idx++) best_T_sl[idx] = DBL_MAX;
	double* c_T_sl = (double*) calloc(6,sizeof(double));
	for (int candidate_T = max_allowed_T; candidate_T > 0; candidate_T--)
	// Condition 1
	//if(!((M/(candidate_T) + ((M%(candidate_T))? 1:0))
	// *(N/(candidate_T) + ((N%(candidate_T))? 1:0)) % active_unit_num))
	{ 
		get_T_slowdowns(c_T_sl, candidate_T); 
		if (c_T_sl[0] < best_T_sl[0]){
			for(int idx = 0; idx < 6; idx++) best_T_sl[idx] = c_T_sl[idx];
			best_T = candidate_T;
		}
	}
	T = best_T;
	set_T_slowdowns(best_T_sl);
	free(best_T_sl);
	free(c_T_sl);
#ifdef PDEBUG
	fprintf(stderr,  "====================================\n");
	fprintf(stderr,  "Predict T=%ld with T_aggregate_sl = %lf, T_remainder_sl= %lf, T_small_sl= %lf, "
	"T_sknum_sl= %lf, T_big_sl = %lf\n", T, T_aggregate_sl, T_remainder_sl, T_small_sl, T_sknum_sl, T_big_sl);
#endif
	timer = csecond() - timer;
#ifdef TEST
	fprintf(stderr,  "Tile selection time:%lf ms\n", timer*1000);
#endif
#ifdef DEBUG
	fprintf(stderr,  "<-----|\n");
#endif
	return timer;
}

/******************************************************************************/
/**************************** Helper Fuctions *********************************/
void ATC::print(){
	//int dev_ids_token = 0;
	//int ctr = 0, itter = 0;
	//if (active_unit_num > 0) for (int i = 0; i < active_unit_num; i++) dev_ids_token+=pow(10,idxize(active_unit_id_list[i]));
	fprintf(stderr, "Autotune controller:\n->T = %ld\n->active_unit_num = %d\n->active_unit_id_list = %s\n->active_unit_score = %s\
	\n->pred_t = %lf\n->comp_task_num = %ld\n->comp_task_per_unit_num = %s\n->task_num = %ld", T, active_unit_num,
	printlist<int>(active_unit_id_list, active_unit_num),
	printlist<double>(active_unit_score, active_unit_num),
	pred_t, comp_task_num,
	printlist<long int>(comp_task_per_unit_num, active_unit_num), task_num);
	fprintf(stderr, "\n");
	return;
}

const char* ATC::print_csv(){
	char* outstring = (char*) malloc(256*sizeof(char));
	int dev_ids_token = 0;
	int ctr = 0, itter = 0;
	if (active_unit_num > 0) for (int i = 0; i < active_unit_num; i++) dev_ids_token+=std::pow(10,active_unit_id_list[i]);
	sprintf(outstring, "%ld,%d,%d,%lld",  T, active_unit_num, dev_ids_token, cache_limit);
	return outstring;
}