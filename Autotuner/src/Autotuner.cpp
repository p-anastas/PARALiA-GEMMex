///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief The Autotuner controller functions.
///

#include "Autotuner.hpp"
#include "Subkernel_distributions.hpp"

#include <cfloat> /// For DBL_MAX
#include <cmath>

double best_grid_edge_bws[64][64];
int best_grid_edge_active[64][64];
int best_grid_edge_replaced[64][64][2];

/********************** Initialization/Modification ***************************/
ATC::ATC(){
	short lvl = 2;
#ifdef DEBUG
		fprintf(stderr,  "|-----> ATC::ATC\n");
#endif
	active_unit_id_list = (int*) malloc(CHL_WORKERS*sizeof(int));
	active_unit_score = (double*) malloc(CHL_WORKERS*sizeof(double));
	Subkernels_per_unit_num = (int*) malloc(CHL_WORKERS*sizeof(int));
	Subkernels_per_unit_list = (int**) malloc(CHL_WORKERS*sizeof(int*));
	for (int d = 0; d < CHL_WORKERS; d++){
		active_unit_score[d] = -42.0;
		Subkernels_per_unit_list[d] = NULL;
		Subkernels_per_unit_num[d] = 0;
	}
	T = active_unit_num = subkernel_num = -1;
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
	free(active_unit_id_list);
	free(active_unit_score);
	free(Subkernels_per_unit_num);
	for (int d = 0; d < CHL_WORKERS; d++) free(Subkernels_per_unit_list[d]);
	free(Subkernels_per_unit_list);
#ifdef DEBUG
	fprintf(stderr,  "<-----|\n");
#endif
}

void ATC::reset(){
	short lvl = 4;
#ifdef DEBUG
		fprintf(stderr,  "|-----> ATC::reset\n");
#endif
	for (int d = 0; d < CHL_WORKERS; d++){
		free(Subkernels_per_unit_list[d]);
		Subkernels_per_unit_list[d] = NULL;
		Subkernels_per_unit_num[d] = 0;
	}
	T = active_unit_num = subkernel_num = -1;
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
	short lvl = 3;
	#ifdef DEBUG
		fprintf(stderr,  "|-----> ATC::mimic_ATC(other_ATC = %p)\n", other_ATC);
	#endif
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
	inter_grid->copy(other_ATC->inter_grid);

	if(subkernel_num != -1){
		for (int d = 0; d < CHL_WORKERS; d++){
			//fprintf(stderr,"Subkernels_per_unit_list[%d] = %p\n", d, Subkernels_per_unit_list[d]);
			//free(Subkernels_per_unit_list[d]);
			//Subkernels_per_unit_num[d] = 0;
			;//Subkernels_per_unit_list[d] = NULL;
		} /// TODO: Got some "double free 2cache" error when used in many different mimicked ATCs ->
			/// potential problem here  ATC::update_sk_num resizing Subkernels_per_unit_list[d] might be solution and/or relevant.
		subkernel_num = -1;
	}

	if (other_ATC->subkernel_num != -1){
		for (int d = 0; d < other_ATC->active_unit_num; d++){
			Subkernels_per_unit_num[d] = other_ATC->Subkernels_per_unit_num[d];
			free(Subkernels_per_unit_list[d]);
			Subkernels_per_unit_list[d] = (int*) malloc(other_ATC->subkernel_num*sizeof(int));
			for (int sk = 0; sk < other_ATC->subkernel_num; sk++)
				Subkernels_per_unit_list[d][sk] = other_ATC->Subkernels_per_unit_list[d][sk];
		}
		subkernel_num = other_ATC->subkernel_num;
	}

#ifdef DEBUG
	fprintf(stderr,  "<-----|\n");
#endif
}

void ATC::update_sk_num(long long int subkernel_num_in){
	short lvl = 3;
	#ifdef DEBUG
		fprintf(stderr,  "|-----> ATC::update_sk_num\n");
	#endif
	int prev_sk_num = subkernel_num;
	subkernel_num = subkernel_num_in;
	if (prev_sk_num == -1)  for (int d = 0; d < CHL_WORKERS; d++){
		Subkernels_per_unit_list[d] = (int*) malloc(subkernel_num*sizeof(int));
		for (int sk = 0; sk < subkernel_num; sk++) Subkernels_per_unit_list[d][sk] = -1;
	}
	else if (prev_sk_num < subkernel_num) for (int d = 0; d < CHL_WORKERS; d++){
		free(Subkernels_per_unit_list[d]);
		Subkernels_per_unit_list[d] = (int*) malloc(subkernel_num*sizeof(int));
		for (int sk = 0; sk < subkernel_num; sk++) Subkernels_per_unit_list[d][sk] = -1;
	}
#ifdef DEBUG
	fprintf(stderr,  "<-----|\n");
#endif
}

void ATC::distribute_subkernels(int D1GridSz, int D2GridSz, int D3GridSz){
	if (!strcmp(DISTRIBUTION, "ROUND-ROBIN"))
		CoCoDistributeSubkernelsRoundRobin(this);
	else if (!strcmp(DISTRIBUTION, "SPLIT-NAIVE"))
		CoCoDistributeSubkernelsNaive(this);
	else if (!strcmp(DISTRIBUTION, "SPLIT-CHUNKS-ROBIN"))
		CoCoDistributeSubkernelsRoundRobinChunk(this, D3GridSz);
	else if (!strcmp(DISTRIBUTION, "SPLIT-CHUNKS-ROBIN-REVERSE"))
		CoCoDistributeSubkernelsRoundRobinChunkReverse(this, D3GridSz);
	else if (!strcmp(DISTRIBUTION, "2D-BLOCK-CYCLIC"))
		CoCoDistributeSubkernels2DBlockCyclic(this, D1GridSz, D2GridSz, D3GridSz);
	else error("ATC::distribute_subkernels: Unknown Subkernel Distribution %s\n", DISTRIBUTION);
	for (int i = 0; i < active_unit_num; i++)
	if(!Subkernels_per_unit_num[i]) {
		free(Subkernels_per_unit_list[i]);
		for (int i_move = i; i_move < active_unit_num - 1; i_move++){
			active_unit_score[i_move] = active_unit_score[i_move+1];
			active_unit_id_list[i_move] = active_unit_id_list[i_move+1];
			Subkernels_per_unit_num[i_move] = Subkernels_per_unit_num[i_move+1];
			Subkernels_per_unit_list[i_move] = Subkernels_per_unit_list[i_move+1];
		}
		i--;
		active_unit_num--;
	}
#ifdef PDEBUG
    print();
#endif
}
/******************************************************************************/
/****************************** Autotuning ************************************/

double ATC::autotune_problem(int A_loc, int B_loc, int C_loc, int D_loc, 
    int M, int N, int K, int elemSize){
	short lvl = 3;
	double cpu_timer = csecond();
#ifdef DEBUG
	fprintf(stderr,  "|-----> ATC::autotune_problem(%s, %p)\n", routine_name,
		initial_problem_wrap);
	print();
#endif
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
	
	if (autotune_eval_devices){
		error("ATC::autotune_problem: Not implemented (yet) for autotune_eval_devices\n");
	}
	else{
		int initial_T = T;
		double tile_selection_t = 0, split_selection_t = 0, best_t = DBL_MAX;
		Gamalg_p test_grid;
		for (int idx = 0; idx < MAX_WORKER_CONFIG; idx++){
			test_grid = new Grid_amalgamation(CHL_INPUT_QUEUES_CASE_IDS[active_unit_num][idx]);
			test_grid->load_edges(CHL_INPUT_QUEUES_CASE_IDS[active_unit_num][idx], 
				CHL_OUTPUT_QUEUES_CASE_IDS[active_unit_num][idx]);
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
			if(temp_t < best_t){
				inter_grid = test_grid;
				best_t = temp_t; 
			}
			else{
				delete test_grid; 
			}
		}
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
	}
#ifdef SDEBUG
	inter_grid->print_edge_active();
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
	update_sk_num(Grid_M*Grid_N*Grid_K);
	distribute_subkernels(Grid_M, Grid_N, Grid_K);

	cpu_timer = csecond() - cpu_timer;
	if(T!=-1){
		if (T_imbalance_sl > 0) 
			warning("ATC::optimize_tile -> T = %d: C1 (NO-imbalance) was not satisfied, estimated sl = %lf\n", 
				T, T_imbalance_sl);
		if (T_remainder_sl > 0) 
			warning("ATC::optimize_tile -> T = %d: C2 (NO-remainder) was not satisfied, estimated sl = %lf\n", 
				T, T_remainder_sl);
		if (T_small_sl > 0) 
			warning("ATC::optimize_tile -> T = %d: C3 (T >= %d) was not satisfied, estimated sl = %lf\n", 
				T, TILE_MIN, T_small_sl);
		double sl_too_many_sk = 0;
		if (subkernel_num/active_unit_num > MAX_DESIRED_SK_DEV){
			sl_too_many_sk = (1.0*subkernel_num/active_unit_num/MAX_DESIRED_SK_DEV)*MAX_DESIRED_SK_DEV_SLOWDOWN;
			warning("ATC::optimize_tile -> T = %d: C4 (SK_DEV <= %d) was not satisfied, estimated sl = %lf\n", 
				T, MIN_DESIRED_SK_DEV, sl_too_many_sk);
		}
		fprintf(stderr, "ATC::optimize_tile -> T = %d: estimated sl from overlap = %lf\n", 
			T, T_sknum_sl - sl_too_many_sk);
		if (T_big_sl > 0 ) 
			warning("ATC::optimize_tile -> T = %d: C5 (T <= %d) was not satisfied, estimated sl = %lf\n", 
				T, TILE_MAX, T_big_sl);
	}
	fprintf(stderr, "====================================\n");
	fprintf(stderr, "ATC::autotune_problem: Autotuning complete-> t_autotune = %lf ms\n", cpu_timer*1000);
	fprintf(stderr, "autotune_controller: T=%ld, T_aggregate_sl=%lf active_unit_num=%d, Problem split = %s -> %s :\n"
		"\t -> pred_t = %lf ms, pred_J = %lf kJ, pred_PDP = %lf Gflops/J, pred_EDP = %lf Gflops^2/J\n"
		"\t -> pred_t_pesimistic = %lf ms, pred_J_pesimistic = %lf kJ, PDP_pesimistic = %lf Gflops/J, EDP_pesimistic = %lf Gflops^2/J\n",
		T, T_aggregate_sl, active_unit_num, printlist<int>(active_unit_id_list, active_unit_num),
		printlist<int>(Subkernels_per_unit_num, active_unit_num), pred_t*1000, pred_J/1000, power_delay, energy_delay,
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
	if(candidate_T > TILE_MAX) slowdown[5]+=candidate_T/TILE_MAX*TILE_MΑΧ_SLOWDOWN;
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
fprintf(stderr,  "|-----> ATC::optimize_tile( autotune_controller{ T=%ld, active_unit_num=%d, Problem split = %s -> %s : t_pred = %lf ms}, unit_modeler_list =%p)\n",
	T, active_unit_num, printlist<int>(active_unit_id_list, active_unit_num),
	printlist<double>(active_unit_score, active_unit_num), pred_t*1000, unit_modeler_list);
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
	\n->pred_t = %lf\n->subkernel_num = %ld\n->Subkernels_per_unit_num = %s\n", T, active_unit_num,
	printlist<int>(active_unit_id_list, active_unit_num),
	printlist<double>(active_unit_score, active_unit_num),
	pred_t, subkernel_num,
	printlist<int>(Subkernels_per_unit_num, active_unit_num));
 	if(subkernel_num != -1){
	for (int d = 0; d < active_unit_num; d++) fprintf(stderr, "Subkernels_per_unit_list[%d] = %s\n", d,
		printlist<int>(Subkernels_per_unit_list[d], subkernel_num));
	}
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

/******************************************************************************/

/* FIXME: DEPRECATED - kept for future model-based T selection update if needed
double ATC::optimize_tile_modelBased(ATC_p autotune_controller, short used_devs, short* used_dev_ids,
	int* dev_idx_ignore, MD_p* unit_modeler_list){
	short lvl = 3;
#ifdef DEBUG
	fprintf(stderr,  "||-----> ATC::optimize_tile_modelBased(used_devs=%d, used_dev_ids= [ ", used_devs);
	for (int i =0; i < used_devs; i++) fprintf(stderr, "%d ", used_dev_ids[i]);
	fprintf(stderr, "]\n");
#endif
#ifdef TEST
	double timer = csecond();
#endif
	short first_model_idx = (used_dev_ids[0] == -1) ? LOC_NUM - 1 : used_dev_ids[0];
	MD_p model = unit_modeler_list[first_model_idx];
	active_unit_score = NULL;
	int best_idx = -1;
	double temp_score = 0;

	long int min_T = 0, max_allowed_T = 0, ctr = 0;
	max_allowed_T = CoCopeLiaMaxT(model);
	min_T = CoCopeLiaMinT(model);
#ifdef PDEBUG
		fprintf(stderr,  "min_T = %ld, max_allowed_T = %ld\n",
			min_T, max_allowed_T);
#endif
	if (min_T > max_allowed_T){
		outparams->T = max_allowed_T;
		// FIXME: Undefined expected performance for tiles < than the smaller microbenchmark
		outparams->pred_t = 0;
#ifdef PDEBUG
		fprintf(stderr,  "min_T = %ld > max_allowed_T = %ld: returning T = %ld",
			min_T, max_allowed_T, max_allowed_T);
#endif
		return outparams;
	}
	double temp_t, min_overlap_t = 10000000;
	long int prev_trial_T = 0;

	int lines = CoCoPeLiaGPUexecGetLines(model);
	for (ctr = 0 ; ctr < lines ; ctr++){
		long int trial_T = CoCoPeLiaGPUexecGetElem(model, ctr);
		if (trial_T > max_allowed_T) break;
		if (trial_T ==  prev_trial_T) continue;

		double temp_overlap_t = 0;
		for(int idx = 0; idx < used_devs; idx++){
			if(dev_idx_ignore[idx]) continue;
			short cur_dev_id = used_dev_ids[idx], cur_dev_idx = (cur_dev_id == -1)? LOC_NUM - 1 : cur_dev_id;
			model = unit_modeler_list[cur_dev_idx];
			ModelType used_model;
			switch(model->problem){
				case BLAS1:
					used_model = COCOPELIA_BIDIRECTIONAL;
					break;
				case BLAS2:
					used_model = COCOPELIA_BIDIRECTIONAL;
					break;
				case BLAS3:
					used_model = COCOPELIA_REUSE;
					break;
				default:
					error("ATC::optimize_tile_modelBased:\
					model->problem switch default reached\n");
			}
				temp_t = CoCoPeLiaModelPredict(model, trial_T, used_model);
				double imb_time_multiplier = 1.0, reduce_time_multiplier = 1.0;
#ifdef TILE_IMBALANCE_PENALTY
					if (model->D1 != -1 && (model->D1%trial_T)) imb_time_multiplier+=TILE_IMBALANCE_PENALTY;
					if (model->D2 != -1 && (model->D2%trial_T)) imb_time_multiplier+=TILE_IMBALANCE_PENALTY;
					if (model->D3 != -1 && (model->D3%trial_T)) imb_time_multiplier+=TILE_IMBALANCE_PENALTY;
#endif
#ifdef REDUCE_PENALTY
					if ((model->D1/trial_T + ((model->D1%trial_T)? 1 : 0))*(model->D2/trial_T + ((model->D2%trial_T)? 1 : 0))
						*(model->D3/trial_T + ((model->D3%trial_T)? 1 : 0)) % used_devs) reduce_time_multiplier+=REDUCE_PENALTY;
#endif
			temp_t *= imb_time_multiplier *reduce_time_multiplier;
			if(temp_t > 0) temp_overlap_t = fmax(temp_overlap_t, temp_t);
			else error("CoCoPeLiaModelPredict(%p(dev_id = %d, (idx = %d )), trial_T = %ld): negative prediction temp_t = %lf\n",
				model, cur_dev_id, cur_dev_idx, trial_T, temp_t);
#ifdef PDEBUG
			fprintf(stderr,  "CoCoPeLiaModelPredict(%p) for dev_id = %d (idx = %d ) with trial_T = %ld: temp_overlap_t = %lf, temp_t = %lf\n",
				model, cur_dev_id, cur_dev_idx, trial_T, temp_overlap_t, temp_t);
#endif
		}
		if (temp_overlap_t < min_overlap_t){
			min_overlap_t = temp_overlap_t;
			min_T = trial_T;
		}
		prev_trial_T = trial_T;
	}
	outparams->T = min_T;
	outparams->pred_t = min_overlap_t;
#ifdef PDEBUG
	fprintf(stderr,  "====================================\n");
	fprintf(stderr,  "Predict T=%ld : t_pred = %lf\n", outparams->T, outparams->pred_t);
#endif
#ifdef TEST
	timer = csecond() - timer;
	fprintf(stderr,  "Tile selection time:%lf ms\n", timer*1000);
	fprintf(stderr,  "<-----|\n");
#endif
#ifdef DEBUG
	fprintf(stderr,  "outparams->T = %ld\n : outparams->pred_t = %lf ms\n", outparams->T, outparams->pred_t);
	fprintf(stderr,  "<-----|\n");
#endif
	return outparams;
}
*/
