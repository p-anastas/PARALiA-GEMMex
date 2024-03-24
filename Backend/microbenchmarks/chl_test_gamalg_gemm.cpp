///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief A transfer microbenchmark aiming to capture how link overlap effects their bandwidth
///

#include <unistd.h>
#include <cassert>

#include <numa.h>

#include "smart_wrappers.hpp"
#include "grid_amalgamation.hpp"
#include "microbenchmarks.hpp"

int main(const int argc, const char *argv[]) {

	int ctr = 1, test_case_id = -1, rev_test_case_id;
	int M = 1, N = 1, K = 1;
	int A_loc, B_loc, C_loc, D_loc = A_loc = B_loc = C_loc = CHL_MEMLOCS - 1;

	switch (argc) {
	case (10):
		M = atoi(argv[ctr++]);
		N = atoi(argv[ctr++]);
		K = atoi(argv[ctr++]);
		test_case_id = atoi(argv[ctr++]);
		rev_test_case_id = atoi(argv[ctr++]);
		A_loc = atoi(argv[ctr++]);
		B_loc = atoi(argv[ctr++]);
		C_loc = atoi(argv[ctr++]);
		D_loc = atoi(argv[ctr++]);
		break;
	default:
		error("Incorrect input arguments. Usage: ./correct_run M N K test_case_id rev_test_case_id Aloc Bloc Cloc Dloc:\n"
		"M, N, K : GEMM problem dimensions\n"
		"test_case_id, rev_test_case_id: The sets of active H2D and D2H devices. -1 for itterative search of best combination.\n"
		"Aloc, Bloc, Cloc, Dloc: location of matrices in MEM_LOC:"
		"-> If set to -1, matrix Xloc is consider equally split in host numa nodes\n"
		"-> If set to -2, matrix Xloc is consider optimaly split in host numa nodes\n");
  	}

    int active_unit_num, active_unit_id_list[CHL_WORKERS] = {-1};
	double active_unit_id_scores[CHL_WORKERS] = {0};
	//C_loc = 0; 
	int elemSize = 8, best_idx = 0; 
	double best_t = DBL_MAX, temp_t; 
	if(test_case_id > 0)
	{	
		Gamalg_p test_grid = new Grid_amalgamation(test_case_id);
		test_grid->load_edges(test_case_id, rev_test_case_id);
		test_grid->print_edge_active();
		test_grid->print_edge_bws();
		test_grid->print_simu_edge_bws();
		translate_binary_to_unit_list(test_grid->active_nodes_id, &active_unit_num, active_unit_id_list);
		for (int idx = 0; idx < active_unit_num; idx++) active_unit_id_scores[idx] = 1.0/active_unit_num;
		long long edge_load[64][64];
		gemm_translate_problem_comm(edge_load, A_loc, B_loc, C_loc, D_loc, M, N, K, elemSize, active_unit_num, active_unit_id_list, active_unit_id_scores);
		test_grid->set_edge_load(edge_load);
		test_grid->print_edge_load();
		test_grid->update_problem_edges();
		test_grid->print_problem_edge_bws();
		test_grid->load_nodes();
		test_grid->print_nodes();
		long long node_ops[CHL_WORKERS], node_mem_ops[CHL_WORKERS];
		gemm_translate_problem_ops(node_ops, node_mem_ops, M, N, K, active_unit_num, active_unit_id_list, active_unit_id_scores);
		test_grid->set_node_load("MM_FP64", node_ops, node_mem_ops);
		test_grid->print_node_load();
		temp_t = test_grid->get_problem_perf_estimation();
		fprintf(stderr, "Performance estimation total_t = %lf ms (%.2lf Tops/s)\n", 1000*temp_t, Gval_per_s(gemm_ops(M, N, K), temp_t));

	}
	else if(test_case_id == -1){
		fprintf(stderr, "ch_load_gamalg: Loading all stored grid amalgamations...\n");
		system_gamalg_init_from_DB();
		fprintf(stderr, "ch_load_gamalg: Loaded %d grid amalgamations\n", system_gamalg_ctr);
		//int 
		for (int ctr = 0; ctr < system_gamalg_ctr; ctr++){
#ifdef CLBEBUG
			fprintf(stderr,"------------------------------------------------------------------------------------------------------------------------------------------------------\n");
			fprintf(stderr, "ch_load_gamalg: Amalgamation %d (active_nodes_id = %d)\n:", ctr, system_gamalg[ctr]->active_nodes_id);
			system_gamalg[ctr]->print_edge_active();
			system_gamalg[ctr]->print_edge_bws();
			system_gamalg[ctr]->print_simu_edge_bws();
#endif
			translate_binary_to_unit_list(system_gamalg[ctr]->active_nodes_id, &active_unit_num, active_unit_id_list);
			for (int idx = 0; idx < active_unit_num; idx++) active_unit_id_scores[idx] = 1.0/active_unit_num;
			long long edge_load[64][64];
			gemm_translate_problem_comm(edge_load, A_loc, B_loc, C_loc, D_loc, M, N, K, elemSize, active_unit_num, active_unit_id_list, active_unit_id_scores);
			system_gamalg[ctr]->set_edge_load(edge_load);
			system_gamalg[ctr]->update_problem_edges();
			system_gamalg[ctr]->load_nodes();
			long long node_ops[CHL_WORKERS], node_mem_ops[CHL_WORKERS];
			gemm_translate_problem_ops(node_ops, node_mem_ops, M, N, K, active_unit_num, active_unit_id_list, active_unit_id_scores);
			system_gamalg[ctr]->set_node_load("MM_FP64", node_ops, node_mem_ops);
			temp_t= system_gamalg[ctr]->get_problem_perf_estimation();
			if(temp_t < best_t){
				best_t = temp_t; 
				best_idx = ctr;
			}
#ifdef CLBEBUG
			system_gamalg[ctr]->print_edge_load();
			system_gamalg[ctr]->print_problem_edge_bws();
			system_gamalg[ctr]->print_nodes();
			system_gamalg[ctr]->print_node_load();
			fprintf(stderr, "Performance estimation total_t = %lf ms (%.2lf Tops/s)\n", 1000*temp_t, Gval_per_s(gemm_ops(M, N, K), temp_t));
			fprintf(stderr,"------------------------------------------------------------------------------------------------------------------------------------------------------\n");
#endif
		}
		fprintf(stderr, "ch_load_gamalg: Selected amalgamation %d (active_nodes_id = %d)\n:", best_idx, system_gamalg[best_idx]->active_nodes_id);
		system_gamalg[best_idx]->print_edge_active();
		system_gamalg[best_idx]->print_edge_bws();
		system_gamalg[best_idx]->print_simu_edge_bws();
		system_gamalg[best_idx]->print_edge_load();
		system_gamalg[best_idx]->print_problem_edge_bws();
		system_gamalg[best_idx]->print_nodes();
		system_gamalg[best_idx]->print_node_load();
		fprintf(stderr, "Performance estimation best_t = %lf ms (%.2lf Tops/s)\n", 1000*best_t, Gval_per_s(gemm_ops(M, N, K), best_t));

	}
	
  	return 0;
}
