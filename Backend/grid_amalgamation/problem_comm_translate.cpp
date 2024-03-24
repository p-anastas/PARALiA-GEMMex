

#include "smart_wrappers.hpp"
#include "grid_amalgamation.hpp"

#include <cmath>

void gemm_translate_problem_comm(long long edge_load[64][64], int A_loc, int B_loc, int C_loc, int D_loc, 
    int M, int N, int K, int elemSize, int active_unit_num, int* active_unit_id_list, double* active_unit_score){
#ifdef PDEBUG
    fprintf(stderr, "gemm_translate_problem_comm():\n-> M = %d, N = %d, K = %d, elemSize =%d\n"
        "-> A_loc = %d, B_loc = %d, C_loc = %d, D_loc = %d\n-> active_unit_id_list = %s, active_unit_score = %s\n",
        M, N, K, elemSize, A_loc, B_loc, C_loc, D_loc, 
        printlist<int>(active_unit_id_list, active_unit_num), printlist<double>(active_unit_score, active_unit_num));
#endif
    int case_id = translate_unit_list_to_binary(active_unit_id_list, active_unit_num);
    for (int idx = 0; idx < CHL_MEMLOCS; idx++) for (int idx1 = 0; idx1 < CHL_MEMLOCS; idx1++) edge_load[idx][idx1] = 0;
    for (int idx = 0; idx < active_unit_num; idx++){
        int dev_id = active_unit_id_list[idx];

        if(A_loc >= 0) edge_load[dev_id][A_loc]+= (long long) (active_unit_score[idx]*M*K*elemSize);
        else error("gemm_translate_problem_comm: A_loc = %d not supported\n", A_loc);
        //else if (A_loc == -1) for(int idx2 = CHL_WORKERS; idx2 < CHL_MEMLOCS - 1; idx2++) edge_load[dev_id][idx2] +=  
        //    (long long) (active_unit_score[idx]*M*K*elemSize/(CHL_MEMLOCS - CHL_WORKERS - 1));
        //else if(A_loc == -2) edge_load[dev_id][CHL_WORKER_CLOSE_TO_MEMLOC[dev_id]] +=  
        //    (long long) (active_unit_score[idx]*M*K*elemSize);

        if(B_loc >= 0) edge_load[dev_id][B_loc]+= (long long) (active_unit_score[idx]*N*K*elemSize);
        else error("gemm_translate_problem_comm: B_loc = %d not supported\n", B_loc);
        //else if (B_loc == -1) for(int idx2 = CHL_WORKERS; idx2 < CHL_MEMLOCS - 1; idx2++) edge_load[dev_id][idx2] +=  
        //    (long long) (active_unit_score[idx]*N*K*elemSize/(CHL_MEMLOCS - CHL_WORKERS - 1));
        //else if(B_loc == -2) edge_load[dev_id][CHL_WORKER_CLOSE_TO_MEMLOC[dev_id]] +=  
        //    (long long) (active_unit_score[idx]*N*K*elemSize);

        if(C_loc >= 0) edge_load[dev_id][C_loc]+= (long long) (active_unit_score[idx]*M*N*elemSize);
        else error("gemm_translate_problem_comm: C_loc = %d not supported\n", C_loc);
        //else if (C_loc == -1) for(int idx2 = CHL_WORKERS; idx2 < CHL_MEMLOCS - 1; idx2++) edge_load[dev_id][idx2] +=  
        //    (long long) (active_unit_score[idx]*M*N*elemSize/(CHL_MEMLOCS - CHL_WORKERS - 1));
        //else if(C_loc == -2) edge_load[dev_id][CHL_WORKER_CLOSE_TO_MEMLOC[dev_id]] +=  
        //    (long long) (active_unit_score[idx]*M*N*elemSize);

        if(D_loc >= 0) edge_load[D_loc][dev_id]+= (long long) (active_unit_score[idx]*M*N*elemSize);
        else error("gemm_translate_problem_comm: D_loc = %d not supported\n", D_loc);
        //else if (D_loc == -1) for(int idx2 = CHL_WORKERS; idx2 < CHL_MEMLOCS - 1; idx2++) edge_load[idx2][dev_id] +=  
        //    (long long) (active_unit_score[idx]*M*N*elemSize/(CHL_MEMLOCS - CHL_WORKERS - 1));
        //else if(D_loc == -2) edge_load[CHL_WORKER_CLOSE_TO_MEMLOC[dev_id]][dev_id] +=  
        //   (long long) (active_unit_score[idx]*M*N*elemSize);
    }

    /// Calculate extra transfers created from internal dims due to multi-unit spliting.
	/// Algorithm may vary for other BLAS3, but not at that bridge yet.
	/// The assumtion for extra transfers is made based on the 2D cyclic distribution,
	/// but the estimation is also useful for other distributions as a best case scenario (worse distributions -> more extra transfers).
	int D1_parts, D2_parts;
    DECOM_2D(active_unit_num, &D1_parts, &D2_parts);
    
    for(int unit_idx = 0; unit_idx < active_unit_num; unit_idx++){
        int dev_id = active_unit_id_list[unit_idx];
        int dev_decom_row = unit_idx/D2_parts, dev_decom_col = unit_idx%D2_parts;
        int list_row_bros[active_unit_num], list_row_bro_ctr = 0, list_col_bros[active_unit_num], list_col_bro_ctr = 0;
        for(int unit_idy = 0; unit_idy < active_unit_num; unit_idy++) if(unit_idy != unit_idx){
            int dev_id_bro = active_unit_id_list[unit_idy];
            int dev_decom_row_bro = unit_idy/D2_parts, dev_decom_col_bro = unit_idy%D2_parts;
            if (dev_decom_row == dev_decom_row_bro)	list_row_bros[list_row_bro_ctr++] = dev_id_bro;	
            if (dev_decom_col == dev_decom_col_bro)	list_col_bros[list_col_bro_ctr++] = dev_id_bro;	
        }
        for (int idx = 0; idx < list_row_bro_ctr; idx ++){
            edge_load[dev_id][list_row_bros[idx]] += (long long) M*K*elemSize*active_unit_score[unit_idx];
        }
        for (int idx = 0; idx < list_col_bro_ctr; idx ++){
            edge_load[dev_id][list_col_bros[idx]] += (long long) N*K*elemSize*active_unit_score[unit_idx];
        }
    }
}

void gemm_translate_problem_ops(long long node_ops[32],
	int M, int N, int K, int active_unit_num, int* active_unit_id_list, double* active_unit_score){
#ifdef PDEBUG
    fprintf(stderr, "gemm_translate_problem_comm():\n-> M = %d, N = %d, K = %d\n"
        "-> active_unit_id_list = %s, active_unit_score = %s\n",
        M, N, K, printlist<int>(active_unit_id_list, active_unit_num), 
        printlist<double>(active_unit_score, active_unit_num));
#endif
    long long dgemm_ops = gemm_ops(M, N, K);
    long long dgemm_mem_ops = gemm_mem_ops(M, N, K);
    for (int idx = 0; idx < CHL_WORKERS; idx++) node_ops[idx] = 0;
    for (int idx = 0; idx < active_unit_num; idx++){
        int dev_id = active_unit_id_list[idx];
        node_ops[dev_id] = dgemm_ops/active_unit_num;
    }
}
