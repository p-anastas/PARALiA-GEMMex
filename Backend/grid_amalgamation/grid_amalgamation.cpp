#include <cmath>

#include "smart_wrappers.hpp"
#include "grid_amalgamation.hpp"

Gamalg_p* system_gamalg = NULL;
int system_gamalg_ctr = 0; 

double best_grid_edge_bws[64][64];
int best_grid_edge_active[64][64];
int best_grid_edge_replaced[64][64][2];

/// TODO: This could be made much more IO-efficient, but not a priority currently.
void system_gamalg_init_from_DB(){
    if(system_gamalg || system_gamalg_ctr) 
        warning("system_gamalg_init_from_DB(): system_gamalg = %p and system_gamalg_ctr = %d, overwriting...\n", 
            system_gamalg, system_gamalg_ctr);
    system_gamalg_ctr = 0; 
    if(system_gamalg){ 
        free(system_gamalg); 
        system_gamalg = NULL;
    }
    int explored_cases = std::pow(2, CHL_WORKERS);
    Gamalg_p temp_gamalgs[explored_cases*64] = {NULL};
    for (int case_idx = 1; case_idx < explored_cases; case_idx++){
        int queue_configuration_list[64][2], queue_configuration_num = 0;
		gamalg_backend_get_configs(case_idx, queue_configuration_list, &queue_configuration_num);
        for (int config_idx = 0; config_idx < queue_configuration_num; config_idx++){
            temp_gamalgs[system_gamalg_ctr] = new Grid_amalgamation(case_idx);
            int load_success = temp_gamalgs[system_gamalg_ctr]->load_edges(queue_configuration_list[config_idx][0], 
				queue_configuration_list[config_idx][1]);
            if (load_success){
//#ifdef CLDEBUG
                fprintf(stderr,"system_gamalg_init_from_DB(): Loaded Grid_amalgamation (in_queue_id = %d, out_queue_id = %d) from DB in system_gamalg_ctr = %d\n",
                    queue_configuration_list[config_idx][0], queue_configuration_list[config_idx][1], system_gamalg_ctr);
//#endif
                system_gamalg_ctr++;
            }
            else{ 
//#ifdef CLDEBUG
                fprintf(stderr,"system_gamalg_init_from_DB(): Combination (in_queue_id = %d, out_queue_id = %d) not found in DB\n", 
                    queue_configuration_list[config_idx][0], queue_configuration_list[config_idx][1]);
//#endif
                delete temp_gamalgs[system_gamalg_ctr];
            }
        }
    }
    if (!system_gamalg_ctr) warning("system_gamalg_init_from_DB(): Found 0 matching configurations in files...something might crash\n");
    system_gamalg = (Gamalg_p*) malloc(system_gamalg_ctr*sizeof(Gamalg_p));
    for (int ctr = 0; ctr < system_gamalg_ctr; ctr++) system_gamalg[ctr] = temp_gamalgs[ctr];
}

void gamalg_backend_get_configs(int case_id, int queue_configuration_list[64][2], int* queue_configuration_num){
    *queue_configuration_num = 1;
    int active_unit_id_list[CHL_WORKERS], active_unit_num;
    translate_binary_to_unit_list(case_id, &active_unit_num, active_unit_id_list);
    //printf("active_unit_id_list = %s\n", printlist(active_unit_id_list, active_unit_num));
    // Configuration with all h<->d links used (e.g. bidirectional, two copy engines per device)
    queue_configuration_list[0][0] = case_id;
    queue_configuration_list[0][1] = case_id;
    // Configuration with all h->d used, half h<-d used
    if(0 == active_unit_num%2){
        queue_configuration_list[*queue_configuration_num][0] = case_id;
        queue_configuration_list[*queue_configuration_num][1] = binary_case_id_split(case_id);
        (*queue_configuration_num)++;

    }
}
 
Grid_amalgamation::Grid_amalgamation(int active_nodes_id_in){
    active_nodes_id = active_nodes_id_in;
    for (int d1 = 0; d1 < 64; d1++)
    for (int d2 = 0; d2 < 64; d2++){
        edge_active[d1][d2] = edge_replaced[d1][d2][0] = edge_replaced[d1][d2][1] = -1;
        edge_bw[d1][d2] = -1;
        simu_edge_bw[d1][d2] = -1;
        problem_edge_bw[d1][d2] = -1;
    }
    for (int d1 = 0; d1 < 32; d1++) node_ops[d1] = node_mem_ops[d1] = 0; 
}

Grid_amalgamation::~Grid_amalgamation(){

}

void Grid_amalgamation::print_edge_active(){
    fprintf(stderr,"\n Grid_amalgamation::print_edge_active():\n\n      |");
    for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
        fprintf(stderr, "  %s  |", mem_name(d2));
    fprintf(stderr, "\n      |");
    for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
        fprintf(stderr, "----------");
    fprintf(stderr, "\n");
    for (int d1 = 0; d1 < CHL_MEMLOCS; d1++){
        fprintf(stderr, "%s | ", mem_name(d1));
        for (int d2 = 0; d2 < CHL_MEMLOCS; d2++){
            fprintf(stderr, "%5d   | ", edge_active[d1][d2]);
        }
        fprintf(stderr, "\n");
    }
}

void Grid_amalgamation::print_edge_replaced(){
    fprintf(stderr,"\n Grid_amalgamation::print_edge_replaced():\n\n      |");
    for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
        fprintf(stderr, "  %s  |", mem_name(d2));
    fprintf(stderr, "\n      |");
    for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
        fprintf(stderr, "----------");
    fprintf(stderr, "\n");
    for (int d1 = 0; d1 < CHL_MEMLOCS; d1++){
        fprintf(stderr, "%s | ", mem_name(d1));
        for (int d2 = 0; d2 < CHL_MEMLOCS; d2++){
            fprintf(stderr, "[%2d, %2d]| ", edge_replaced[d1][d2][0], edge_replaced[d1][d2][1]);
        }
        fprintf(stderr, "\n");
    }
}


void Grid_amalgamation::print_edge_bws(){
    fprintf(stderr,"\n Grid_amalgamation::print_edge_bws():\n\n      |");
    for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
        fprintf(stderr, "  %s  |", mem_name(d2));
    fprintf(stderr, "\n      |");
    for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
        fprintf(stderr, "----------");
    fprintf(stderr, "\n");
    for (int d1 = 0; d1 < CHL_MEMLOCS; d1++){
        fprintf(stderr, "%s | ", mem_name(d1));
        for (int d2 = 0; d2 < CHL_MEMLOCS; d2++){
            fprintf(stderr, "%6.2lf  | ", edge_bw[d1][d2]);
        }
        fprintf(stderr, "\n");
    }
}

void Grid_amalgamation::print_simu_edge_bws(){
    fprintf(stderr,"\n Grid_amalgamation::print_simu_edge_bws():\n\n      |");
    for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
        fprintf(stderr, "  %s  |", mem_name(d2));
    fprintf(stderr, "\n      |");
    for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
        fprintf(stderr, "----------");
    fprintf(stderr, "\n");
    for (int d1 = 0; d1 < CHL_MEMLOCS; d1++){
        fprintf(stderr, "%s | ", mem_name(d1));
        for (int d2 = 0; d2 < CHL_MEMLOCS; d2++){
            fprintf(stderr, "%6.2lf  | ", simu_edge_bw[d1][d2]);
        }
        fprintf(stderr, "\n");
    }
}

void Grid_amalgamation::print_edge_load(){
    fprintf(stderr,"\n Grid_amalgamation::print_edge_load(in MB):\n\n      |");
    for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
        fprintf(stderr, "  %s  |", mem_name(d2));
    fprintf(stderr, "\n      |");
    for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
        fprintf(stderr, "----------");
    fprintf(stderr, "\n");
    for (int d1 = 0; d1 < CHL_MEMLOCS; d1++){
        fprintf(stderr, "%s | ", mem_name(d1));
        for (int d2 = 0; d2 < CHL_MEMLOCS; d2++){
            fprintf(stderr, "%7lld | ", edge_load[d1][d2]/(1024*1024));
        }
        fprintf(stderr, "\n");
    }
}

void Grid_amalgamation::print_problem_edge_bws(){
    fprintf(stderr,"\n Grid_amalgamation::print_problem_edge_bws():\n\n      |");
    for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
        fprintf(stderr, "  %s  |", mem_name(d2));
    fprintf(stderr, "\n      |");
    for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
        fprintf(stderr, "----------");
    fprintf(stderr, "\n");
    for (int d1 = 0; d1 < CHL_MEMLOCS; d1++){
        fprintf(stderr, "%s | ", mem_name(d1));
        for (int d2 = 0; d2 < CHL_MEMLOCS; d2++){
            fprintf(stderr, "%6.2lf  | ", problem_edge_bw[d1][d2]);
        }
        fprintf(stderr, "\n");
    }    
}

void Grid_amalgamation::copy(class Grid_amalgamation* gamalg){

}

void Grid_amalgamation::reset_problem_edge_bws(){
      for (int d1 = 0; d1 < CHL_MEMLOCS; d1++) 
        for (int d2 = 0; d2 < CHL_MEMLOCS; d2++) 
            problem_edge_bw[d1][d2] = -1;
}

/*    
void load_results(char* filename, int numDev, int* case_id_list, double case_id_bw_list[][3], int &case_id_lists_len){
	FILE* fp = fopen(filename,"w");
	if (!fp) error("changelink_select_grids: log_results: LogFile failed to open");
	for (int idx = 0; idx < case_id_lists_len; idx++)
		fprintf(fp,"%d, %lf,%lf,%lf\n", case_id_list[idx], case_id_bw_list[idx][0], case_id_bw_list[idx][1], case_id_bw_list[idx][2]);
    fclose(fp);
}

void log_results_bid(char* filename, int numDev, int* case_id_list, int* rev_case_id_list, double case_id_bw_list_bid[][3],  double case_id_bw_list_h2d[][3],  double case_id_bw_list_d2h[][3], int &case_id_lists_len){
	FILE* fp = fopen(filename,"w");
	if (!fp) error("changelink_select_grids: log_results: LogFile failed to open");
	for (int idx = 0; idx < case_id_lists_len; idx++)
		fprintf(fp,"%d,%d, %lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", case_id_list[idx], rev_case_id_list[idx], case_id_bw_list_bid[idx][0], case_id_bw_list_bid[idx][1], case_id_bw_list_bid[idx][2]
		, case_id_bw_list_h2d[idx][0], case_id_bw_list_h2d[idx][1], case_id_bw_list_h2d[idx][2]
		, case_id_bw_list_d2h[idx][0], case_id_bw_list_d2h[idx][1], case_id_bw_list_d2h[idx][2]);
    fclose(fp);
}*/

int Grid_amalgamation::load_edges(int case_id, int rev_case_id){
    char *filename = (char *) malloc(1024 * sizeof(char));
#ifdef CLDEBUG
    fprintf(stderr, "Grid_amalgamation::load_edges(%d,%d)\n", case_id, rev_case_id);
#endif
    int active_unit_num, active_unit_id_list[CHL_WORKERS];
    int rev_active_unit_num, rev_active_unit_id_list[CHL_WORKERS];
    if ((!case_id || !rev_case_id) || !is_subset(rev_case_id, case_id)) return 0; 

    sprintf(filename, "%s/Database/chl_bw_grid_%d_%d.log", DEPLOYDB, case_id, rev_case_id);
    FILE* fp = fopen(filename, "r");
    if(!fp){
        warning("Grid_amalgamation::load_edges(%d,%d): File %s not found\n", case_id, rev_case_id, filename);
        return 0;
    }
    int tmp_worker, temp_memloc;
    
    fscanf(fp, "CHL_WORKERS = %d\nCHL_MEMLOCS = %d\nOne-directional:\n", &tmp_worker, &temp_memloc);
    if (tmp_worker != CHL_WORKERS){
        warning("Grid_amalgamation::load_edges: Loaded different CHL_WORKERS"
            " = %d (instead of %d) from %s...probably wrong system logs used\n", tmp_worker, CHL_WORKERS, filename);
        return 0;
    }
    if (temp_memloc != CHL_MEMLOCS){
        warning("Grid_amalgamation::load_edges: Loaded different CHL_MEMLOCS"
            " = %d (instead of %d) from %s...probably wrong system logs used\n", temp_memloc, CHL_MEMLOCS, filename);
        return 0;
    }

    for(int idx = 0; idx < CHL_MEMLOCS; idx++){
        for(int idx1 = 0; idx1 < CHL_MEMLOCS; idx1++){
            fscanf(fp, "%lf", &(edge_bw[idx][idx1]));
        }
        fscanf(fp, "\n");
    }
    fscanf(fp, "Bidirectional:\n");
    int last_idx_hostsrc = -42, last_idx_hostdest = -42;
    for(int idx = 0; idx < CHL_MEMLOCS; idx++){
        for(int idx1 = 0; idx1 < CHL_MEMLOCS; idx1++){
            fscanf(fp, "%lf", &(simu_edge_bw[idx][idx1]));
            if (simu_edge_bw[idx][idx1] != -1.0){
                edge_active[idx][idx1] = 1;
                if (idx >= CHL_WORKERS && idx1 < CHL_WORKERS) last_idx_hostdest = idx1;
                if (idx1 >= CHL_WORKERS && idx < CHL_WORKERS) last_idx_hostsrc = idx;
            }
            else if(idx != idx1){
                if(idx >= CHL_WORKERS){
                    if (idx1 < CHL_WORKERS){
                        edge_replaced[idx][idx1][0] = idx;
                        edge_replaced[idx][idx1][1] = last_idx_hostdest;
                    }
                }
                else if (idx1 >= CHL_WORKERS){
                    if (idx < CHL_WORKERS){
                        edge_replaced[idx][idx1][0] = last_idx_hostsrc;
                        edge_replaced[idx][idx1][1] = idx1;
                    }
                }
                else error("Grid_amalgamation::load_edges() simu_edge_bw[%d][%d] was -1...incompatible with backend\n", 
                    idx, idx1, simu_edge_bw[idx][idx1]);
            } 
        }
        fscanf(fp, "\n");
    }
    
    fclose(fp);
    free(filename);
    return 1;
}

/// Load the characteristics of the worker 'nodes' (e.g. devices) 
/// Also check file layout in case microbenchmarks messed or some file has been edited by hand incorrectly
void Grid_amalgamation::load_nodes(){
    sprintf(filename, "%s/Database/chl_worker_grid_%d.log", DEPLOYDB, CHL_WORKERS);
    FILE* fp = fopen(filename, "r");
    if(!fp){
        error("Grid_amalgamation::load_nodes(): File %s not found\n", filename);
        return 0;
    }
    for (int it = 0; it < 100; it++) fscanf(fp, "=");
    int tmp_worker, tmp_dtype_lines;

    massert(fscanf(fp, "\nCHL_WORKERS = %d\n\n", &tmp_worker), "Grid_amalgamation::load_nodes(): "
        "%s -> Wrong worker grid file layout at CHL_WORKERS\n", filename);
    massert(tmp_worker == CHL_WORKERS, "Grid_amalgamation::load_nodes: Loaded different CHL_WORKERS"
            " = %d (instead of %d) from %s\n", tmp_worker, CHL_WORKERS, filename);

    massert(fscanf(fp, "WORKER_GOPS: %d\n", &tmp_dtype_lines), "Grid_amalgamation::load_nodes(): "
        "%s -> Wrong worker grid file layout at WORKER_GOPS\n", filename);
    massert(tmp_dtype_lines == DTYPE_NUM, "Grid_amalgamation::load_nodes: Loaded different DTYPE_NUM"
            " = %d (instead of %d) from %s\n", tmp_dtype_lines, DTYPE_NUM, filename);

    for (int dtidx = 0; dtidx < DTYPE_NUM; dtidx++){
        dtype_name[dtidx] = malloc (256*sizeof(char));
        massert(fscanf(fp, "%s:", &(dtype_name[dtidx])), "Grid_amalgamation::load_nodes(): "
        "%s -> Wrong worker grid file layout at dtype naming - dtidx = %d\n", filename, dtidx);
        for (int widx = 0; widx < CHL_WORKERS; widx++)
            massert(fscanf(fp, " %d", &(node_Gops_s[dtidx][widx])), "Grid_amalgamation::load_nodes(): "
            "%s -> Wrong worker grid file layout at node_Gops_s[%d][%d]\n", filename, dtidx, widx);
        fscanf(fp, "\n");
    }

    fscanf(fp, "\nWORKER_MEMOPS:\nGB_s:");
    for (int widx = 0; widx < CHL_WORKERS; widx++){
            massert(fscanf(fp, " %d", &(node_mem_Gb_s[widx])), "Grid_amalgamation::load_nodes(): "
            "%s -> Wrong worker grid file layout at node_mem_Gb_s[%d]\n", filename, widx);
    }

    fscanf(fp, "\n\nWORKER_POWER:\nWATTS:");
    for (int widx = 0; widx < CHL_WORKERS; widx++){
            massert(fscanf(fp, " %d", &(node_watts[widx])), "Grid_amalgamation::load_nodes(): "
            "%s -> Wrong worker grid file layout at node_watts[%d]\n", filename, widx);
    }

    //for (int it = 0; it < 100; it++) fscanf(fp, "=");
}

void Grid_amalgamation::print_nodes(){
    fprintf(stderr,"\n Grid_amalgamation::print_nodes():\n\n               |");
    for (int d2 = 0; d2 < CHL_WORKERS; d2++)
        fprintf(stderr, "  %s   |", mem_name(d2));
    fprintf(stderr, "\n               |");
    for (int d2 = 0; d2 < CHL_WORKERS; d2++)
        fprintf(stderr, "-----------");
    fprintf(stderr, "\n node_Gops_s   | ");
    for (int d2 = 0; d2 < CHL_WORKERS; d2++){
        fprintf(stderr, "%7.0lf  | ", node_Gops_s[d2]);
    }
    fprintf(stderr, "\n node_mem_Gb_s | ");
    for (int d2 = 0; d2 < CHL_WORKERS; d2++){
        fprintf(stderr, "%7.0lf  | ", node_mem_Gb_s[d2]);
    }
    fprintf(stderr, "\n node_watts    | ");
    for (int d2 = 0; d2 < CHL_WORKERS; d2++){
        fprintf(stderr, "%7.0lf  | ", node_watts[d2]);
    }
    fprintf(stderr, "\n");  
}

void Grid_amalgamation::set_edge_load(long long edge_load_in[64][64]){
    for (int idx = 0; idx < CHL_MEMLOCS; idx++) for (int idx1 = 0; idx1 < CHL_MEMLOCS; idx1++) 
        edge_load[idx][idx1] = edge_load_in[idx][idx1];
}

void Grid_amalgamation::update_problem_edges(){
    long long pure_load[CHL_MEMLOCS][CHL_MEMLOCS], simu_load[CHL_MEMLOCS][CHL_MEMLOCS];
    for (int idx = 0; idx < CHL_MEMLOCS; idx++){
        for (int idx1 = 0; idx1 < CHL_MEMLOCS; idx1++){
            pure_load[idx][idx1] = simu_load[idx][idx1] = 0; 
            if (edge_load[idx][idx1]){
                simu_load[idx][idx1] = std::min(edge_load[idx][idx1], edge_load[idx1][idx]);
                pure_load[idx][idx1] = edge_load[idx][idx1] - simu_load[idx][idx1];
                problem_edge_bw[idx][idx1] = 1.0*edge_load[idx][idx1]/(simu_load[idx][idx1]/simu_edge_bw[idx][idx1] + 
                    pure_load[idx][idx1]/edge_bw[idx][idx1]);
            }
            else problem_edge_bw[idx][idx1] = edge_bw[idx][idx1];
        }
    }
}

void Grid_amalgamation::set_node_load(long long node_ops_in[32], long long node_mem_ops_in[32]){
    for (int idx = 0; idx < CHL_WORKERS; idx++){
        node_ops[idx] = node_ops_in[idx];
        node_mem_ops[idx] = node_mem_ops_in[idx];
    }
}

void Grid_amalgamation::print_node_load(){
    fprintf(stderr,"\n Grid_amalgamation::print_node_load():\n\n                     |");
    for (int d2 = 0; d2 < CHL_WORKERS; d2++)
        fprintf(stderr, "  %s   |", mem_name(d2));
    fprintf(stderr, "\n                     |");
    for (int d2 = 0; d2 < CHL_WORKERS; d2++)
        fprintf(stderr, "-----------");
    fprintf(stderr, "\n node_ops (Gops)     | ");
    for (int d2 = 0; d2 < CHL_WORKERS; d2++){
        fprintf(stderr, "%7lld  | ", node_ops[d2]/(1024*1024*1024));
    }
    fprintf(stderr, "\n node_mem_ops (Mops) | ");
    for (int d2 = 0; d2 < CHL_WORKERS; d2++){
        fprintf(stderr, "%7lld  | ", node_mem_ops[d2]/(1024*1024));
    }
    fprintf(stderr, "\n");  
}

double Grid_amalgamation::get_problem_perf_estimation(){
    double exec_t = 0, mem_t = 0, comm_t = 0; 
    long long total_ops = 0, total_mem_ops = 0, total_comm_bytes = 0; 
    int active_unit_num, active_unit_id_list[CHL_WORKERS];
    translate_binary_to_unit_list(active_nodes_id, &active_unit_num, active_unit_id_list);
    for (int idx = 0; idx < CHL_WORKERS; idx++){
        if(is_in_list(idx, active_unit_id_list, active_unit_num)){
            exec_t = std::max(node_ops[idx]/(node_Gops_s[idx]*1e9), exec_t);
            total_ops+= node_ops[idx];
            mem_t = std::max(node_mem_ops[idx]/(node_mem_Gb_s[idx]*1e9), mem_t);
            total_mem_ops+= node_mem_ops[idx]; 
        }
        for (int idx1 = 0; idx1 < CHL_WORKERS; idx1++){
            if (idx == idx1) continue;
            if(simu_edge_bw[idx][idx1]!= -1 && edge_load[idx][idx1]){
                total_comm_bytes+= edge_load[idx][idx1];
                comm_t = std::max(edge_load[idx][idx1]/(problem_edge_bw[idx][idx1]*1e9), comm_t);
            }
        }
        double comm_temp = 0, comp_temp_rev = 0;
        for (int idx1 = CHL_WORKERS; idx1 < CHL_MEMLOCS; idx1++){
            if(simu_edge_bw[idx][idx1]!= -1 && edge_load[idx][idx1]){
                total_comm_bytes+= edge_load[idx][idx1];
                comm_temp += edge_load[idx][idx1]/(problem_edge_bw[idx][idx1]*1e9);
            }
            if(simu_edge_bw[idx1][idx]!= -1 && edge_load[idx1][idx]){
                total_comm_bytes+= edge_load[idx1][idx];
                comp_temp_rev += edge_load[idx1][idx]/(problem_edge_bw[idx1][idx]*1e9);
            }
        }
        comm_t = std::max(std::max(comm_temp, comp_temp_rev), comm_t);
    }

    double total_t = std::max(std::max(exec_t, mem_t), comm_t);
#ifdef CLDEBUG
    fprintf(stderr,"\nGrid_amalgamation::get_problem_perf_estimation(): Node num = %d, total_t = %lf ms (%lf Top/s), \n-> exec_t = %lf ms (%lf Top/s)\n"
    "-> mem_t = %lf ms (%lf Gmem_op/s)\n-> comm_t = %lf ms (%lf Gb/s)\n",
        active_unit_num, total_t*1000, Gval_per_s(total_ops, total_t)/1024, exec_t*1000, Gval_per_s(total_ops, exec_t), mem_t*1000, Gval_per_s(total_mem_ops, mem_t), 
        comm_t*1000, Gval_per_s(total_comm_bytes, comm_t));
#endif
    return total_t;
}