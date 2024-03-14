///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief The headers for functions for general use throught CHL
///

#ifndef LINKMAP_H
#define LINKMAP_H

#include<iostream>
#include <string>
#include <cstdlib>
#include <cstring>
#include <atomic>

#include <stdio.h>
#include <stdarg.h>

/*****************************************************/
/// LinkRoute stuff

// Struct for multi-hop optimized transfers
typedef class LinkRoute{
public:
	int hop_num;
	int hop_uid_list[64];
	void* hop_buf_list[64];
	int hop_ldim_list[64];
	int starting_hop = 0;
	int streaming_workers = 8;

	CQueue_p hop_cqueue_list[64-1];
	Event_p hop_event_list[64-1];
	
	/// Target: 0/42 -> 2 [+1] wrapper
	long double optimize(int* loc_map, long int size, int update_flag);
	// Specific implementations
	long double optimize_p2p_init(int* loc_map, long int size);
	long double optimize_p2p_serial(int* loc_map, long int size);
	long double optimize_p2p_distance(int* loc_map, long int size);
	long double optimize_chain_serial(int* loc_map, long int size);
	long double optimize_chain_random(int* loc_map, long int size);
	long double optimize_chain_time(int* loc_map, long int size);
	long double optimize_chain_ETA(int* loc_map, long int size, int update_flag);

	/// Target: 42 -> 0 wrapper
	long double optimize_reverse(int* loc_map, long int size);
	// Specific implementations
	long double optimize_reverse_p2p_init(int* loc_map, long int size);
	
	/// Target: SRC -> ? -> dest 
	//long double optimize_hop_route(void* transfer_tile_wrapped, int update_ETA_flag, int dest_loc, int src_loc); 
	void print();
}* LinkRoute_p;

// A memcpy implementation using multiple units as intermendiate hops for a better transfer bandwidth
void FasTCHLMemcpy2DAsync(LinkRoute_p roadMap, long int rows, long int cols, short elemSize);
// Print and log bandwidths and links used with FasTCHLMemcpy2DAsync. Unusable with TTEST flag
void HopMemcpyPrint();

/*****************************************************/
/// LinkMap stuff

#define MAX_ALLOWED_HOPS 1
#define MAX_HOP_ROUTES 1
#define HOP_PENALTY 0.2

typedef class Grid_amalgamation{
	public:

		/********************** System-wide ***************************/
		// A corresponding node 
		int active_nodes_id = -1;

		/// Input function that loads edge_active, edge_bw and simu_edge_bw for a certain case_id, rev_case_id combination.
		int load_edges(int case_id, int rev_case_id);

		// Define which connections are active in the grid. 
		int edge_active[64][64] = {{0}};
		int edge_replaced[64][64][2] = {{0}};
		/// Empirically obtained bandwidths for all edges.
		double edge_bw[64][64] = {{0}};

		// Empirically obtained bandwidths for all edges if used simultaneously
		double simu_edge_bw[64][64] = {{0}};

		void print_edge_active();
		void print_edge_replaced();
		void print_edge_bws();
		void print_simu_edge_bws();

		double node_Gops_s[32], node_mem_Gb_s[32], node_watts[32];
		void load_nodes();

		void print_nodes();

		/******************************************************************************/
		/********************** Problem-specific ***************************/

		/// The estimated load of each edge in bytes.
		long long edge_load[64][64] = {{0}};
		/// Loads the estimated load of each edge in bytes. 
		void set_edge_load(long long edge_load_in[64][64]);

		// Estimated bandwidths for a given problem based on its edge_bw and simu_edge_bw being under edge_load.
		double problem_edge_bw[64][64] = {{0}};
		void update_problem_edges();

		void print_edge_load();
		void print_problem_edge_bws();

		long long node_ops[32], node_mem_ops[32];
		void set_node_load(long long node_ops_in[32], long long node_mem_ops_in[32]);

		void print_node_load();

		double get_problem_perf_estimation();

		/******************************************************************************/
		/********************** Initialization/Modification ***************************/
		Grid_amalgamation(int active_nodes_id_in);
		~Grid_amalgamation();
		/******************************************************************************/
		/**************************** Helper Fuctions *********************************/
		void copy(class Grid_amalgamation* gamalg);
		void reset_problem_edge_bws(); // Resets problem_edge_bw

		/******************************************************************************/
		/************************ Class main Functions ********************************/
}* Gamalg_p;

extern Gamalg_p* system_gamalg;
extern int system_gamalg_ctr; 
void system_gamalg_init_from_DB();

void gemm_translate_problem_comm(long long edge_load[64][64], int A_loc, int B_loc, int C_loc, int D_loc, 
    int M, int N, int K, int elemSize, int active_unit_num, int* active_unit_id_list, double* active_unit_score);

void gemm_translate_problem_ops(long long node_ops[32], long long node_mem_ops[32], 
	int M, int N, int K, int active_unit_num, int* active_unit_id_list, double* active_unit_score);

extern double best_grid_edge_bws[64][64];
extern int best_grid_edge_active[64][64];
extern int best_grid_edge_replaced[64][64][2];

#endif