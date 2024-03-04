///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief The core header of the autotuner and its components.
///
#ifndef AUTOTUNER_H
#define AUTOTUNER_H

#include <cstdio>

#include "chl_smart_wrappers.hpp"
#include "chl_grid_amalgamation.hpp"

#ifdef PDEBUG
#ifndef SDEBUG
#define SDEBUG
#endif
#endif

enum TileTaskType{
    FETCH = 1,
    COMPUTE = 2,
    WRITEBACK = 3
};

typedef class TileTask{
public:
    TileTaskType type;
    long long DecomIdx, TileIdx, TileIdy;
	int op_id;
    LinkRoute_p predef_route;
}* Ttask_p;

typedef class ATC{
	public:
		long int M, N, K, T; /// The problem dims and tiling size used for 1D/2D Data split to tiles.
		int A_loc, B_loc, C_loc, D_loc; /// The initial locations of the matrices
		long int Grid_M, Grid_N, Grid_K; /// The resulting grid sizes from the tiling.
		/// Slowdowns for the selected T that can be used for model adjustment
		double T_aggregate_sl, T_imbalance_sl, T_remainder_sl, T_small_sl, T_sknum_sl, T_big_sl;
		int* active_memlocs; // The memlocs that are utilized (workers + inputs + outputs) for a problem.
		int active_memloc_num; // The number of memlocs of the problem.
		int active_unit_num; /// The number of units that will be used in the involving operation.
		int* active_unit_id_list;	/// The list of ids of said units.
		double* active_unit_score; /// The 'score' of each said units relative to the total task completion.
		short split_homogeneously = 0; /// A flag that disables workload ratio selection
		double pred_t; /// The predicted seconds the whole operation will require using the above parameters.
		double pred_J; /// The predicted Joules the whole operation will require using the above parameters.
		double power_delay, energy_delay; /// The predicted power and energy delay products using the above parameters.

		double pred_t_pesimistic; /// The predicted seconds the whole operation will require if all overlap fails.
		double pred_J_pesimistic; /// The predicted Joules the whole operation will require if all overlap fails.
		double power_delay_pesimistic, energy_delay_pesimistic; /// The predicted power and energy delay products if all overlap fails.
	
		long int comp_task_num; // The total number of compute tasks created by the distribution.
		long int* comp_task_per_unit_num;  // The number of compute tasks fired per unit.
		// The unit id for each compute task. Two lists for unit -> task and task -> unit fast search.
		int* comp_task_unit_list, **comp_task_per_unit_list;  

		long int task_num; /// The number of tasks to be executed.
		Ttask_p* task_list; 
		// loc_map: A runtime representation of Tile availability in each device: 
		// - not available = -42
		// - available in location = 42 (exhept initial)
		// - initial location = 0
		// - priority target loc = 2, 
		// - other target loc(s) = 1
    	//int loc_map[64]; 
   		//void set_loc_idx(int loc_idx, int val); // Set the loc_idx element of loc_map to val.
    	//void try_set_loc_idx(int loc_idx, int val); // Similar but can only set uninitialized values (-42).
		int*** A_tile_loc_map = NULL; 
		int*** B_tile_loc_map = NULL; 
		int*** C_tile_loc_map = NULL; 

		long long cache_limit; /// The 'cache' size allocation limit for all devices in bytes, IF any.
		Gamalg_p inter_grid; /// The LinkMap representation of the system memory interconnection.
/********************** Initialization/Modification ***************************/
	ATC();	/// Constructor
	~ATC(); /// Destructor
	void reset(); /// Resets controller to default parameters (untuned).
	int diff_intialized_params_ATC(class ATC* other_ATC); /// Rerurns the number of parameters defined in other_ATC that defer from caller.
	void mimic_ATC(class ATC* other_ATC); /// Copy all characteristics of another autotune controller, using its modelers.

/******************************************************************************/
/********************** Tile & device autotuning ******************************/
	double autotune_problem(int A_loc, int B_loc, int C_loc, int D_loc, 
    int M, int N, int K, int elemSize); /// Fire the autotuner for a given problem.
	double optimize_tile(); ///  Predicts the best tile T for a multi-unit problem
	void get_T_slowdowns(double* slowdowns, int candidate_T);
	void set_T_slowdowns(double* slowdowns);
/******************************************************************************/
/********************** Route & distribution autotuning ***********************/
	void update_comp_task_num(long long int task_num_in); /// Updates the autotuner lists for a given number of tasks.
	/// 2D block cyclic distribution is prefered
	void distribute_comp_tasks();
	void initialize_tasks();
	void optimize_tasks();
/******************************************************************************/
/**************************** Helper Fuctions *********************************/
	void print(); /// Print the characteristics of the autotune controller to stderr
	const char* print_csv(); /// Print the basic characteristics of the autotune controller to a string in csv-friendly format (X,Y,...)
/******************************************************************************/

}* ATC_p;

#endif
