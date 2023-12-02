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

typedef class ATC{
	public:
		long int M, N, K, T; /// The problem dims and tiling size used for 1D/2D Data split to tiles.
		long int Grid_M, Grid_N, Grid_K; /// The resulting grid sizes from the tiling. 
		/// Slowdowns for the selected T that can be used for model adjustment
		double T_aggregate_sl, T_imbalance_sl, T_remainder_sl, T_small_sl, T_sknum_sl, T_big_sl;
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
	
		long int subkernel_num; /// The number of subkernels.
		int* Subkernels_per_unit_num; /// The number of subkernels derived from a unit's score that that unit unit will fire.
		int** Subkernels_per_unit_list; /// The sk_id ids of said sub-kernels, IF they are predefined and not dynamic.
		long long cache_limit; /// The 'cache' size allocation limit for all devices in bytes, IF any.
		Gamalg_p inter_grid; /// The LinkMap representation of the system memory interconnection.
/********************** Initialization/Modification ***************************/
	ATC();	/// Constructor
	~ATC(); /// Destructor
	void reset(); /// Resets controller to default parameters (untuned).
	int diff_intialized_params_ATC(class ATC* other_ATC); /// Rerurns the number of parameters defined in other_ATC that defer from caller.
	void mimic_ATC(class ATC* other_ATC); /// Copy all characteristics of another autotune controller, using its modelers.
	void update_sk_num(long long int subkernel_num_in); /// Updates the autotuner for a given number of subkernels.

	/// 2D block cyclic distribution
	void distribute_subkernels(int D1GridSz, int D2GridSz, int D3GridSz);
/******************************************************************************/
/****************************** Autotuning ************************************/
	double autotune_problem(int A_loc, int B_loc, int C_loc, int D_loc, 
    int M, int N, int K, int elemSize); /// Fire the autotuner for a given problem.
	double optimize_tile(); ///  Predicts the best tile T for a multi-unit problem
	void get_T_slowdowns(double* slowdowns, int candidate_T);
	void set_T_slowdowns(double* slowdowns);
/******************************************************************************/
/**************************** Helper Fuctions *********************************/
	void print(); /// Print the characteristics of the autotune controller to stderr
	const char* print_csv(); /// Print the basic characteristics of the autotune controller to a string in csv-friendly format (X,Y,...)
/******************************************************************************/

}* ATC_p;

#endif
