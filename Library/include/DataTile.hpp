///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief Generalised DataTiles for 1D and 2D tile definition, transfer and sharing. 
///

#ifndef DATATILE_H
#define DATATILE_H

#include "DataCaching.hpp"

enum dtype_enum{
		DOUBLE,
		FLOAT,
        HALF
};

enum WR_properties{
		RONLY = 1,
		WONLY = 2,
        WR = 3,
        WR_LAZY = 4,
        W_REDUCE = 5
};

typedef class Tile2D{
public:
    //----------------------------------------------General class-----------------------------------------//

    dtype_enum dtype; // An enum of the datatype (used to select the operation)
    int get_dtype_size();

    int id, GridId1, GridId2; // The tile location in a 1D/2D grid (1D = order of initialization). 
    int dim1, dim2, ldim[64]; // The tile dimensions and its decomposer's leading dimension.
    int size();

	// Constructor : Initializes a tile with certain dimensions on a grid (used by the decomposer)
	Tile2D(void* tile_addr, int T1tmp, int T2tmp,
			int ldim, int inGrid1, int inGrid2, dtype_enum dtype_in, int init_loc, Buffer_p* init_loc_cache_p);
	~Tile2D(); 	//Destructor
    
    short get_initial_location();
    void reset(void* new_adrr, int new_ldim, Buffer_p* init_loc_cache_p);

    //----------------------------------------------Tile caching------------------------------------------//
    
    CBlock_p StoreBlock[64];// The softcache blocks that store this Tile in each device.
    int Block_reuses[64];// The remaining comp tasks that need to use this block in each location

    void fetch(LinkRoute_p in_route); // Fetch block to a list of locations using a predefined route.

    //--------------------------------------------Tile properties-----------------------------------------//
    WR_properties WRP; // Defines an enum for the type of the tile. 
    WR_properties get_WRP();
    const char* get_WRP_string();
    void set_WRP(WR_properties inprop);

    //--------------------------------------------WTile properties-----------------------------------------//
    // Note: Only relevant for output tiles (e.g. not RONLY)
    int W_init_loc, W_op_dev_id, W_op_queue_ctr, W_op_num, W_op_fired;
    Event_p W_op_complete, W_wb_complete, W_ready;

    void** W_op_params;
	const char* W_op_name;
    void run_operation(int W_op_id, LinkRoute_p lazy_route); 

    void writeback(LinkRoute_p out_route); // Write back block to initial location using a predefined route.

    /// Only applicable for ALGO_WR_LAZY/ALGO_WREDUCE. 
    /// For ALGO_WR_LAZY: C = reduce_mult * C' + C (axpy)
    /// For ALGO_WREDUCE: C = 1.0 * C' + reduce_mult * C (axpby)
    void* reduce_mult; 
    CBlock_p backup_C; 

    /// Only for ALGO_WR_LAZY. Fetch C0 to a temp Cblock and perform C = reduce_mult * C' + C
    /// Must be called after W_op_fired = W_op_num and uses the related W_master_backend_ctr queue.
    void WR_lazy_combine(LinkRoute_p lazy_route); 
    
    /// Only for ALGO_WREDUCE
    int backup_C_ldim; 
    void WReduce_backup_C(); // Store the pointer of C0 to backup_C. 
    void WReduce_combine(); // After WB, perform C = 1.0 * C' + reduce_mult * C and restore StoreBlock[init].

    //------------------------------------PARALia 2.0 - timed queues and blocks----------------------------//
    //void ETA_add_task(long double task_duration, int dev_id);
    //void ETA_set(long double new_workload_t, int dev_id);
    //long double ETA_get(int dev_id);
    //long double ETA_fetch_estimate(int target_id); 
    //long double block_ETA[64]; 

}* Tile2D_p;

void set_val(dtype_enum dtype, void** wrap_ptr, double value);
/// Defines if the task scheduler has to make software-buffer blocks AVAILABLE
/// when their corresponding DataTiles are no longer needed in some location
extern int conserve_memory_curr; 

#endif