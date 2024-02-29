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
		FLOAT
};

enum WR_properties{
		RONLY = 1,
		WONLY = 2,
        WR = 3,
        WR_LAZY = 4,
        W_REDUCE = 5
};

enum TileTaskType{
    FETCH = 1,
    COMPUTE = 2,
    WRITEBACK = 3
};

typedef class TileTask{
    TileTaskType type;
    long long tile_id; 
    LinkRoute_p predef_route;
}* Ttask_p;

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
			int ldim, int inGrid1, int inGrid2, dtype_enum dtype_in, CBlock_p init_loc_block_p);
	~Tile2D(); 	//Destructor
    
    short get_initial_location();
    void reset(void* new_adrr, int new_ldim, CBlock_p new_init_loc_block_p);

    //----------------------------------------------Tile caching------------------------------------------//
    
    CBlock_p StoreBlock[64];// The softcache blocks that store this Tile in each device.
 
    // loc_map: A runtime representation of Tile availability in each device: 
    // - not available = -42
    // - available in location = 42 (exhept initial)
    // - initial location = 0
    // - priority target loc = 2, 
    // - other target loc(s) = 1
    int loc_map[64]; 

    void set_loc_idx(int loc_idx, int val); // Set the loc_idx element of loc_map to val.
    void try_set_loc_idx(int loc_idx, int val); // Similar but can only set uninitialized values (-42).

    LinkRoute_p in_route;
    void fetch(); // Fetch block to a list of locations using a predefined route.

    //--------------------------------------------Tile properties-----------------------------------------//
    WR_properties WRP; // Defines an enum for the type of the tile. 
    WR_properties get_WRP();
    const char* get_WRP_string();
    void set_WRP(WR_properties inprop);

    //--------------------------------------------WTile properties-----------------------------------------//
    // Note: Only relevant for output tiles (e.g. not RONLY)
    int W_init_loc, W_op_dev_id, W_op_queue_ctr, W_op_num, W_op_fired;
    Event_p W_op_complete, W_wb_complete, W_ready;

    Event_p W_op_dependencies[4] = {NULL};
    int W_op_dep_num;
    void** W_op_params;
	const char* W_op_name;
    void run_operation(int W_op_id); 

    LinkRoute_p out_route;
    void writeback(); // Write back block to initial location using a predefined route.

    /// Only applicable for ALGO_WR_LAZY/ALGO_WREDUCE. 
    /// For ALGO_WR_LAZY: C = reduce_mult * C' + C (axpy)
    /// For ALGO_WREDUCE: C = 1.0 * C' + reduce_mult * C (axpby)
    double reduce_mult; 
    CBlock_p backup_C; 

    /// Only for ALGO_WR_LAZY. Fetch C0 to a temp Cblock and perform C = reduce_mult * C' + C
    /// Must be called after W_op_fired = W_op_num and uses the related W_master_backend_ctr queue.
    void WR_lazy_combine(); 
    
    /// Only for ALGO_WREDUCE
    int backup_C_ldim; 
    void WReduce_backup_C(); // Store the pointer of C0 to backup_C. 
    void WReduce_combine(); // After WB, perform C = 1.0 * C' + reduce_mult * C and restore StoreBlock[init].

    //------------------------------------PARALia 2.0 - timed queues and blocks----------------------------//
    void ETA_add_task(long double task_duration, int dev_id);
    void ETA_set(long double new_workload_t, int dev_id);
    long double ETA_get(int dev_id);
    long double ETA_fetch_estimate(int target_id); 
    long double block_ETA[64]; 

}* Tile2D_p;

#endif