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

typedef class DataTile{
public:
    dtype_enum dtype;
    int id, GridId1, GridId2;
    int dim1, dim2;
    WR_properties WRP;
    int W_master = -42;
    int W_pending = -42;
    int fired_times = 0; 
    Event_p W_complete = NULL, W_reduce = NULL; 
    int W_master_backend_ctr = -42;
    double reduce_mult; 
    // loc_map values mean: 
    // - not available = -42
    // - available in location = 42 (exhept initial)
    // - initial location = 0
    // - priority target loc = 2, 
    // - other target loc(s) = 1
    int loc_map[64]; 
    CBlock_p StoreBlock[64];

    // General Functions
    int get_dtype_size();
    short get_initial_location();
    WR_properties get_WRP();
    const char* get_WRP_string();
    int size();
    virtual long get_chunk_size(int loc_idx);
    virtual void set_chunk_size(int loc_idx, long value);

    void set_loc_idx(int loc_idx, int val);
    /// Same functionality with the exheption that can only set uninitialized values. 
    void try_set_loc_idx(int loc_idx, int val);
    void set_WRP(WR_properties inprop);

    LinkRoute_p fetch(CBlock_p target_block, int priority_loc_id, LinkRoute_p in_route);
    LinkRoute_p writeback(CBlock_p WB_block, LinkRoute_p in_route);
    void operations_complete(CQueue_p assigned_exec_queue, 
        LinkRoute_p* in_route_p, LinkRoute_p* out_route_p);


    void reset(void* new_adrr, int new_ldim, CBlock_p new_init_loc_block_p);
    
    /*****************************************************/
    /// PARALia 2.0 - timed queues and blocks
    void ETA_add_task(long double task_duration, int dev_id);
    void ETA_set(long double new_workload_t, int dev_id);
    long double ETA_get(int dev_id);
    long double ETA_fetch_estimate(int target_id); 
    long double block_ETA[64]; 

}* DataTile_p;

class Tile1D : public DataTile {
public:
    int inc[64];

    // Constructor
    Tile1D(void* tile_addr, int T1tmp,
        int inc, int inGrid1, dtype_enum dtype_in, CBlock_p init_loc_block_p);
    //Destructor
    ~Tile1D();

    long get_chunk_size(int loc_idx); 
    void set_chunk_size(int loc_idx, long value);

};

class Tile2D : public DataTile {
public:    
    int ldim[64];
	// Constructor
	Tile2D(void* tile_addr, int T1tmp, int T2tmp,
			int ldim, int inGrid1, int inGrid2, dtype_enum dtype_in, CBlock_p init_loc_block_p);
	//Destructor
	~Tile2D();

    long get_chunk_size(int loc_idx);
    void set_chunk_size(int loc_idx, long value);
};

#endif