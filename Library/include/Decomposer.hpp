/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief The header containing the "Asset" definition for data scheduling and management in heterogeneous multi-device systems.
///

#ifndef DECOM_H
#define DECOM_H

#include<iostream>
#include <string>
#include <mutex> // std::mutex

//#include "linkmap.hpp"
#include <atomic>

#include "DataTile.hpp"

typedef class Decomposer
{
	// Variables
	private:
	public:
	dtype_enum dtype;
	void *adrs;
	int id;
	char transpose;
	int GridSz1, GridSz2;
	int loc;
	int dim1, dim2;
	short pin_internally;
	Tile2D_p *Tile_map;

	// General Functions
	virtual long get_chunk_size();
	virtual	void set_chunk_size(long value); 
	virtual long get_mem_size(); 
	virtual void InitTileMap(int T1, int T2, Buffer_p* init_loc_cache_p, WR_properties prop);
	virtual void Reset(void* new_adrs, int T1, int T2, long new_chunk_size, Buffer_p* init_loc_cache_p);

	void WBTileMap();
	void SyncTileMap();
	void DestroyTileMap();
	int dtypesize() {
			if (dtype == DOUBLE) return sizeof(double);
			else if (dtype == FLOAT) return sizeof(float);
			else if (dtype == HALF) return sizeof(__half);
			else error("dtypesize: Unknown type"); return 0;}
	int size() { return dtypesize()*dim1*dim2; }
	void DrawTileMap();

	// Backend Functions
	void prepareAsync();
	void resetProperties();
}* Decomposer_p;


typedef class Decom2D : public Decomposer {
	public:
	int ldim;

	long get_chunk_size(); 
	void set_chunk_size(long value); 
	long get_mem_size(); 
	// Constructor, sets dim1, dim2, ldim, adrs and derives loc from get_loc(adr)
	Decom2D(void* adrr, int in_dim1, int in_dim2, int in_ldim, char transpose, dtype_enum dtype_in);
	void InitTileMap(int T1, int T2, Buffer_p* init_loc_cache_p, WR_properties prop);
	void Reset(void* new_adrs, int T1, int T2, long new_ldim, Buffer_p* init_loc_cache_p);
	void MatrixReset(void* new_adrs, long new_ldim);
	Tile2D_p getTile(int iloc1, int iloc2);


}* Decom2D_p;

#endif
