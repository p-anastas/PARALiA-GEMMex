/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief The header containing the "Subkernel" definition for data scheduling and management in heterogeneous multi-device systems.
///

#ifndef Subkernel_H
#define Subkernel_H

#include<iostream>
#include <string>

#include "DataTile.hpp"

class Subkernel
{
	private:
	public:
		int id, iloc1, iloc2, iloc3;
		int run_dev_id, launch_order, launched;
		int TileNum;
		Tile2D_p* TileList;
		LinkRoute_p* predef_route;

#ifdef STEST
		double req_in_ts = 0, req_out_ts = 0;
		double reqT_fire_ts[3] = {0}, reqT_start_ts[3] = {0}, reqT_end_ts[3] = {0};
		double op_fire_ts = 0, op_start_ts = 0, op_end_ts = 0;
		double wbT_fire_ts[3] = {0}, wbT_start_ts[3] = {0}, wbT_end_ts[3] = {0};
		long long bytes_in[3] = {0}, bytes_out[3] = {0}, flops = 0;
		int dev_in_from[3], dev_in_to[3], dev_out_from[3], dev_out_to[3];
#endif
		//Event* operation_complete;
		void* operation_params;
		const char* op_name;

		/// Constructors
		Subkernel(int TileNum, const char* name);

		/// Destructors
		~Subkernel();

		void sync_request_data();

		/// Functions
		void prepare_launch(int dev_id);
		void request_data();
		void run_operation();
		int check_ready(); 
		void run_ready_operation();

		void reset();

		/*****************************************************/
		/// PARALia 2.0 - timed queues and blocks
		long double fetch_ETA, run_op_est_t;
		//long double run_op_estimate(ATC_p autotuner); 
};

//Subkernel** CoCoAsignTilesToSubkernelsGemm(Decom2D<double>* A_asset, Decom2D<double>* B_asset,
//Decom2D<double>* C_asset, int T, int* kernelNum);

int SKConfigResources();
void SKInitResources(int* dev_list, int dev_num);
void SKInitWS(int* dev_list, int dev_num);
void SKFreeResources();
void SKCleanResources();

void sync_recv_queues();

int SubkernelPrefetchCheapRONLYTiles(int numTiles, int dev_id, Subkernel** Subkernel_list, long Subkernel_list_len);

Subkernel* SubkernelSelect(int dev_id, Subkernel** Subkernel_list, long Subkernel_list_len);

void sync_request_paired(int dev_id);

#ifdef STEST
void STEST_print_SK(kernel_pthread_wrap_p* thread_dev_data_list, double routine_entry_ts, int dev_num);
#endif

extern CQueue_p recv_queues[64][64];
extern CQueue_p wb_queues[64][64];

extern CQueue_p exec_queue[32][MAX_BACKEND_L];
extern int exec_queue_ctr[32]; 

extern CQueue_p reduce_queue[64][REDUCE_WORKERS_PERDEV];
extern int reduce_queue_ctr[64]; 
extern int reduce_loc; 

#endif
