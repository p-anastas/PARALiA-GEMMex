/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief The header containing the "Subkernel" definition for data scheduling and management in heterogeneous multi-device systems.
///

#ifndef Subkernel_H
#define Subkernel_H

#include<iostream>
#include <string>

//#include "chl_smart_wrappers.hpp"
#include "DataTile.hpp"

int RMConfigResources();
void RMInitResources(int* dev_list, int dev_num);
void RMInitWS(int* dev_list, int dev_num);
void RMFreeResources();
void RMCleanResources();
void RMSyncRecvQueues();

extern CQueue_p recv_queues[64][64];
extern CQueue_p wb_queues[64][64];

extern CQueue_p exec_queue[32][MAX_BACKEND_L_IN];
extern int exec_queue_ctr[32]; 

extern CQueue_p reduce_queue[64][REDUCE_WORKERS_PERDEV];
extern int reduce_queue_ctr[64]; 
extern int reduce_loc; 

#endif
