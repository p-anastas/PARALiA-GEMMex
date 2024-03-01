///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief The subkernel distributions implemented for PARALiA
///
#ifndef SUBKERNEL_DIST_H
#define SUBKERNEL_DIST_H

#include "Autotuner.hpp"

void DistributeCompTasksRoundRobin(ATC_p autotune_controller);
void DistributeCompTasksNaive(ATC_p autotune_controller);
void DistributeCompTasksRoundRobinChunk(ATC_p autotune_controller,  int Chunk_size);
void DistributeCompTasksRoundRobinChunkReverse(ATC_p autotune_controller,  int Chunk_size);
void DistributeCompTasks2DBlockCyclic(ATC_p autotune_controller, int D1GridSz, int D2GridSz, int D3GridSz);
#endif
