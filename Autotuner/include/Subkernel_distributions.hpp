///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief The subkernel distributions implemented for PARALiA
///
#ifndef SUBKERNEL_DIST_H
#define SUBKERNEL_DIST_H

#include "Autotuner.hpp"

void CoCoDistributeTasksRoundRobin(ATC_p autotune_controller);
void CoCoDistributeTasksNaive(ATC_p autotune_controller);
void CoCoDistributeTasksRoundRobinChunk(ATC_p autotune_controller,  int Chunk_size);
void CoCoDistributeTasksRoundRobinChunkReverse(ATC_p autotune_controller,  int Chunk_size);
void CoCoDistributeTasks2DBlockCyclic(ATC_p autotune_controller, int D1GridSz, int D2GridSz, int D3GridSz);
#endif
