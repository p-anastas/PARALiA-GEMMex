///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
/// \author Theodoridis Aristomenis (theodoridisaristomenis@gmail.com)
///
/// \brief The header containing the caching functions for data scheduling and management in heterogeneous multi-device systems.
///

#ifndef DATACACHNING_H
#define DATACACHNING_H

#include<iostream>
#include <string>
#include <mutex>
#include <atomic>

#include "chl_smart_wrappers.hpp"
#include "chl_grid_amalgamation.hpp"

enum state{
	NATIVE = 0, /// Buffer Block is native in memory, Should never be scheduled out or have its Adrs freed
	EXCLUSIVE = 1,  /// is being modified locally.
	SHARABLE = 2,  /// is available for sharing only but cannot be modified.
	AVAILABLE = 3, /// Buffer Block is available for usage
	INVALID = 4 /// Buffer block is not ready for usage. Might not be needed in simplified caching
};

const char* print_state(state in_state);

typedef class Buffer* Buffer_p;
typedef class BufferBlock* CBlock_p;

typedef struct Node_LL* Node_LL_p;
typedef class LinkedList* LinkedList_p;

// A class for each Buffer block.
typedef class BufferBlock{
	private:
	public:
		int id; // A unique per DevBuffer id for each block
		Buffer_p Parent;  // Is this needed?
		long long Size; // Included here but should be available at parent DevBuffer (?)

		void* Adrs;
		state State; // The current state of the block.
		Event_p Available;
		
		//Constructor
		BufferBlock(int id, Buffer_p Parent, long long Size);
		//Destructor
		~BufferBlock();

		// Functions
		void draw_block();
		void allocate();
		void reset(bool forceReset=false);  // Cleans a block to be given to someone else

}* CBlock_p;

/// Device-wise software Buffer class declaration
typedef class Buffer{
	private:
	public:
		int id; // A unique id per Buffer
		int dev_id; /// Pressumably this should be sufficient for current use cases instead of id, since we use only 1 Buffer/dev
		long long Size; // The sum of a Buffer's CBlock_sizes.
		void* cont_buf_head; /// Used only if ENABLE_BUFFER_CONTINUOUS_ALLOC
		long long cont_buf_head_sz;
		
		int SerialCtr; // Number of blocks currently in buffer.
		int BlockNum; // Number of Blocks the buffer holds
		long long BlockSize; // Size allocated for each block, including padding. 
		CBlock_p* Blocks;

		//Constructor
		Buffer(int dev_id, long long block_num, long long block_size);
		//Destructor
		~Buffer();

		// Functions
		void draw_buffer(bool print_blocks=true);
		void allocate();
		void reset(bool forceReset=false);
		CBlock_p assign_Cblock(state start_state=AVAILABLE);

}* Buffer_p;

typedef struct CBlock_wrap{
	CBlock_p CBlock;
}* CBlock_wrap_p;

void* CBlock_AVAIL_wrap(void* CBlock_wraped);

extern Buffer_p current_SAB[64];
#endif
