///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
/// \author Theodoridis Aristomenis (theodoridisaristomenis@gmail.com)
///
/// \brief The PARALiA caching functions.
///

//#include <cassert>
//#include <atomic>

#include "DataCaching.hpp"

Buffer_p current_SAB[64] = {NULL};
int CBlock_ctr[64] = {0};
int DevBuffer_ctr = 0;

#ifdef ENABLE_BUFFER_CONTINUOUS_ALLOC
void* buffer_backup[64] = {NULL};
long long buffer_backup_sz[64] = {0};
#endif

const char* print_state(state in_state){
	switch(in_state){
		case(NATIVE):
			return "NATIVE";
		case(EXCLUSIVE):
			return "EXCLUSIVE";
		case(SHARABLE):
			return "SHARABLE";
		case(AVAILABLE):
			return "AVAILABLE";
		case(INVALID):
			return "INVALID";
		default:
			error("DataCaching -> print_state: Unknown state\n");
	}
	return "";
}

/***************************
 ** Buffer Block Functions **
 ***************************/

BufferBlock::BufferBlock(int block_id, Buffer_p block_parent, long long block_size){
	// Constructor for Blocks.
	// Args:
	// - block_id: Id of the block.
	// - block_parent: Buffer that the block belongs to.
	// - block_size: Size of usable memory in the block.
#ifndef PRODUCTION
	if(block_parent!=NULL && block_id>=0 && block_id<block_parent->BlockNum && block_size>0){
#endif
#ifdef CDEBUG
		fprintf(stderr, "|-----> [dev_id=%d] BufferBlock::BufferBlock(id=%d, Parent_id=%d, Size=%llu)\n", block_parent->dev_id, block_id, block_parent->id, block_size);
#endif
		id = block_id;
		Parent = block_parent;
		Size = block_size;

		Adrs = NULL; // Will be set by cache allocate.
		State = INVALID;
		Available = new Event();
#ifdef CDEBUG
		fprintf(stderr, "<-----| [dev_id=%d] BufferBlock::BufferBlock()\n", Parent->dev_id);
#endif
#ifndef PRODUCTION
	}
	else{
		if(block_parent==NULL)
			error("BufferBlock::BufferBlock(): Constructor called with no buffer to belong.\n");
		else if(block_id<0 && block_id>=block_parent->BlockNum)
			error("[dev_id=%d] BufferBlock::BufferBlock(): Constructor called with invalid id=%d.\n", block_parent->dev_id, block_id);
		else
			error("[dev_id=%d] BufferBlock::BufferBlock(): Constructor called with invalid mem size=%llu.\n", block_parent->dev_id, block_size);
	}
#endif
}

BufferBlock::~BufferBlock(){
	short lvl = 2;
	// Destructor of the block.
#ifdef CDEBUG
	fprintf(stderr, "|-----> [dev_id=%d] BufferBlock::~BufferBlock()\n", Parent->dev_id);
#endif
	if(State != NATIVE){
		if(Adrs) CHLFree(Adrs, Size, Parent->dev_id);
		delete Available;
#ifdef CDEBUG
		fprintf(stderr, "------- [dev_id=%d] BufferBlock::~BufferBlock(): Deleting non-NATIVE block id =%d\n",
			Parent->dev_id, id);
#endif
	}
#ifdef CDEBUG
	else fprintf(stderr, "------- [dev_id=%d] BufferBlock::~BufferBlock(): Refrain from deleting NATIVE block id =%d\n",
		Parent->dev_id, id);
	fprintf(stderr, "<-----| [dev_id=%d] BufferBlock::~BufferBlock()\n", Parent->dev_id);
#endif
}

void BufferBlock::draw_block(){
	// Draws the block for debugging purposes.
	fprintf(stderr, " Block:   \
		\n_________________________________________\
		\n|  Parent dev_id  | %d\
		\n| - - - - - - - - - - - - - - - - - - - -\
		\n|       Id        | %d\
		\n| - - - - - - - - - - - - - - - - - - - -\
		\n|      Size       | %llu\
		\n| - - - - - - - - - - - - - - - - - - - -\
		\n|      State      | %s \
		\n|________________________________________\
		\n", Parent->dev_id, id, Size, print_state(State));
	
}

// Old one. Hasn't changed.
void* CBlock_AVAIL_wrap(void* CBlock_wraped){
	CBlock_wrap_p CBlock_unwraped = (CBlock_wrap_p) CBlock_wraped;
#ifndef PRODUCTION
	if(CBlock_unwraped->CBlock->State == NATIVE) error("[dev_id=%d] CBlock_AVAIL_wrap: block_id=%d is NATIVE\n", 
		CBlock_unwraped->CBlock->Parent->dev_id, CBlock_unwraped->CBlock->id); 
#endif	
	// TODO: This event was never destroyed using cuda because that would be illegal in a host_func...
	// Creating a new event is legal here because lazy events do not call cuda funcs in the constructor
	free(CBlock_unwraped->CBlock->Available);
	CBlock_unwraped->CBlock->Available = new Event();
	CBlock_unwraped->CBlock->State = AVAILABLE;
#ifdef DEBUG
	fprintf(stderr, "|-----> [dev_id=%d] CBlock_AVAIL_wrap: block_id=%d was made AVAILABLE\n", 
		CBlock_unwraped->CBlock->Parent->dev_id, CBlock_unwraped->CBlock->id);
#endif
	//CBlock_unwraped->CBlock->draw_block();
	free(CBlock_unwraped);
	return NULL;
}

void BufferBlock::reset(bool forceReset){
	// Resets block attibutes if it's AVAILABLE to be used again.
#ifdef CDEBUG
	fprintf(stderr, "|-----> [dev_id=%d] BufferBlock::reset(block_id=%d)\n", Parent->dev_id, id);
#endif
	if(State==INVALID || State==AVAILABLE || forceReset){
		Available->reset();
		if(forceReset && State==NATIVE){
			Adrs = NULL;
		}
		State = AVAILABLE;
#ifdef CDEBUG
		if(forceReset)
			fprintf(stderr, "------- [dev_id=%d] BufferBlock::reset(): Block with id=%d forced to be reset.\n", Parent->dev_id, id);
		else
			fprintf(stderr, "------- [dev_id=%d] BufferBlock::reset(): Block with id=%d was reset.\n", Parent->dev_id, id);
#endif
	}
	else error("[dev_id=%d] BufferBlock::reset(): Reset was called on a %s block.\n", Parent->dev_id, print_state(State));
#ifdef CDEBUG
	fprintf(stderr, "<-----| [dev_id=%d] BufferBlock::reset(block_id=%d)\n", Parent->dev_id, id);
#endif
}

void BufferBlock::allocate(){
	// Allocates a buffer block if not already pointing to some memory (not null!)
#ifdef CDEBUG
	fprintf(stderr, "|-----> [dev_id=%d] BufferBlock::allocate(block_id=%d)\n", Parent->dev_id, id);
#endif
	if(Adrs == NULL){
		Adrs = CHLMalloc(Size, Parent->dev_id, 1);
#ifndef CDEBUG
	}
#else
		fprintf(stderr, "------- [dev_id=%d] BufferBlock::allocate(block_id=%d): Allocated Adrs = %p\n", Parent->dev_id, id, Adrs);
	}
	else fprintf(stderr, "------- [dev_id=%d] BufferBlock::allocate(block_id=%d) -> "
		"Supposedly already allocated block...", Parent->dev_id, id);
	fprintf(stderr, "<-----| [dev_id=%d] BufferBlock::allocate(block_id=%d)\n", Parent->dev_id, id);
#endif
}

/*********************
 ** Buffer Functions **
 *********************/

Buffer::Buffer(int dev_id_in, long long block_num, long long block_size){
#ifdef CDEBUG
	fprintf(stderr, "|-----> [dev_id=%d] Buffer::Buffer(block_num = %lld, block_size = %lld)\n", dev_id_in, block_num, block_size);
#endif
	id = DevBuffer_ctr++;
	dev_id = dev_id_in;
	BlockSize = block_size;
	SerialCtr = 0;
	BlockNum = block_num;
	Size = BlockSize*BlockNum;
	Blocks =  (CBlock_p*) malloc (BlockNum * sizeof(CBlock_p));
	for (int idx = 0; idx < BlockNum; idx++) Blocks[idx] = new BufferBlock(idx, this, BlockSize);
	//buffer_backup[dev_id] = NULL;
	//buffer_backup_sz[dev_id] = 0; 
#ifdef CDEBUG
	fprintf(stderr, "<-----| [dev_id=%d] Buffer::Buffer()\n", dev_id_in);
#endif
}

Buffer::~Buffer(){
#ifdef CDEBUG
	fprintf(stderr, "|-----> [dev_id=%d] Buffer::~Buffer()\n", dev_id);
#endif
	DevBuffer_ctr--;
#ifdef ENABLE_BUFFER_CONTINUOUS_ALLOC
#ifndef BUFFER_REUSE_ENABLE
	if(buffer_backup[dev_id]){
		CHLFree(buffer_backup[dev_id], buffer_backup_sz[dev_id], dev_id);
		buffer_backup[dev_id] = NULL;
		buffer_backup_sz[dev_id] = 0;
	}
#endif
	for (int idx = 0; idx < BlockNum; idx++) if(Blocks[idx]!=NULL)
		Blocks[idx]->Adrs = NULL;
#endif
	for (int idx = 0; idx < BlockNum; idx++) delete Blocks[idx];
	free(Blocks);
#ifdef CDEBUG
	fprintf(stderr, "<-----| [dev_id=%d] Buffer::~Buffer()\n", dev_id);
#endif
	return ;
}

void Buffer::reset(bool forceReset){
#ifdef CDEBUG
	fprintf(stderr, "|-----> [dev_id=%d] Buffer::reset()\n", dev_id);
#endif
	for (int idx = 0; idx < BlockNum; idx++) Blocks[idx]->reset(forceReset);
	SerialCtr = 0;
#ifdef CDEBUG
	fprintf(stderr, "<-----| [dev_id=%d] Buffer::reset()\n", dev_id);
#endif
	return;
}

void Buffer::draw_buffer(bool print_blocks){
	fprintf(stderr, " Buffer:\
		\n==========================================\
		\n||      Buffer Id       | %d\
		\n|| - - - - - - - - - - - - - - - - - - - -\
		\n||      Device Id      | %d\
		\n|| - - - - - - - - - - - - - - - - - - - -\
		\n||        Size         | %llu\
		\n|| - - - - - - - - - - - - - - - - - - - -\
		\n||  Number of blocks   | %d\
		\n|| - - - - - - - - - - - - - - - - - - - -\
		\n||   Size of blocks    | %llu\
		\n|| - - - - - - - - - - - - - - - - - - - -", id, dev_id,Size, BlockNum, BlockSize);
	fprintf(stderr, "\n==========================================");

	if(print_blocks){
		fprintf(stderr, "======================================\
			\n|| Start of Blocks in Buffer ||\
			\n==============================\n");
		for(int i=0; i<BlockNum; i++)
			if(Blocks[i]!=NULL)
				Blocks[i]->draw_block();
		fprintf(stderr, "============================\
			\n|| End of Blocks in Buffer ||\
			\n================================================================================\n");
	}
	fprintf(stderr, "\n");
}

#ifdef ENABLE_BUFFER_CONTINUOUS_ALLOC
/// Allocates all buffer blocks in buffer, but in a single continuous piece of memory.

void Buffer::allocate(){
#ifdef CDEBUG
	fprintf(stderr, "|-----> [dev_id=%d] Buffer::allocate-continuous()\n", dev_id);
#endif
	long long total_sz = 0, total_offset = 0;
	for(int i=0; i<BlockNum; i++) if(Blocks[i]!=NULL && Blocks[i]->Adrs==NULL) total_sz+= Blocks[i]->Size;
	if(total_sz){
		if(!buffer_backup[dev_id]){
			buffer_backup[dev_id] = CHLMalloc(total_sz, dev_id, 1);
			buffer_backup_sz[dev_id] = total_sz;
		}
		else if(buffer_backup[dev_id] && total_sz > buffer_backup_sz[dev_id]){
			CHLFree(buffer_backup[dev_id], buffer_backup_sz[dev_id], dev_id);
			buffer_backup[dev_id] = CHLMalloc(total_sz, dev_id, 1);
			buffer_backup_sz[dev_id] = total_sz;	
		}
	}
	for(int i=0; i<BlockNum; i++)
		if(Blocks[i]!=NULL){
			if(Blocks[i]->Adrs==NULL){
				Blocks[i]->Adrs = buffer_backup[dev_id] + total_offset;
				total_offset+=Blocks[i]->Size;
			}
		}
		else error("[dev_id=%d] Buffer::allocate-continuous() -> Blocks[%d] was NULL\n", dev_id, i);
#ifdef CDEBUG
	fprintf(stderr, "<-----| [dev_id=%d] Buffer::allocate-continuous()\n", dev_id);
#endif
}

#else

void Buffer::allocate(){
	// Allocates all buffer blocks in buffer
#ifdef CDEBUG
	fprintf(stderr, "|-----> [dev_id=%d] Buffer::allocate()\n", dev_id);
#endif
	for(int i=0; i<BlockNum; i++)
		if(Blocks[i]!=NULL) Blocks[i]->allocate();
		else error("[dev_id=%d] Buffer::allocate() -> Blocks[%d] was NULL\n", dev_id, i);
#ifdef CDEBUG
	fprintf(stderr, "<-----| [dev_id=%d] Buffer::allocate()\n", dev_id);
#endif
}

#endif

CBlock_p Buffer::assign_Cblock(state start_state){
	// Assigns a block from buffer to be used for memory.
#ifdef CDEBUG
	fprintf(stderr, "|-----> [dev_id=%d] Buffer::assign_Cblock()\n", dev_id);
#endif
#ifndef PRODUCTION
	if(start_state==INVALID)
		error("[dev_id=%d] Buffer::assign_Cblock(): New block called to be initialized as invalid\n", dev_id);
#endif
	CBlock_p result = NULL;
	if (SerialCtr >= BlockNum){
		int remove_block_idx = -42;
		while(remove_block_idx < 0){
			for (int idx = 0; idx < BlockNum; idx++)
			if(Blocks[idx]->State == INVALID || Blocks[idx]->State == AVAILABLE){
				remove_block_idx = idx;
				break;
			}
		}
#ifdef CDEBUG
		fprintf(stderr, "|-----> [dev_id=%d] Buffer::assign_Cblock(): Selecting block %d\n", dev_id, remove_block_idx);
#endif	
		result = Blocks[remove_block_idx];
	}
	else result = Blocks[SerialCtr++];
	result->reset(false);
	result->State = start_state;
#ifdef CDEBUG
	fprintf(stderr, "<-----| [dev_id=%d] Buffer::assign_Cblock()\n", dev_id);
#endif
  return result;
}

/*********************
 ** Other Functions **
 *********************/

int BufferSelectBlockToRemove_naive(Buffer_p buffer){
	if (buffer == NULL) error("BufferSelectBlockToRemove_naive(): Called on empty buffer\n");
#ifdef CDEBUG
	fprintf(stderr, "|-----> [dev_id=%d] BufferSelectBlockToRemove_naive()\n",buffer->dev_id);
#endif
	int result_idx = -1;
	for (int idx = 0; idx < buffer->BlockNum; idx++){ // Iterate through buffer serially.
		state tmp_state = buffer->Blocks[idx]->State; // Update all events etc for idx.
		if(tmp_state == INVALID || tmp_state == AVAILABLE){ // Indx can be removed if there are no pending events.
			result_idx = idx;
			buffer->Blocks[idx]->State = INVALID;
#ifdef CDEBUG
			fprintf(stderr, "------- [dev_id=%d] BufferSelectBlockToRemove_naive(): Found available block. Invalidated.\n",buffer->dev_id);
#endif
			break;
		}
	}
#ifdef CDEBUG
	fprintf(stderr, "<-----| [dev_id=%d] BufferSelectBlockToRemove_naive()\n",buffer->dev_id);
#endif
	return result_idx;
}