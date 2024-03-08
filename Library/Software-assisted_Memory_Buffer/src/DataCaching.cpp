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
		if(Adrs) CHLFree(Adrs, Parent->dev_id, Size);
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
		\n|       Id        | %d\
		\n| - - - - - - - - - - - - - - - - - - - -\
		\n|      Size       | %llu\
		\n| - - - - - - - - - - - - - - - - - - - -\
		\n|      State      | %s \
		\n|________________________________________\
		\n", id, Size, print_state(State));
	
}

// Old one. Hasn't changed.
void* CBlock_INV_wrap(void* CBlock_wraped){
	CBlock_wrap_p CBlock_unwraped = (CBlock_wrap_p) CBlock_wraped;
	CBlock_unwraped->CBlock->Available->soft_reset();
	CBlock_unwraped->CBlock->set_state(INVALID, CBlock_unwraped->lockfree);
	free(CBlock_unwraped);
	return NULL;
}

void BufferBlock::set_owner(void** owner_adrs, ){
	short lvl = 2;
	#ifdef CDEBUG
		fprintf(stderr, "|-----> BufferBlock::set_owner(owner_adrs=%p)\n", owner_adrs);
	#endif

	
	Owner_p = owner_adrs;
	
#ifdef CDEBUG
	fprintf(stderr, "<-----| BufferBlock::set_owner()\n");
#endif
}

void BufferBlock::reset(bool forceReset){
	// Resets block attibutes if it's AVAILABLE to be used again.
	short lvl = 2;
#ifdef CDEBUG
	fprintf(stderr, "|-----> [dev_id=%d] BufferBlock::reset(block_id=%d)\n", Parent->dev_id, id);
#endif
	if(State==INVALID || State==AVAILABLE || forceReset){
		if(!lockfree){
		#if defined(FIFO) || defined(MRU) || defined(LRU)
			Parent->InvalidQueue->lock();
			Parent->Queue->lock();
		#endif
			lock();
		}
		PendingReaders = 0;
		PendingWriters = 0;
		free(WritebackData_p);
		WritebackData_p = NULL;
		Available->reset();

		if(forceReset && State==NATIVE){
			Adrs = NULL;
			State = INVALID;
		#if defined(FIFO) || defined(MRU) || defined(LRU)
			Node_LL_p node = Parent->Queue->remove(Parent->Hash[id], true);
			Parent->InvalidQueue->put_last(node, true);
			node->valid=false;
		#endif
		}
		else{
			set_state(INVALID, true);
		}

		if(Owner_p){
			*Owner_p = NULL;
			Owner_p = NULL;
		}

		if(!lockfree){
			unlock();
		#if defined(FIFO) || defined(MRU) || defined(LRU)
			Parent->InvalidQueue->unlock();
			Parent->Queue->unlock();
		#endif
		}
		#ifdef CDEBUG
			if(forceReset)
				fprintf(stderr, "------- [dev_id=%d] BufferBlock::reset(): Block with id=%d forced to be reseted.\n", Parent->dev_id, id);
			else
				fprintf(stderr, "------- [dev_id=%d] BufferBlock::reset(): Block with id=%d reseted.\n", Parent->dev_id, id);
		#endif
	}
	else
		error("[dev_id=%d] BufferBlock::reset(): Reset was called on a %s block.\n", Parent->dev_id, print_state(State));
#ifdef CDEBUG
	fprintf(stderr, "<-----| [dev_id=%d] BufferBlock::reset(block_id=%d)\n", Parent->dev_id, id);
#endif
}

void BufferBlock::allocate(){
	// Allocates a buffer block if not already pointing to some memory (not null!)
	short lvl = 2;
#ifdef CDEBUG
	fprintf(stderr, "|-----> [dev_id=%d] BufferBlock::allocate(block_id=%d)\n", Parent->dev_id, id);
#endif
	if(Adrs == NULL){
		Adrs = CHLMalloc(Size, Parent->dev_id, 1);
#ifdef CDEBUG
		fprintf(stderr, "------- [dev_id=%d] BufferBlock::allocate(block_id=%d): Allocated Adrs = %p\n", Parent->dev_id, id, Adrs);
#endif
	}
	else{
		#ifdef CDEBUG
			fprintf(stderr, "------- [dev_id=%d] BufferBlock::allocate(block_id=%d) -> Supposedly already allocated block...", Parent->dev_id, id);
		#endif
	}
#ifdef CDEBUG
	fprintf(stderr, "<-----| [dev_id=%d] BufferBlock::allocate(block_id=%d)\n", Parent->dev_id, id);
#endif
}

state BufferBlock::get_state(){
	if(id < 0 || id >= Parent->BlockNum)
		error("[dev_id=%d] BufferBlock::get_state(): Invalid block id=%d\n", Parent->dev_id, id);
	return State;
}

state BufferBlock::set_state(state new_state, ){
	// Forces a new state.
	short lvl = 2;
#ifdef CDEBUG
	fprintf(stderr, "|-----> [dev_id=%d] BufferBlock::set_state(block_id=%d, prior_state=%s)\n", Parent->dev_id, id, print_state(State));
#endif
	if(id < 0 || id >= Parent->BlockNum)
		error("[dev_id=%d] BufferBlock::set_state(%s): Invalid block id=%d\n", Parent->dev_id, print_state(new_state), id);
	if(!lockfree){
	#if defined(FIFO) || defined(MRU) || defined(LRU)
		Parent->InvalidQueue->lock();
		Parent->Queue->lock();
	#endif
		lock();
	}

	state old_state = State;
	if(old_state == NATIVE){;
#ifdef CDEBUG
	fprintf(stderr, "<-----| [dev_id=%d] BufferBlock::set_state(block_id=%d, new_state=%s):\
	Tried to set state of NATIVE block, ignoring...\n", Parent->dev_id, id, print_state(new_state));
#endif
	}
	else State = new_state;
#if defined(FIFO) || defined(MRU) || defined(LRU)
	if(State == INVALID && old_state != INVALID && Parent->Hash[id]->valid){
		Node_LL_p node = Parent->Queue->remove(Parent->Hash[id], true);
		Parent->InvalidQueue->put_last(node, true);
		node->valid=false;
	}

#endif

	if(!lockfree){
		unlock();
	#if defined(FIFO) || defined(MRU) || defined(LRU)
		Parent->Queue->unlock();
		Parent->InvalidQueue->unlock();
	#endif
	}

#ifdef CDEBUG
	fprintf(stderr, "<-----| [dev_id=%d] BufferBlock::set_state(block_id=%d, new_state=%s)\n", Parent->dev_id, id, print_state(State));
#endif
	return old_state;
}

int BufferBlock::update_state(){
	// Updates the state of the block. It cannot raise the state but only lower.
	short lvl = 2;
#ifdef CDEBUG
	fprintf(stderr, "|-----> [dev_id=%d] BufferBlock::update_state(block_id=%d)\n", Parent->dev_id, id);
#endif
	int ret = 0;
	if(!lockfree){
	#if defined(FIFO) || defined(MRU) || defined(LRU)
		Parent->InvalidQueue->lock();
		Parent->Queue->lock();
	#endif
		lock();
	}

	state prev_state = State;
	if(PendingWriters > 0){
		if(State!=EXCLUSIVE && State!=NATIVE)
			error("[dev_id=%d] BufferBlock::update_state(): Block has writers but state was %s.\n", Parent->dev_id, print_state(State));
	}
	else if(PendingReaders > 0){
		if(State==EXCLUSIVE){
			; // Do nothing - Not allowed to schedule out EXCLUSIVE blocks, unless we implement a writeback-to-native mechanism
		}
		else if(State!=SHARABLE && State!=NATIVE)
			error("[dev_id=%d] BufferBlock::update_state(): Block has readers but state was %s.\n", Parent->dev_id, print_state(State));
	}
	else if(State == SHARABLE){
		set_state(AVAILABLE, true);
		ret = 1;
	}
	else if(State == EXCLUSIVE){
		; // Do nothing - Not allowed to schedule out EXLUSIVE blocks, unless we implement a writeback-to-native mechanism
	}
#ifdef CDEBUG
	if(ret==1)
		fprintf(stderr, "------- [dev_id=%d] BufferBlock::update_state(block_id=%d): Block state was changed %s -> %s \n", Parent->dev_id, id, print_state(prev_state), print_state(State));
	else
		fprintf(stderr, "------- [dev_id=%d] BufferBlock::update_state(block_id=%d): Block state is still %s \n", Parent->dev_id, id, print_state(State));
#endif
	if(!lockfree){
		unlock();
	#if defined(FIFO) || defined(MRU) || defined(LRU)
		Parent->Queue->unlock();
		Parent->InvalidQueue->unlock();
	#endif
	}
#ifdef CDEBUG
	fprintf(stderr, "<-----| [dev_id=%d] BufferBlock::update_state(block_id=%d)\n", Parent->dev_id, id);
#endif
	return ret;
}

/*********************
 ** Buffer Functions **
 *********************/

Buffer::Buffer(int dev_id_in, long long block_num, long long block_size){
	// Constructor for buffers
	short lvl = 2;
#ifdef CDEBUG
	fprintf(stderr, "|-----> [dev_id=%d] Buffer::Buffer(block_num = %lld, block_size = %lld)\n", dev_id_in, block_num, block_size);
#endif
  Lock = 0;
	id = DevBuffer_ctr++;
	dev_id = dev_id_in;
	BlockSize = block_size;
	SerialCtr = 0;
	BlockNum = block_num;
	Size = BlockSize*BlockNum;
	Blocks =  (CBlock_p*) malloc (BlockNum * sizeof(CBlock_p));
	for (int idx = 0; idx < BlockNum; idx++) Blocks[idx] = new BufferBlock(idx, this, BlockSize); // Or NULL here and initialize when requested? not sure
	cont_buf_head = NULL;

	#if defined(FIFO) || defined(MRU) || defined(LRU)
	Hash = (Node_LL_p*) malloc(BlockNum * sizeof(Node_LL_p));
	InvalidQueue = new LinkedList(this, "InvalidQueue");
	for(int idx = 0; idx < BlockNum; idx++){
		InvalidQueue->push_back(idx, true);
		Hash[idx] = InvalidQueue->end;
		Hash[idx]->valid=false;
	}
	Queue = new LinkedList(this, "Queue");
	#endif
#ifdef CDEBUG
	fprintf(stderr, "<-----| [dev_id=%d] Buffer::Buffer()\n", dev_id_in);
#endif
}

Buffer::~Buffer(){
	short lvl = 2;
#ifdef CDEBUG
	fprintf(stderr, "|-----> [dev_id=%d] Buffer::~Buffer()\n", dev_id);
#endif
	lock();
	DevBuffer_ctr--;
#ifdef ENABLE_BUFFER_CONTINUOUS_ALLOC
	long long buff_size = 0;
	for (int idx = 0; idx < BlockNum; idx++) if(Blocks[idx]!=NULL)
		Blocks[idx]->Adrs = NULL;
	//if(cont_buf_head)
	if(cont_buf_head_sz) CHLFree(cont_buf_head, dev_id, cont_buf_head_sz);
	cont_buf_head = NULL;
#endif
	for (int idx = 0; idx < BlockNum; idx++) delete Blocks[idx];
	free(Blocks);
#if defined(FIFO)
	free(Hash);
	delete InvalidQueue;
	delete Queue;
#endif
	unlock();
#ifdef CDEBUG
	fprintf(stderr, "<-----| [dev_id=%d] Buffer::~Buffer()\n", dev_id);
#endif
	return ;
}

void Buffer::reset(, bool forceReset){
	short lvl = 2;
#ifdef CDEBUG
	fprintf(stderr, "|-----> [dev_id=%d] Buffer::reset()\n", dev_id);
#endif
	
	for (int idx = 0; idx < BlockNum; idx++) Blocks[idx]->reset(lockfree, forceReset);
#ifdef STEST
	timer = 0; // Keeps total time spend in buffer operations-code
#endif
	SerialCtr = 0;
#if defined(FIFO) || defined(MRU) || defined(LRU)
	if(!lockfree){
		InvalidQueue->lock();
		Queue->lock();
	}
	Node_LL_p node = Queue->start_iterration();
	while(node->idx != -1){//i < Queue->length){
		node->valid=false;
		node = Queue->next_in_line();
	}
	if(!InvalidQueue->is_empty(true)){// || InvalidQueue->length>0){
		if(!Queue->is_empty(true)){//Queue->length>0){
		InvalidQueue->end->next = Queue->start;
		Queue->start->previous = InvalidQueue->end;
		InvalidQueue->end = Queue->end;
		InvalidQueue->length += Queue->length;
		}
	}
	else{
		InvalidQueue->start = Queue->start;
		InvalidQueue->end = Queue->end;
		InvalidQueue->length = Queue->length;
	}
	Queue->start = NULL;
	Queue->end = NULL;
	Queue->length = 0;
	if(!lockfree){
		InvalidQueue->unlock();
		Queue->unlock();
	}
#endif
	
#ifdef CDEBUG
	fprintf(stderr, "<-----| [dev_id=%d] Buffer::reset()\n", dev_id);
#endif
	return ;
}

void Buffer::draw_buffer(bool print_blocks, bool print_queue, ){
	short lvl = 0;

	
	fprintf(stderr, " Buffer:\
						\n==========================================\
						\n||      Buffer Id       | %d\
						\n|| - - - - - - - - - - - - - - - - - - - -\
						\n||      Device Id      | %d\
						\n|| - - - - - - - - - - - - - - - - - - - -\
						\n||        Name         | %s\
						\n|| - - - - - - - - - - - - - - - - - - - -\
						\n||        Size         | %llu\
						\n|| - - - - - - - - - - - - - - - - - - - -\
						\n||  Number of blocks   | %d\
						\n|| - - - - - - - - - - - - - - - - - - - -\
						\n||   Size of blocks    | %llu\
						\n|| - - - - - - - - - - - - - - - - - - - -", id, dev_id, Name.c_str(), Size, BlockNum, BlockSize);
#if defined(NAIVE)
	fprintf(stderr,"\n||  Scheduling policy  | NAIVE");
#elif defined(FIFO)
	fprintf(stderr,"\n||  Scheduling policy  | FIFO");
#elif defined(MRU)
	fprintf(stderr,"\n||  Scheduling policy  | MRU");
#elif defined(LRU)
	fprintf(stderr,"\n||  Scheduling policy  | LRU");
#endif

	fprintf(stderr, "\n==========================================");

	if(print_blocks){
		fprintf(stderr, "======================================\
							\n|| Start of Blocks in Buffer ||\
							\n==============================\n");
		for(int i=0; i<BlockNum; i++)
			if(Blocks[i]!=NULL)
				Blocks[i]->draw_block(lockfree);
		fprintf(stderr, "============================\
							\n|| End of Blocks in Buffer ||\
							\n================================================================================\n");

	}
	if(print_queue){
		#if defined(NAIVE)
		fprintf(stderr, "There is no Queue.\n");
		#elif defined(FIFO) || defined(MRU) || defined(LRU)
		fprintf(stderr, "|| Start of Queue with Invalid Blocks ||\
							\n=======================================\n");
		InvalidQueue->draw_queue(lockfree);
		fprintf(stderr, "============================\
							\n|| End of Queue with Invalid Blocks ||\
							\n================================================================================\n");
		fprintf(stderr, "|| Start of Queue with Valid Blocks ||\
							\n=======================================\n");
		Queue->draw_queue(lockfree);
		fprintf(stderr, "============================\
							\n|| End of Queue with Valid Blocks ||\
							\n================================================================================\n");
		#endif
	}
	fprintf(stderr, "\n");
	
}

#ifdef ENABLE_BUFFER_CONTINUOUS_ALLOC
/// Allocates all bufferblocks in buffer, but in a single continuous piece of memory.
void Buffer::allocate(){
	short lvl = 2;
#ifdef CDEBUG
	fprintf(stderr, "|-----> [dev_id=%d] Buffer::allocate-continuous()\n", dev_id);
#endif
	long long total_sz = 0, total_offset = 0;
	for(int i=0; i<BlockNum; i++) if(Blocks[i]!=NULL && Blocks[i]->Adrs==NULL) total_sz+= Blocks[i]->Size;
	if(!cont_buf_head && total_sz) cont_buf_head = CHLMalloc(total_sz, dev_id, 1);
	cont_buf_head_sz = total_sz;
	for(int i=0; i<BlockNum; i++)
		if(Blocks[i]!=NULL){
			if(Blocks[i]->Adrs==NULL){
				Blocks[i]->Adrs = cont_buf_head + total_offset;
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
	// Allocates all bufferblocks in buffer
	short lvl = 2;
#ifdef CDEBUG
	fprintf(stderr, "|-----> [dev_id=%d] Buffer::allocate()\n", dev_id);
#endif
	for(int i=0; i<BlockNum; i++)
		if(Blocks[i]!=NULL) Blocks[i]->allocate(lockfree);
		else error("[dev_id=%d] Buffer::allocate() -> Blocks[%d] was NULL\n", dev_id, i);
#ifdef CDEBUG
	fprintf(stderr, "<-----| [dev_id=%d] Buffer::allocate()\n", dev_id);
#endif
}
#endif



CBlock_p Buffer::assign_Cblock(state start_state, ){
	// Assigns a block from buffer to be used for memory.
	// State options are:
	// - INVALID: Raise error.
	// - NATIVE:
	// - EXCLUSIVE: Will add one writer.
	// - SHARABLE: Will add one reader.
	// - AVAILABLE: Will be initialized as AVAILABLE and someone might take it.

	short lvl = 2;
#ifdef CDEBUG
	fprintf(stderr, "|-----> [dev_id=%d] Buffer::assign_Cblock()\n", dev_id);
#endif
	CBlock_p result = NULL;
	 // Lock buffer
#if defined(NAIVE)
	if (SerialCtr >= BlockNum){
#endif
		int remove_block_idx = -42;
	#if defined(NAIVE)
		remove_block_idx = BufferSelectBlockToRemove_naive(this, false);
	#elif defined(FIFO) || defined(MRU) || defined(LRU)
		Node_LL_p remove_block;
		remove_block = BufferSelectBlockToRemove_fifo_mru_lru(this, false);
		remove_block_idx = remove_block->idx;
	#endif
		if(remove_block_idx < 0){ // Check again
		#if defined(NAIVE)
			remove_block_idx = BufferSelectBlockToRemove_naive(this, false);
		#elif defined(FIFO) || defined(MRU) || defined(LRU)
			remove_block = BufferSelectBlockToRemove_fifo_mru_lru(this, false);
			remove_block_idx = remove_block->idx;
		#endif
			if(remove_block_idx < 0){ // Check for exclusive
			#if defined(NAIVE)
				remove_block_idx = BufferSelectExclusiveBlockToRemove_naive(this, false);
			#elif defined(FIFO) || defined(MRU) || defined(LRU)
				remove_block = BufferSelectExclusiveBlockToRemove_fifo_mru_lru(this, false);
				remove_block_idx = remove_block->idx;
			#endif
			}
		}
		if(remove_block_idx >= 0){
	#if defined(FIFO)
			Queue->put_last(remove_block, false);
	#elif defined(MRU)
			Queue->put_first(remove_block, false);
	#elif defined(LRU)
			Queue->put_last(remove_block, false);
	#endif
			result = Blocks[remove_block_idx];
			if(!lockfree){
			#if defined(FIFO) || defined(MRU) || defined(LRU)
			 	InvalidQueue->lock();
				Queue->lock();
			#endif
				result->lock();
			}
			result->reset(true,false);
	#ifdef CDUBUG
		fprintf(stderr,"------ [dev_id=%d] Buffer::assign_Cblock(): Block with id=%d reseted.\n", dev_id, remove_block_idx);
	#endif
			// Set state
			if(start_state==INVALID)
				error("[dev_id=%d] Buffer::assign_Cblock(): New block called to be initialized as invalid\n", dev_id);
			else if(start_state==NATIVE)
				result->set_state(NATIVE, true);
			else if(start_state==EXCLUSIVE){
				result->set_state(EXCLUSIVE, true);
				result->add_writer(true);
			}
			else if(start_state==SHARABLE){
				result->set_state(SHARABLE, true);
				result->add_reader(true);
			}
			else if(start_state==AVAILABLE){
				result->set_state(AVAILABLE, true);
			}
			else
				error("[dev_id=%d] Buffer::assign_Cblock(): Uknown state(%s)\n", dev_id, print_state(start_state));
		#if defined(FIFO) || defined(MRU) || defined(LRU)
			Hash[result->id]->valid = true;
		#endif
			if(!lockfree){
				result->unlock();
			#if defined(FIFO) || defined(MRU) || defined(LRU)
				Queue->unlock();
				InvalidQueue->unlock();
			#endif
				unlock(); // Unlock buffer
			}
		}
		else{
	#if defined(FIFO) || defined(MRU) || defined(LRU)
			delete(remove_block);
	#endif

			if(!lockfree)  unlock(); // Unlock buffer
			result = assign_Cblock(start_state, lockfree);
		}
#if defined(NAIVE)
	}
	else{
		result = Blocks[SerialCtr];
		if(!lockfree)
			result->lock();

		result->reset(true,false);
		SerialCtr++;
		// Set state
		if(start_state==INVALID)
			error("[dev_id=%d] Buffer::assign_Cblock(): New block called to be initialized as invalid\n", dev_id);
		else if(start_state==NATIVE)
			result->set_state(NATIVE, true);
		else if(start_state==EXCLUSIVE){
			result->set_state(EXCLUSIVE, true);
			result->add_writer(true);
		}
		else if(start_state==SHARABLE){
			result->set_state(SHARABLE, true);
			result->add_reader(true);
		}
		else if(start_state==AVAILABLE){
			result->set_state(AVAILABLE, true);
		}
		else
			error("[dev_id=%d] Buffer::assign_Cblock(): Uknown state(%s)\n", dev_id, print_state(start_state));

		if(!lockfree){
			result->unlock();
			unlock(); // Unlock buffer
		}
	}
#endif

#ifdef CDEBUG
	fprintf(stderr, "<-----| [dev_id=%d] Buffer::assign_Cblock()\n", dev_id);
#endif
  return result;
}

void Buffer::lock(){
	while(__sync_lock_test_and_set(&Lock, 1));
	// Lock++;
	// Lock.lock();
}

void Buffer::unlock(){
	__sync_lock_release(&Lock);
	// Lock--;
}

bool Buffer::is_locked(){
	if(Lock==0)
		return false;
	return true;
}

/*********************
 ** Other Functions **
 *********************/

#if defined(NAIVE)
int BufferSelectBlockToRemove_naive(Buffer_p buffer, ){
	short lvl = 2;

	if (buffer == NULL)
		error("BufferSelectBlockToRemove_naive(): Called on empty buffer\n");

#ifdef CDEBUG
	fprintf(stderr, "|-----> [dev_id=%d] BufferSelectBlockToRemove_naive()\n",buffer->dev_id);
#endif

	int result_idx = -1;

	for (int idx = 0; idx < buffer->BlockNum; idx++){ // Iterate through buffer serially.
		if(!lockfree)
			buffer->Blocks[idx]->lock();
		buffer->Blocks[idx]->update_state(true);
		state tmp_state = buffer->Blocks[idx]->get_state(); // Update all events etc for idx.
		if(tmp_state == INVALID || tmp_state == AVAILABLE){ // Indx can be removed if there are no pending events.
			result_idx = idx;
			buffer->Blocks[idx]->set_state(INVALID, true);
			if(!lockfree)
				buffer->Blocks[idx]->unlock();
		#ifdef CDEBUG
			fprintf(stderr, "------- [dev_id=%d] BufferSelectBlockToRemove_naive(): Found available block. Invalidated.\n",buffer->dev_id);
		#endif
			break;
		}
		if(!lockfree)
			buffer->Blocks[idx]->unlock();
	}
	#ifdef CDEBUG
		fprintf(stderr, "<-----| [dev_id=%d] BufferSelectBlockToRemove_naive()\n",buffer->dev_id);
	#endif
	return result_idx;
}

int BufferSelectExclusiveBlockToRemove_naive(Buffer_p buffer, ){
	short lvl = 2;

	if (buffer == NULL)
		error("BufferSelectExclusiveBlockToRemove_naive(): Called on empty buffer\n");

#ifdef CDEBUG
	fprintf(stderr, "|-----> [dev_id=%d] BufferSelectExclusiveBlockToRemove_naive()\n",buffer->dev_id);
#endif

	int result_idx = -1;

	for (int idx = 0; idx < buffer->BlockNum; idx++){ // Iterate through buffer serially.
		if(!lockfree)
			buffer->Blocks[idx]->lock();
		buffer->Blocks[idx]->update_state(true);
		state tmp_state = buffer->Blocks[idx]->get_state(); // Update all events etc for idx.
		if(tmp_state == EXCLUSIVE){
			CBlock_p native_block = buffer->Blocks[idx]->WritebackData_p->Native_block;
			if(!lockfree)
				native_block->lock();
			if(buffer->Blocks[idx]->PendingReaders==0 && buffer->Blocks[idx]->PendingWriters==0){
				result_idx = idx;
				native_block->add_writer(true);
				buffer->Blocks[idx]->add_reader(true);
				if(!lockfree){
					native_block->unlock();
					buffer->Blocks[idx]->unlock();
				}
				buffer->Blocks[idx]->write_back(true);
				if(!lockfree){
					buffer->Blocks[idx]->lock();
					native_block->lock();
				}
				native_block->remove_writer(true);
				// buffer->Blocks[idx]->write_back(true);
				buffer->Blocks[idx]->set_state(INVALID, true);
				if(!lockfree) buffer->Blocks[idx]->unlock();

			#ifdef CDEBUG
				fprintf(stderr, "------- [dev_id=%d] BufferSelectExclusiveBlockToRemove_naive(): Found exclusive block with no pernding operations on it. Invalidated.\n",buffer->dev_id);
			#endif
			}
			if(!lockfree)
				native_block->unlock();
		}
		if(!lockfree)
			buffer->Blocks[idx]->unlock();
	}
	#ifdef CDEBUG
		fprintf(stderr, "<-----| [dev_id=%d] BufferSelectExclusiveBlockToRemove_naive()\n",buffer->dev_id);
	#endif
	return result_idx;
}

#elif defined(FIFO) || defined(MRU) || defined(LRU)
Node_LL_p BufferSelectBlockToRemove_fifo_mru_lru(Buffer_p buffer, ){
	short lvl = 2;
	if (buffer == NULL)
		error("BufferSelectBlockToRemove_fifo_mru_lru(): Called on empty buffer\n");

#ifdef CDEBUG
	fprintf(stderr, "|-----> [dev_id=%d] BufferSelectBlockToRemove_fifo_mru_lru()\n", buffer->dev_id);
#endif
	Node_LL_p result_node;
	if(!lockfree){
		buffer->InvalidQueue->lock();
		buffer->Queue->lock();
	}
	if(!buffer->InvalidQueue->is_empty(true)){//buffer->InvalidQueue->length > 0){
		result_node = buffer->InvalidQueue->remove(buffer->InvalidQueue->start, true);
		if(!lockfree) buffer->Blocks[result_node->idx]->lock();
		buffer->Blocks[result_node->idx]->set_state(INVALID, true);
		if(!lockfree) buffer->Blocks[result_node->idx]->unlock();
	}
	else{
		result_node = new Node_LL();
		result_node->idx = -1;
		state tmp_state = INVALID;
		// if(!lockfree){
		// 	buffer->InvalidQueue->lock();
		// 	buffer->Queue->lock();
		// }
		Node_LL_p node = buffer->Queue->start_iterration();
		// int i=0;
		if(node->idx >= 0){
			if(!lockfree)
				buffer->Blocks[node->idx]->lock();
			tmp_state = buffer->Blocks[node->idx]->get_state(); // Update all events etc for idx.
			while(tmp_state != AVAILABLE){
				if(!lockfree)
					buffer->Blocks[node->idx]->unlock();
				node = buffer->Queue->next_in_line();
				if(node->idx >= 0){// && i < buffer->Queue->length){
					if(!lockfree)
						buffer->Blocks[node->idx]->lock();
					tmp_state = buffer->Blocks[node->idx]->get_state(); // Update all events etc for idx.
					// i++;
				}
				else
					break;
			}
		}
		if(node->idx >=0){// && i < buffer->Queue->length){
			if(tmp_state == AVAILABLE){
				delete(result_node);
				buffer->Blocks[node->idx]->set_state(INVALID, true);
				result_node = buffer->InvalidQueue->remove(node, true);
			#ifdef CDEBUG
				fprintf(stderr, "------- [dev_id=%d] BufferSelectBlockToRemove_fifo_mru_lru(): Found available block. Invalidated.\n",buffer->dev_id);
			#endif
			}
			if(!lockfree)
				buffer->Blocks[result_node->idx]->unlock();
		}
		else{//  if(i >= buffer->Queue->length){
			result_node = new Node_LL();
			result_node->idx = -1;
			result_node->valid = false;
		}
	}
	if(!lockfree){
		buffer->Queue->unlock();
		buffer->InvalidQueue->unlock();
	}
#ifdef CDEBUG
	fprintf(stderr, "<-----| [dev_id=%d] BufferSelectBlockToRemove_fifo_mru_lru()\n",buffer->dev_id);
#endif
	return result_node;
}

Node_LL_p BufferSelectExclusiveBlockToRemove_fifo_mru_lru(Buffer_p buffer, ){
	short lvl = 2;
	if (buffer == NULL)
		error("BufferSelectExclusiveBlockToRemove_fifo_mru_lru(): Called on empty buffer\n");

#ifdef CDEBUG
	fprintf(stderr, "|-----> [dev_id=%d] BufferSelectExclusiveBlockToRemove_fifo_mru_lru()\n", buffer->dev_id);
#endif
	Node_LL_p result_node;
	if(!lockfree){
		buffer->InvalidQueue->lock();
		buffer->Queue->lock();
	}
	if(!buffer->InvalidQueue->is_empty(true)){
		result_node = buffer->InvalidQueue->remove(buffer->InvalidQueue->start, true);
		if(!lockfree) buffer->Blocks[result_node->idx]->lock();
		buffer->Blocks[result_node->idx]->set_state(INVALID, true);
		if(!lockfree) buffer->Blocks[result_node->idx]->unlock();
	}
	else{
		result_node = new Node_LL();
		result_node->idx = -1;
		state tmp_state = INVALID;
		// if(!lockfree){
		// 	buffer->InvalidQueue->lock();
		// 	buffer->Queue->lock();
		// }
		Node_LL_p node = buffer->Queue->start_iterration();
		// int i=0;
		if(node->idx >= 0){
			if(!lockfree)
				buffer->Blocks[node->idx]->lock();
			tmp_state = buffer->Blocks[node->idx]->get_state(); // Update all events etc for idx.
			while(tmp_state != EXCLUSIVE || buffer->Blocks[node->idx]->PendingReaders>0 || buffer->Blocks[node->idx]->PendingWriters>0){
				if(!lockfree)
					buffer->Blocks[node->idx]->unlock();
				node = buffer->Queue->next_in_line();
				if(node->idx >= 0){// && i < buffer->Queue->length){
					if(!lockfree)
						buffer->Blocks[node->idx]->lock();
					tmp_state = buffer->Blocks[node->idx]->get_state(); // Update all events etc for idx.
					// i++;
				}
				else
					break;
			}
		}
		if(node->idx >=0){// && i < buffer->Queue->length){
			if(tmp_state == EXCLUSIVE){
				CBlock_p native_block = buffer->Blocks[node->idx]->WritebackData_p->Native_block;
				if(!lockfree)
					native_block->lock();
				if(buffer->Blocks[node->idx]->PendingReaders==0 && buffer->Blocks[node->idx]->PendingWriters==0){
					delete(result_node);
					native_block->add_writer(true);
					buffer->Blocks[node->idx]->add_reader(true);
					if(!lockfree){
						native_block->unlock();
						buffer->Blocks[node->idx]->unlock();
						buffer->Queue->unlock();
						buffer->InvalidQueue->unlock();
					}
					buffer->Blocks[node->idx]->write_back(true);
					if(!lockfree){
						buffer->InvalidQueue->lock();
						buffer->Queue->lock();
						buffer->Blocks[node->idx]->lock();
						native_block->lock();
					}
					native_block->remove_writer(true);
					buffer->Blocks[node->idx]->set_state(INVALID, true);
					result_node = buffer->InvalidQueue->remove(node, true);
				#ifdef CDEBUG
					fprintf(stderr, "------- [dev_id=%d] BufferSelectExclusiveBlockToRemove_fifo_mru_lru(): Found exclusive block with no pernding operations on it. Invalidated.\n",buffer->dev_id);
				#endif
				}
				if(!lockfree)
					native_block->unlock();
			}
			if(!lockfree)
				buffer->Blocks[result_node->idx]->unlock();
		}
		else{// if(i >= buffer->Queue->length){
			result_node = new Node_LL();
			result_node->idx = -1;
			result_node->valid = false;
		}
	}
	if(!lockfree){
		buffer->Queue->unlock();
		buffer->InvalidQueue->unlock();
	}
#ifdef CDEBUG
	fprintf(stderr, "<-----| [dev_id=%d] BufferSelectExclusiveBlockToRemove_fifo_mru_lru()\n",buffer->dev_id);
#endif
	return result_node;
}
#endif