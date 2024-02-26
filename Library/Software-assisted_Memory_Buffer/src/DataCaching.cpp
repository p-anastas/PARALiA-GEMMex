///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
/// \author Theodoridis Aristomenis (theodoridisaristomenis@gmail.com)
///
/// \brief The CHLpeLia caching functions.
///

//#include <cassert>
//#include <atomic>

#include "DataCaching.hpp"
//#include "backend_wrappers.hpp"


//Buffer_p Global_Buffer[64] = {NULL};
//Buffer_p Global_Buffer_1D[64] = {NULL};
//Buffer_p Global_Buffer_2D[64] = {NULL};
Buffer_p current_SAB[64] = {NULL};
int CBlock_ctr[64] = {0};
int DevBuffer_ctr = 0;

int globalock = 0;

#if defined(FIFO) || defined(MRU) || defined(LRU)

/**************************
 ** LinkedList Functions **
 **************************/

LinkedList::LinkedList(Buffer_p buffer, std::string name){
	short lvl = 2;
	if(buffer==NULL)
		error("LinkedList::LinkedList(): Creating buffer that doesn't belong to a buffer.\n");
#ifdef CDEBUG
	lprintf(lvl-1, "|-----> [dev_id=%d] LinkedList::LinkedList(name=%s):\n", buffer->dev_id, name.c_str());
#endif
	lock_ll = 0;
	Parent = buffer;
	Name = name;
	start = NULL;
	end = NULL;
	length = 0;
	iter = NULL;
#ifdef CDEBUG
	lprintf(lvl-1, "<-----| [dev_id=%d] LinkedList::LinkedList(name=%s)\n", Parent->dev_id, Name.c_str());
#endif
}

// Rewrited
LinkedList::~LinkedList(){
	short lvl = 2;
#ifdef CDEBUG
lprintf(lvl-1, "|-----> [dev_id=%d] LinkedList::~LinkedList(name=%s)\n", Parent->dev_id, Name.c_str());
#endif
	lock();
	iter = end;
	while(iter!=NULL && length>0){
		end = iter->previous;
		delete(iter);
		iter = end;
		length--;
	}
	unlock();
#ifdef CDEBUG
lprintf(lvl-1, "<-----| [dev_id=%d] LinkedList::~LinkedList(name=%s)\n", Parent->dev_id, Name.c_str());
#endif
}

// Old one. Hasn't changed.
void LinkedList::draw_queue(bool lockfree){
	short lvl = 1;
	if(!lockfree)
		lock();
	int count = 1;
#if defined(FIFO)
	lprintf(lvl-1, " FIFO");
#elif defined(MRU)
	lprintf(lvl-1, " MRU");
#elif defined(LRU)
	lprintf(lvl-1, " LRU");
#endif
	lprintf(lvl-1, " Queue:\
						\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
						\n||      Buffer Id       | %d\
						\n|| - - - - - - - - - - - - - - - - - - - -\
						\n||      Device Id      | %d\
						\n|| - - - - - - - - - - - - - - - - - - - -\
						\n||        Name         | %s\
						\n|| - - - - - - - - - - - - - - - - - - - -\
						\n||       Length        | %d\
						\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n", Parent->id, Parent->dev_id, Name.c_str(), length);
	if(length>0){
	iter = start;
	lprintf(lvl-1, " Position: %d\
						\n_________________________________________\
						\n|       Idx       | %d\
						\n| - - - - - - - - - - - - - - - - - - - -\
						\n|      Valid      | %s\
						\n|________________________________________\
						\n", count, iter->idx, iter->valid ? "True" : "False");
		iter = iter->next;
		count++;
		while(iter!=NULL && count<=length){
			lprintf(lvl-1, " Position: %d\
								\n_________________________________________\
								\n|       Idx       | %d\
								\n| - - - - - - - - - - - - - - - - - - - -\
								\n|      Valid      | %s\
								\n|________________________________________\
								\n", count, iter->idx, iter->valid ? "True" : "False");
			iter = iter->next;
			count++;
		}
	}
	if(!lockfree)
		unlock();
}

//  Rewrited
void LinkedList::invalidate(Node_LL_p node, bool lockfree){
	short lvl = 2;
#ifdef CDEBUG
	lprintf(lvl-1, "|-----> [dev_id=%d] LinkedList::invalidate(name=%s, node_id=%d)\n", Parent->dev_id, Name.c_str(), node->idx);
#endif
	if(!lockfree)
		lock();
	node->valid = false;
	if(!lockfree)
		unlock();
#ifdef CDEBUG
	lprintf(lvl-1, "<-----| [dev_id=%d] LinkedList::invalidate(name=%s, node_id=%d)\n", Parent->dev_id, Name.c_str(), node->idx);
#endif
}

bool LinkedList::is_empty(bool lockfree){
	short lvl = 2;
#ifdef CDEBUG
lprintf(lvl-1, "|-----> [dev_id=%d] LinkedList::is_empty(name=%s)\n", Parent->dev_id, Name.c_str());
#endif
	bool ret = false;
	if(!lockfree) lock();
	if(start==NULL && end==NULL)
		ret = true;
	else if(start!=NULL && end!=NULL)
		ret = false;
	else // Should never happen
		error("[dev_id=%d] LinkedList::is_empty(name=%s): One of the start and end pointers is null and the other not, which cannot happen.\n", Parent->dev_id, Name.c_str());
	if(!lockfree) unlock();
#ifdef CDEBUG
	if(ret)
		lprintf(lvl-1, "<-----| [dev_id=%d] LinkedList::is_empty(name=%s): List is empty.\n", Parent->dev_id, Name.c_str());
	else
		lprintf(lvl-1, "<-----| [dev_id=%d] LinkedList::is_empty(name=%s): List not empty.\n", Parent->dev_id, Name.c_str());
#endif
	return ret;
}

// Rewritten
void LinkedList::push_back(int idx, bool lockfree){
	// Pushes a new element in the back of the queue.
	short lvl = 2;
#ifdef CDEBUG
lprintf(lvl-1, "|-----> [dev_id=%d] LinkedList::push_back(name=%s, idx=%d)\n", Parent->dev_id, Name.c_str(), idx);
#endif
	if(idx >= Parent->BlockNum || idx < 0)
		error("[dev_id=%d] LinkedList::push_back(name=%s): Tried to push back an idx=%d that is not valid because there are %d blocks in memory.", Parent->dev_id, Name.c_str(), idx, Parent->BlockNum);
	else if(length > Parent->BlockNum)
		error("[dev_id=%d] LinkedList::push_back(name=%s): Called to put another element but max length is reached with BlockNum=%d and length=%d\n", Parent->dev_id, Name.c_str(), Parent->BlockNum, length);
	else{
		if(!lockfree)
			lock();
		// Set new node.

		Node_LL_p node = new Node_LL();
		node->idx = idx;
		node->valid = false;
		// Attach it to the end of the list.
		node->next = NULL;
		if(is_empty(true)){
			node->previous = NULL;
			start = node;
			end = node;
		}
		else{
			node->previous = end;
			end->next = node;
			end = node;
		}
		length++;
		if(!lockfree)
			unlock();
	}
#ifdef CDEBUG
lprintf(lvl-1, "<-----| [dev_id=%d] LinkedList::push_back(name=%s, idx=%d)\n", Parent->dev_id, Name.c_str(), idx);
#endif

}

// Rewrittern
Node_LL_p LinkedList::start_iterration(){
	// Returns the first valid element without removing it.
	short lvl = 2;
#ifdef CDEBUG
lprintf(lvl-1, "|-----> [dev_id=%d] LinkedList::start_iterration(name=%s)\n", Parent->dev_id, Name.c_str());
#endif
	Node_LL_p tmp_node;
	tmp_node = new Node_LL();
	tmp_node->idx = -1;
	tmp_node->valid = false;
	if(!is_empty(true)){
		iter = start;
		while(iter != NULL && !iter->valid)
			iter = iter->next;
		if(iter != NULL){
			tmp_node = iter;
		}
	#ifdef CDEBUG
	lprintf(lvl-1, "<-----| [dev_id=%d] LinkedList::start_iterration(name=%s)\n", Parent->dev_id, Name.c_str());
	#endif
		return tmp_node;
	}
	// else{
		// tmp_node = new Node_LL();
		// tmp_node->idx = -1;
		// tmp_node->valid = false;
	// }
	// error("[dev_id=%d] LinkedList::start_iterration(name=%s): Called to iterrate on empty list.\n", Parent->dev_id, Name.c_str());
#ifdef CDEBUG
lprintf(lvl-1, "<-----| [dev_id=%d] LinkedList::start_iterration(name=%s): Iterration not started.\n", Parent->dev_id, Name.c_str());
#endif
	return tmp_node;
}

// Rewritten
Node_LL_p LinkedList::next_in_line(){
	// Returns next element in iterration.
	short lvl = 2;
#ifdef CDEBUG
lprintf(lvl-1, "|-----> [dev_id=%d] LinkedList::next_in_line(name=%s)\n", Parent->dev_id, Name.c_str());
#endif
	Node_LL_p tmp_node;
	tmp_node = new Node_LL();
	tmp_node->idx = -1;
	tmp_node->valid = false;
	if(iter != NULL){
		iter = iter->next;
		while(iter != NULL && !iter->valid)
			iter = iter->next;
		if(iter != NULL){
			tmp_node = iter;
		}
	#ifdef CDEBUG
	lprintf(lvl-1, "<-----| [dev_id=%d] LinkedList::next_in_line(name=%s)\n", Parent->dev_id, Name.c_str());
	#endif
		return tmp_node;
	}

	warning("[dev_id=%d] LinkedList::start_iterration(name=%s): Called to iterrate on empty list.\n", Parent->dev_id, Name.c_str());
	// #ifdef CDEBUG
	// 	lprintf(lvl-1, "<-----| [dev_id=%d] LinkedList::next_in_line(name=%s): Iterration not started.\n", Parent->dev_id, Name.c_str());
	// #endif
	return tmp_node;
}
// Rewritten
Node_LL_p LinkedList::remove(Node_LL_p node, bool lockfree){
	short lvl = 2;

	#ifdef CDEBUG
		if(node!=NULL)
		lprintf(lvl-1, "|-----> [dev_id=%d] LinkedList::remove(name=%s, idx=%d)\n", Parent->dev_id, Name.c_str(), node->idx);
	#endif
	if(node == NULL) // Node not found.
		error("[dev_id=%d] LinkedList::remove(name=%s): Input node not found. Can't pop element.\n", Parent->dev_id, Name.c_str());
	else if(node->idx >= Parent->BlockNum)
		error("[dev_id=%d] LinkedList::remove(name=%s): Index of given node(%d) is larger than the number of blocks(%d).\n", Parent->dev_id, Name.c_str(), node->idx, Parent->BlockNum);
	else if(node->idx < 0)
		error("[dev_id=%d] LinkedList::remove(name=%s): Index of given node(%d) is negative.\n", Parent->dev_id, Name.c_str(), node->idx);
	else{
		if(!lockfree)
			lock();
		if(node->previous != NULL && node->next != NULL){
			node->previous->next = node->next;
			node->next->previous = node->previous;
		}// If last node
		else if(node->previous != NULL && node->next == NULL){
			node->previous->next = NULL;
			end = node->previous;
		}// If first node
		else if(node->previous == NULL && node->next != NULL){
			node->next->previous = NULL;
			start = node->next;
		}
		else if(node->previous == NULL && node->next == NULL){
			start = NULL;
			end = NULL;
		}
		else{ // Should NEVER happen.
			if(is_empty(true)) error("[dev_id=%d] LinkedList::remove(name=%s): Node not from list. List is empty.\n", Parent->dev_id, Name.c_str());
			error("[dev_id=%d] LinkedList::remove(name=%s): Internal.\n", Parent->dev_id, Name.c_str());
		}

		node->next = NULL;
		node->previous = NULL;
		length--;
		if(!lockfree)
			unlock();
	}

#ifdef CDEBUG
lprintf(lvl-1, "<-----| [dev0_id=%d] LinkedList::remove(name=%s, idx=%d)\n", Parent->dev_id, Name.c_str(), node->idx);
#endif
	return node;
}

void LinkedList::put_first(Node_LL* node, bool lockfree){
	// Puts the element first on the list.
	short lvl = 2;

#ifdef CDEBUG
lprintf(lvl-1, "|-----> [dev_id=%d] LinkedList::put_first(name=%s, idx=%d)\n", Parent->dev_id, Name.c_str(), node->idx);
#endif
	if(node == NULL)
		error("[dev_id=%d] LinkedList::put_first(name=%s): Called to put first NULL.\n", Parent->dev_id, Name.c_str());
	else if(node->idx < 0)
		error("[dev_id=%d] LinkedList::put_first(name=%s): Invalid element with id < 0.\n", Parent->dev_id, Name.c_str());
	else{
		if(!lockfree)
			lock();
		if(is_empty(true)){
			end = node;
		}
		else{
			node->next = start;
			start->previous = node;
		}
		start = node;
		length++;
		if(!lockfree)
			unlock();
	}
#ifdef CDEBUG
lprintf(lvl-1, "<-----| [dev_id=%d] LinkedList::put_first(name=%s, idx=%d)\n", Parent->dev_id, Name.c_str(), node->idx);
#endif
}

void LinkedList::put_last(Node_LL_p node, bool lockfree){
	// Takes a node and puts it in the end of the list.
	short lvl = 2;

#ifdef CDEBUG
lprintf(lvl-1, "|-----> [dev_id=%d] LinkedList::put_last(name=%s, idx=%d)\n", Parent->dev_id, Name.c_str(), node->idx);
#endif
	if(node == NULL)
		error("[dev_id=%d] LinkedList::put_last(name=%s): Called to put last NULL.\n", Parent->dev_id, Name.c_str());
	else if(node->idx < 0)
		error("[dev_id=%d] LinkedList::put_first(name=%s): Invalid element with id < 0.\n", Parent->dev_id, Name.c_str());
	else{
		if(!lockfree)
			lock();
		if(is_empty(true)){
			start = node;
		}
		else{
			node->previous = end;
			end->next = node;
		}
		end = node;
		length++;
		if(!lockfree)
			unlock();
	}
#ifdef CDEBUG
lprintf(lvl-1, "<-----| [dev_id=%d] LinkedList::put_last(name=%s, idx=%d)\n", Parent->dev_id, Name.c_str(), node->idx);
#endif
}

// Old one. Hasn't changed.
void LinkedList::lock(){
	// Int lock
	while(__sync_lock_test_and_set(&lock_ll, 1));
	// Mutex lock
	// lock_ll.lock();

}

// Old one. Hasn't changed.
void LinkedList::unlock(){
	// Int lock
   __sync_lock_release(&lock_ll);
	// Mutex lock
	// lock_ll.unlock();
}

// Old one. Hasn't changed.
bool LinkedList::is_locked(){
	if(lock_ll==0)
		return false;
	return true;
}

#endif

// Old one. Hasn't changed.
const char* print_state(state in_state){
	switch(in_state){
		case(INVALID):
			return "INVALID";
		case(NATIVE):
			return "NATIVE";
		case(EXCLUSIVE):
			return "EXCLUSIVE";
		case(SHARABLE):
			return "SHARABLE";
		case(AVAILABLE):
			return "AVAILABLE";
		default:
			error("print_state: Unknown state\n");
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

	short lvl = 2;

	if(block_parent!=NULL && block_id>=0 && block_id<block_parent->BlockNum && block_size>0){
	#ifdef CDEBUG
		lprintf(lvl-1, "|-----> [dev_id=%d] BufferBlock::BufferBlock(id=%d, Parent_id=%d, Size=%llu)\n", block_parent->dev_id, block_id, block_parent->id, block_size);
	#endif
		id = block_id;
		Parent = block_parent;
		Size = block_size;
		Owner_p = NULL;
		WritebackData_p = NULL;

		PendingReaders = 0;
		PendingWriters = 0;

		Adrs = NULL; // Will be set by cache allocate.
		State = INVALID;
		Lock = 0;

		Available = new Event(Parent->dev_id);
	#ifdef CDEBUG
		lprintf(lvl-1, "<-----| [dev_id=%d] BufferBlock::BufferBlock()\n", Parent->dev_id);
	#endif
	}
	else{
		if(block_parent==NULL)
			error("BufferBlock::BufferBlock(): Constructor called with no buffer to belong.\n");
		else if(block_id<0 && block_id>=block_parent->BlockNum)
			error("[dev_id=%d] BufferBlock::BufferBlock(): Constructor called with invalid id=%d.\n", block_parent->dev_id, block_id);
		else
			error("[dev_id=%d] BufferBlock::BufferBlock(): Constructor called with invalid mem size=%llu.\n", block_parent->dev_id, block_size);
	}
}

BufferBlock::~BufferBlock(){
	short lvl = 2;
	// Destructor of the block.
#ifdef CDEBUG
	lprintf(lvl-1, "|-----> [dev_id=%d] BufferBlock::~BufferBlock()\n", Parent->dev_id);
#endif
	lock();
	if(Owner_p){
		*Owner_p = NULL;
		Owner_p = NULL;
	}
	free(WritebackData_p);
	if(State != NATIVE){
		CHLFree(Adrs, Parent->dev_id, Size);
		delete Available;
#ifdef CDEBUG
		lprintf(lvl-1, "------- [dev_id=%d] BufferBlock::~BufferBlock(): Deleting non-NATIVE block id =%d\n",
			Parent->dev_id, id);
#endif
	}
	else{;
#ifdef CDEBUG
		lprintf(lvl-1, "------- [dev_id=%d] BufferBlock::~BufferBlock(): Refrain from deleting NATIVE block id =%d\n",
			Parent->dev_id, id);
#endif
	}
	unlock();
#ifdef CDEBUG
	lprintf(lvl-1, "<-----| [dev_id=%d] BufferBlock::~BufferBlock()\n", Parent->dev_id);
#endif
}

void BufferBlock::draw_block(bool lockfree){
	// Draws the block for debugging purposes.
	short lvl=0;

	if(!lockfree)
		lock();
	lprintf(lvl-1, " Block:   \
						\n_________________________________________\
						\n|       Id        | %d\
						\n| - - - - - - - - - - - - - - - - - - - -\
						\n|      Name       | %s\
						\n| - - - - - - - - - - - - - - - - - - - -\
						\n|      Size       | %llu\
						\n| - - - - - - - - - - - - - - - - - - - -\
						\n|      State      | %s \
						\n| - - - - - - - - - - - - - - - - - - - -\
						\n| Pending Readers | %d\
						\n| - - - - - - - - - - - - - - - - - - - -\
						\n| Pending Writers | %d\
						\n|________________________________________\
						\n", id, Name.c_str(), Size, print_state(State), PendingReaders.load(), PendingWriters.load());
	if(!lockfree)
		unlock();
}

void BufferBlock::add_reader(bool lockfree){
	short lvl = 2;
#ifdef CDEBUG
	lprintf(lvl-1, "|-----> [dev_id=%d] BufferBlock::add_reader(block_id=%d)\n", Parent->dev_id, id);
#endif
	if(!lockfree){
	#if defined(FIFO) || defined(LRU) || defined(MRU)
		Parent->InvalidQueue->lock();
		Parent->Queue->lock();
	#endif
		lock();
	}
	update_state(true);
	if(State == INVALID)
		error("[dev_id=%d] CacheBlock::add_reader(block_id=%d): Called to put reader on INVALID block\n", Parent->dev_id, id);
	else if(State == AVAILABLE){
		PendingReaders++;
		set_state(SHARABLE, true);
	}
	else if(State == SHARABLE || State == EXCLUSIVE || State == NATIVE)
		PendingReaders++;
	else
		error("[dev_id=%d] BufferBlock::add_reader(): Can't add reader. Block has State=%s\n", Parent->dev_id, print_state(State));
	if(!lockfree){
		unlock();
		#if defined(FIFO) || defined(LRU) || defined(MRU)
		Parent->Queue->unlock();
		Parent->InvalidQueue->unlock();
		#endif
	}
#ifdef CDEBUG
	lprintf(lvl-1, "<-----| [dev_id=%d] BufferBlock::add_reader(block_id=%d)\n", Parent->dev_id, id);
#endif
}

void BufferBlock::add_writer(bool lockfree){
	short lvl = 2;

#ifdef CDEBUG
	lprintf(lvl-1, "|-----> [dev_id=%d] CacheBlock::add_writer(block_id=%d)\n", Parent->dev_id, id);
#endif
	if(!lockfree){
		#if defined(FIFO) || defined(LRU) || defined(MRU)
		Parent->InvalidQueue->lock();
		Parent->Queue->lock();
		#endif
		lock();
	}
#ifdef CDEBUG
	lprintf(lvl-1, "|-----> [dev_id=%d] BufferBlock::add_writer(block_id=%d)\n", Parent->dev_id, id);
#endif
	if(State==EXCLUSIVE || State==NATIVE)
		PendingWriters++;
	else if(State==AVAILABLE || State==SHARABLE){
		PendingWriters++;
		set_state(EXCLUSIVE, true);
	}
	else
		error("[dev_id=%d] BufferBlock::add_reader(): Can't add reader. Block has State=%s\n", Parent->dev_id, print_state(State));
	if(!lockfree){
		unlock();
		#if defined(FIFO) || defined(LRU) || defined(MRU)
		Parent->Queue->unlock();
		Parent->InvalidQueue->unlock();
		#endif
	}
#ifdef CDEBUG
	lprintf(lvl-1, "<-----| [dev_id=%d] BufferBlock::add_writer(block_id=%d)\n", Parent->dev_id, id);
#endif
}

void BufferBlock::remove_reader(bool lockfree){
	short lvl = 2;
#ifdef CDEBUG
	lprintf(lvl-1, "|-----> [dev_id=%d] BufferBlock::remove_reader(block_id=%d)\n", Parent->dev_id, id);
#endif
	if(!lockfree){
		#if defined(FIFO) || defined(LRU) || defined(MRU)
		Parent->InvalidQueue->lock();
		Parent->Queue->lock();
		#endif
		lock();
	}
	// update_state(true);
	// if(State == SHARABLE || State == EXCLUSIVE || State == NATIVE){
	if(PendingReaders.load()>0){
		PendingReaders--;
	}
	else
		error("[dev_id=%d] BufferBlock::remove_reader(): Can't remove reader. There are none.\n", Parent->dev_id);
	update_state(true);
#if defined(MRU)
	Node_LL_p node = Parent->Queue->remove(Parent->Hash[id], true);
	Parent->Queue->put_first(node, true);
#elif defined(LRU)
	Node_LL_p node = Parent->Queue->remove(Parent->Hash[id], true);
	Parent->Queue->put_last(node, true);
#endif
	if(!lockfree){
		unlock();
		#if defined(FIFO) || defined(LRU) || defined(MRU)
		Parent->Queue->unlock();
		Parent->InvalidQueue->unlock();
		#endif
	}
#ifdef CDEBUG
	lprintf(lvl-1, "<-----| [dev_id=%d] BufferBlock::remove_reader(block_id=%d)\n", Parent->dev_id, id);
#endif
}

void BufferBlock::remove_writer(bool lockfree){
	short lvl = 2;
#ifdef CDEBUG
	lprintf(lvl-1, "|-----> [dev_id=%d] BufferBlock::remove_writer(block_id=%d)\n", Parent->dev_id, id);
#endif
	if(!lockfree){
#if defined(FIFO) || defined(MRU) || defined(LRU)
		Parent->InvalidQueue->lock();
		Parent->Queue->lock();
#endif
		lock();
	}
	// if(State == EXCLUSIVE || State == NATIVE){
	if(PendingWriters.load()>0)
		PendingWriters--;
	else
		error("[dev_id=%d] BufferBlock::remove_writer(block_id=%d): Can't remove writer. There are none.\n", Parent->dev_id, id);
	update_state(true);
#if defined(MRU)
	Node_LL_p node = Parent->Queue->remove(Parent->Hash[id], true);
	Parent->Queue->put_first(node, true);
#elif defined(LRU)
	Node_LL_p node = Parent->Queue->remove(Parent->Hash[id], true);
	Parent->Queue->put_last(node, true);
#endif
	if(!lockfree){
		unlock();
#if defined(FIFO) || defined(MRU) || defined(LRU)
		Parent->Queue->unlock();
		Parent->InvalidQueue->unlock();
#endif
	}
#ifdef CDEBUG
	lprintf(lvl-1, "<-----| [dev_id=%d] BufferBlock::remove_writer(block_id=%d)\n", Parent->dev_id, id);
#endif
}

// Old one. Hasn't changed.
void* CBlock_RR_wrap(void* CBlock_wraped){
	CBlock_wrap_p CBlock_unwraped = (CBlock_wrap_p) CBlock_wraped;
	CBlock_unwraped->CBlock->remove_reader(CBlock_unwraped->lockfree);
	free(CBlock_unwraped);
	return NULL;
}

// Old one. Hasn't changed.
void* CBlock_RW_wrap(void* CBlock_wraped){
	CBlock_wrap_p CBlock_unwraped = (CBlock_wrap_p) CBlock_wraped;
	//if(!CBlock_unwraped->CBlock->PendingWriters.load())
	//	printf("|-----> CBlock_RW_wrap: suspicious PendingWriters = %d\n", CBlock_unwraped->CBlock->PendingWriters.load());
	CBlock_unwraped->CBlock->remove_writer(CBlock_unwraped->lockfree);
	free(CBlock_unwraped);
	return NULL;
}

// Old one. Hasn't changed.
void* CBlock_INV_wrap(void* CBlock_wraped){
	CBlock_wrap_p CBlock_unwraped = (CBlock_wrap_p) CBlock_wraped;
	CBlock_unwraped->CBlock->Available->soft_reset();
	CBlock_unwraped->CBlock->set_state(INVALID, CBlock_unwraped->lockfree);
	free(CBlock_unwraped);
	return NULL;
}

// Old one. Hasn't changed.
void* CBlock_RR_INV_wrap(void* CBlock_wraped){
	CBlock_wrap_p CBlock_unwraped = (CBlock_wrap_p) CBlock_wraped;
	if(!CBlock_unwraped->lockfree){
	#if defined(FIFO) || defined(MRU) || defined(LRU)
		CBlock_unwraped->CBlock->Parent->InvalidQueue->lock();
		CBlock_unwraped->CBlock->Parent->Queue->lock();
	#endif
		CBlock_unwraped->CBlock->lock();
	}
	CBlock_unwraped->CBlock->remove_reader(true);
	CBlock_unwraped->CBlock->Available->soft_reset();
	CBlock_unwraped->CBlock->set_state(INVALID, true);

	if(!CBlock_unwraped->lockfree){
		CBlock_unwraped->CBlock->unlock();
	#if defined(FIFO) || defined(MRU) || defined(LRU)
		CBlock_unwraped->CBlock->Parent->Queue->unlock();
		CBlock_unwraped->CBlock->Parent->InvalidQueue->unlock();
	#endif
	}
	free(CBlock_unwraped);
	return NULL;
}

// Old one. Hasn't changed.
void* CBlock_RW_INV_wrap(void* CBlock_wraped){
	CBlock_wrap_p CBlock_unwraped = (CBlock_wrap_p) CBlock_wraped;
	if(!CBlock_unwraped->lockfree){
	#if defined(FIFO) || defined(MRU) || defined(LRU)
		CBlock_unwraped->CBlock->Parent->InvalidQueue->lock();
		CBlock_unwraped->CBlock->Parent->Queue->lock();
	#endif
		CBlock_unwraped->CBlock->lock();
	}

	//if(!CBlock_unwraped->CBlock->PendingWriters.load())
	//	printf("|-----> CBlock_RW_INV_wrap: suspicious PendingWriters = %d\n", CBlock_unwraped->CBlock->PendingWriters.load());
	CBlock_unwraped->CBlock->remove_writer(true);
	CBlock_unwraped->CBlock->Available->soft_reset();
	CBlock_unwraped->CBlock->set_state(INVALID, true);

	if(!CBlock_unwraped->lockfree){
		CBlock_unwraped->CBlock->unlock();
	#if defined(FIFO) || defined(MRU) || defined(LRU)
		CBlock_unwraped->CBlock->Parent->Queue->unlock();
		CBlock_unwraped->CBlock->Parent->InvalidQueue->unlock();
	#endif
	}
	free(CBlock_unwraped);
	return NULL;
}

void BufferBlock::set_owner(void** owner_adrs, bool lockfree){
	short lvl = 2;
	#ifdef CDEBUG
		lprintf(lvl-1, "|-----> BufferBlock::set_owner(owner_adrs=%p)\n", owner_adrs);
	#endif

	if(!lockfree)
		lock();
	Owner_p = owner_adrs;
	if(!lockfree)
		unlock();
#ifdef CDEBUG
	lprintf(lvl-1, "<-----| BufferBlock::set_owner()\n");
#endif
}

void BufferBlock::reset(bool lockfree, bool forceReset){
	// Resets block attibutes if it's AVAILABLE to be used again.
	short lvl = 2;
#ifdef CDEBUG
	lprintf(lvl-1, "|-----> [dev_id=%d] BufferBlock::reset(block_id=%d)\n", Parent->dev_id, id);
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
				lprintf(lvl-1, "------- [dev_id=%d] BufferBlock::reset(): Block with id=%d forced to be reseted.\n", Parent->dev_id, id);
			else
				lprintf(lvl-1, "------- [dev_id=%d] BufferBlock::reset(): Block with id=%d reseted.\n", Parent->dev_id, id);
		#endif
	}
	else
		error("[dev_id=%d] BufferBlock::reset(): Reset was called on a %s block.\n", Parent->dev_id, print_state(State));
#ifdef CDEBUG
	lprintf(lvl-1, "<-----| [dev_id=%d] BufferBlock::reset(block_id=%d)\n", Parent->dev_id, id);
#endif
}

void BufferBlock::init_writeback_info(CBlock_p WB_block, int* RW_master_p,
	int dim1, int dim2, int ldim, int ldim_wb, int dtype_sz, CQueue_p wb_queue, bool lockfree){
	short lvl = 2;
#ifdef CDEBUG
	lprintf(lvl-1, "|-----> [dev_id=%d] BufferBlock::init_writeback_info(block_id=%d)\n", Parent->dev_id, id);
#endif

	if(WritebackData_p!=NULL) error("[dev_id=%d] BufferBlock::init_writeback_info(block_id=%d):\
		Called with a previously allocated WritebackData_p\n", Parent->dev_id, id);

	if(!lockfree){
#if defined(FIFO) || defined(MRU) || defined(LRU)
		Parent->InvalidQueue->lock();
		Parent->Queue->lock();
#endif
		lock();
	}

	WritebackData_p = (writeback_info_p) malloc(sizeof(struct writeback_info));
	WritebackData_p->Native_block = WB_block;
	WritebackData_p->WB_master_p = RW_master_p;
	WritebackData_p->dim1 = dim1;
	WritebackData_p->dim2 = dim2;
	WritebackData_p->ldim = ldim;
	WritebackData_p->ldim_wb = ldim_wb;
	WritebackData_p->dtype_sz = dtype_sz;
	WritebackData_p->wb_queue = wb_queue;
	if(!lockfree){
#if defined(FIFO) || defined(MRU) || defined(LRU)
		Parent->InvalidQueue->unlock();
		Parent->Queue->unlock();
#endif
		unlock();
	}
	#ifdef CDEBUG
	lprintf(lvl-1, "<-----| [dev_id=%d] BufferBlock::init_writeback_info(block_id=%d)\n", Parent->dev_id, id);
	#endif
}

void BufferBlock::write_back(bool lockfree){
	short lvl = 2;

	#ifdef CDEBUG
	lprintf(lvl-1, "|-----> [dev_id=%d] BufferBlock::write_back(block_id=%d)\n", Parent->dev_id, id);
	#endif

	if(WritebackData_p->Native_block==NULL)
	error("[dev_id=%d] BufferBlock::write_back(block_id=%d): Can't write back. Native block is NULL.\n", Parent->dev_id, id);
	else{
		if(!lockfree){
			// #if defined(FIFO) || defined(MRU) || defined(LRU)
			// Parent->InvalidQueue->lock();
			// Parent->Queue->lock();
			// #endif
			lock();
			WritebackData_p->Native_block->lock();
		}
		CBlock_p Write_back_Native_block = WritebackData_p->Native_block;
		/// We always lock for now, since we can't lock externally (since reset also resets WritebackData_p)
		WritebackData_p->wb_queue->memcpy2DAsync(WritebackData_p->Native_block->Adrs, WritebackData_p->ldim_wb, Adrs, WritebackData_p->ldim,
			WritebackData_p->dim1, WritebackData_p->dim2, WritebackData_p->dtype_sz,
			WritebackData_p->Native_block->Parent->dev_id, Parent->dev_id, 1);
		*(WritebackData_p->WB_master_p) = WritebackData_p->Native_block->Parent->dev_id;
		if(!lockfree){
			add_reader(true);
			Write_back_Native_block->add_writer(true);
			Write_back_Native_block->unlock();
			unlock();
		}
		WritebackData_p->wb_queue->sync_barrier();
#ifdef CDEBUG
		lprintf(lvl-1, "------- [dev_id=%d] BufferBlock::write_back(block_id=%d): Writeback complete\n", Parent->dev_id, id);
#endif
		if(!lockfree)
			Write_back_Native_block->remove_writer();
		reset(lockfree, true);
#ifdef CDEBUG
		lprintf(lvl-1, "------- [dev_id=%d] BufferBlock::write_back(block_id=%d): Reset block complete\n", Parent->dev_id, id);
#endif
		if(!lockfree){
			// Write_back_Native_block->unlock();
			// unlock();
			// #if defined(FIFO) || defined(MRU) || defined(LRU)
			// Parent->Queue->unlock();
			// Parent->InvalidQueue->unlock();
			// #endif
		}
	}

	#ifdef CDEBUG
	lprintf(lvl-1, "<-----| [dev_id=%d] BufferBlock::write_back(block_id=%d)\n", Parent->dev_id, id);
	#endif
}

void BufferBlock::allocate(bool lockfree){
	// Allocates a buffer block if not already pointing to some memory (not null!)
	short lvl = 2;
#ifdef CDEBUG
	lprintf(lvl, "|-----> [dev_id=%d] BufferBlock::allocate(block_id=%d)\n", Parent->dev_id, id);
#endif
	if(Adrs == NULL){
		Adrs = CHLMalloc(Size, Parent->dev_id, 1);
#ifdef CDEBUG
		lprintf(lvl, "------- [dev_id=%d] BufferBlock::allocate(block_id=%d): Allocated Adrs = %p\n", Parent->dev_id, id, Adrs);
#endif
	}
	else{
		#ifdef CDEBUG
			lprintf(lvl, "------- [dev_id=%d] BufferBlock::allocate(block_id=%d) -> Supposedly already allocated block...", Parent->dev_id, id);
		#endif
	}
#ifdef CDEBUG
	lprintf(lvl, "<-----| [dev_id=%d] BufferBlock::allocate(block_id=%d)\n", Parent->dev_id, id);
#endif
}

state BufferBlock::get_state(){
	if(id < 0 || id >= Parent->BlockNum)
		error("[dev_id=%d] BufferBlock::get_state(): Invalid block id=%d\n", Parent->dev_id, id);
	return State;
}

state BufferBlock::set_state(state new_state, bool lockfree){
	// Forces a new state.
	short lvl = 2;
#ifdef CDEBUG
	lprintf(lvl-1, "|-----> [dev_id=%d] BufferBlock::set_state(block_id=%d, prior_state=%s)\n", Parent->dev_id, id, print_state(State));
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
	lprintf(lvl-1, "<-----| [dev_id=%d] BufferBlock::set_state(block_id=%d, new_state=%s):\
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
	lprintf(lvl-1, "<-----| [dev_id=%d] BufferBlock::set_state(block_id=%d, new_state=%s)\n", Parent->dev_id, id, print_state(State));
#endif
	return old_state;
}

int BufferBlock::update_state(bool lockfree){
	// Updates the state of the block. It cannot raise the state but only lower.
	short lvl = 2;
#ifdef CDEBUG
	lprintf(lvl-1, "|-----> [dev_id=%d] BufferBlock::update_state(block_id=%d)\n", Parent->dev_id, id);
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
		lprintf(lvl-1, "------- [dev_id=%d] BufferBlock::update_state(block_id=%d): Block state was changed %s -> %s \n", Parent->dev_id, id, print_state(prev_state), print_state(State));
	else
		lprintf(lvl-1, "------- [dev_id=%d] BufferBlock::update_state(block_id=%d): Block state is still %s \n", Parent->dev_id, id, print_state(State));
#endif
	if(!lockfree){
		unlock();
	#if defined(FIFO) || defined(MRU) || defined(LRU)
		Parent->Queue->unlock();
		Parent->InvalidQueue->unlock();
	#endif
	}
#ifdef CDEBUG
	lprintf(lvl-1, "<-----| [dev_id=%d] BufferBlock::update_state(block_id=%d)\n", Parent->dev_id, id);
#endif
	return ret;
}

void BufferBlock::lock(){
	while(__sync_lock_test_and_set(&Lock, 1));
}

void BufferBlock::unlock(){
	__sync_lock_release(&Lock);
}

bool BufferBlock::is_locked(){
	if(Lock==0)
		return false;
	return true;
}

/*********************
 ** Buffer Functions **
 *********************/

Buffer::Buffer(int dev_id_in, long long block_num, long long block_size){
	// Constructor for buffers
	short lvl = 2;
#ifdef CDEBUG
	lprintf(lvl-1, "|-----> [dev_id=%d] Buffer::Buffer(block_num = %lld, block_size = %lld)\n", dev_id_in, block_num, block_size);
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
	lprintf(lvl-1, "<-----| [dev_id=%d] Buffer::Buffer()\n", dev_id_in);
#endif
}

Buffer::~Buffer(){
	short lvl = 2;
#ifdef CDEBUG
	lprintf(lvl-1, "|-----> [dev_id=%d] Buffer::~Buffer()\n", dev_id);
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
	lprintf(lvl-1, "<-----| [dev_id=%d] Buffer::~Buffer()\n", dev_id);
#endif
	return ;
}

void Buffer::reset(bool lockfree, bool forceReset){
	short lvl = 2;
#ifdef CDEBUG
	lprintf(lvl-1, "|-----> [dev_id=%d] Buffer::reset()\n", dev_id);
#endif
	if(!lockfree) lock();
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
	if(!lockfree) unlock();
#ifdef CDEBUG
	lprintf(lvl-1, "<-----| [dev_id=%d] Buffer::reset()\n", dev_id);
#endif
	return ;
}

void Buffer::draw_buffer(bool print_blocks, bool print_queue, bool lockfree){
	short lvl = 0;

	if(!lockfree)
		lock();
	lprintf(lvl-1, " Buffer:\
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
	lprintf(lvl-1,"\n||  Scheduling policy  | NAIVE");
#elif defined(FIFO)
	lprintf(lvl-1,"\n||  Scheduling policy  | FIFO");
#elif defined(MRU)
	lprintf(lvl-1,"\n||  Scheduling policy  | MRU");
#elif defined(LRU)
	lprintf(lvl-1,"\n||  Scheduling policy  | LRU");
#endif

	lprintf(lvl-1, "\n==========================================");

	if(print_blocks){
		lprintf(lvl-1, "======================================\
							\n|| Start of Blocks in Buffer ||\
							\n==============================\n");
		for(int i=0; i<BlockNum; i++)
			if(Blocks[i]!=NULL)
				Blocks[i]->draw_block(lockfree);
		lprintf(lvl-1, "============================\
							\n|| End of Blocks in Buffer ||\
							\n================================================================================\n");

	}
	if(print_queue){
		#if defined(NAIVE)
		lprintf(lvl-1, "There is no Queue.\n");
		#elif defined(FIFO) || defined(MRU) || defined(LRU)
		lprintf(lvl-1, "|| Start of Queue with Invalid Blocks ||\
							\n=======================================\n");
		InvalidQueue->draw_queue(lockfree);
		lprintf(lvl-1, "============================\
							\n|| End of Queue with Invalid Blocks ||\
							\n================================================================================\n");
		lprintf(lvl-1, "|| Start of Queue with Valid Blocks ||\
							\n=======================================\n");
		Queue->draw_queue(lockfree);
		lprintf(lvl-1, "============================\
							\n|| End of Queue with Valid Blocks ||\
							\n================================================================================\n");
		#endif
	}
	lprintf(lvl-1, "\n");
	if(!lockfree)
		unlock();
}

#ifdef ENABLE_BUFFER_CONTINUOUS_ALLOC
/// Allocates all bufferblocks in buffer, but in a single continuous piece of memory.
void Buffer::allocate(bool lockfree){
	short lvl = 2;
#ifdef CDEBUG
	lprintf(lvl-1, "|-----> [dev_id=%d] Buffer::allocate-continuous()\n", dev_id);
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
	lprintf(lvl-1, "<-----| [dev_id=%d] Buffer::allocate-continuous()\n", dev_id);
#endif
}

#else

void Buffer::allocate(bool lockfree){
	// Allocates all bufferblocks in buffer
	short lvl = 2;
#ifdef CDEBUG
	lprintf(lvl-1, "|-----> [dev_id=%d] Buffer::allocate()\n", dev_id);
#endif
	for(int i=0; i<BlockNum; i++)
		if(Blocks[i]!=NULL) Blocks[i]->allocate(lockfree);
		else error("[dev_id=%d] Buffer::allocate() -> Blocks[%d] was NULL\n", dev_id, i);
#ifdef CDEBUG
	lprintf(lvl-1, "<-----| [dev_id=%d] Buffer::allocate()\n", dev_id);
#endif
}
#endif



CBlock_p Buffer::assign_Cblock(state start_state, bool lockfree){
	// Assigns a block from buffer to be used for memory.
	// State options are:
	// - INVALID: Raise error.
	// - NATIVE:
	// - EXCLUSIVE: Will add one writer.
	// - SHARABLE: Will add one reader.
	// - AVAILABLE: Will be initialized as AVAILABLE and someone might take it.

	short lvl = 2;
#ifdef CDEBUG
	lprintf(lvl-1, "|-----> [dev_id=%d] Buffer::assign_Cblock()\n", dev_id);
#endif
	CBlock_p result = NULL;
	if(!lockfree) lock(); // Lock buffer
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
		lprintf(lvl-1,"------ [dev_id=%d] Buffer::assign_Cblock(): Block with id=%d reseted.\n", dev_id, remove_block_idx);
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
	lprintf(lvl-1, "<-----| [dev_id=%d] Buffer::assign_Cblock()\n", dev_id);
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
int BufferSelectBlockToRemove_naive(Buffer_p buffer, bool lockfree){
	short lvl = 2;

	if (buffer == NULL)
		error("BufferSelectBlockToRemove_naive(): Called on empty buffer\n");

#ifdef CDEBUG
	lprintf(lvl-1, "|-----> [dev_id=%d] BufferSelectBlockToRemove_naive()\n",buffer->dev_id);
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
			lprintf(lvl-1, "------- [dev_id=%d] BufferSelectBlockToRemove_naive(): Found available block. Invalidated.\n",buffer->dev_id);
		#endif
			break;
		}
		if(!lockfree)
			buffer->Blocks[idx]->unlock();
	}
	#ifdef CDEBUG
		lprintf(lvl-1, "<-----| [dev_id=%d] BufferSelectBlockToRemove_naive()\n",buffer->dev_id);
	#endif
	return result_idx;
}

int BufferSelectExclusiveBlockToRemove_naive(Buffer_p buffer, bool lockfree){
	short lvl = 2;

	if (buffer == NULL)
		error("BufferSelectExclusiveBlockToRemove_naive(): Called on empty buffer\n");

#ifdef CDEBUG
	lprintf(lvl-1, "|-----> [dev_id=%d] BufferSelectExclusiveBlockToRemove_naive()\n",buffer->dev_id);
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
				lprintf(lvl-1, "------- [dev_id=%d] BufferSelectExclusiveBlockToRemove_naive(): Found exclusive block with no pernding operations on it. Invalidated.\n",buffer->dev_id);
			#endif
			}
			if(!lockfree)
				native_block->unlock();
		}
		if(!lockfree)
			buffer->Blocks[idx]->unlock();
	}
	#ifdef CDEBUG
		lprintf(lvl-1, "<-----| [dev_id=%d] BufferSelectExclusiveBlockToRemove_naive()\n",buffer->dev_id);
	#endif
	return result_idx;
}

#elif defined(FIFO) || defined(MRU) || defined(LRU)
Node_LL_p BufferSelectBlockToRemove_fifo_mru_lru(Buffer_p buffer, bool lockfree){
	short lvl = 2;
	if (buffer == NULL)
		error("BufferSelectBlockToRemove_fifo_mru_lru(): Called on empty buffer\n");

#ifdef CDEBUG
	lprintf(lvl-1, "|-----> [dev_id=%d] BufferSelectBlockToRemove_fifo_mru_lru()\n", buffer->dev_id);
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
				lprintf(lvl-1, "------- [dev_id=%d] BufferSelectBlockToRemove_fifo_mru_lru(): Found available block. Invalidated.\n",buffer->dev_id);
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
	lprintf(lvl-1, "<-----| [dev_id=%d] BufferSelectBlockToRemove_fifo_mru_lru()\n",buffer->dev_id);
#endif
	return result_node;
}

Node_LL_p BufferSelectExclusiveBlockToRemove_fifo_mru_lru(Buffer_p buffer, bool lockfree){
	short lvl = 2;
	if (buffer == NULL)
		error("BufferSelectExclusiveBlockToRemove_fifo_mru_lru(): Called on empty buffer\n");

#ifdef CDEBUG
	lprintf(lvl-1, "|-----> [dev_id=%d] BufferSelectExclusiveBlockToRemove_fifo_mru_lru()\n", buffer->dev_id);
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
					lprintf(lvl-1, "------- [dev_id=%d] BufferSelectExclusiveBlockToRemove_fifo_mru_lru(): Found exclusive block with no pernding operations on it. Invalidated.\n",buffer->dev_id);
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
	lprintf(lvl-1, "<-----| [dev_id=%d] BufferSelectExclusiveBlockToRemove_fifo_mru_lru()\n",buffer->dev_id);
#endif
	return result_node;
}
#endif


// ///
// /// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
// /// \author Theodoridis Aristomenis (theodoridisaristomenis@gmail.com)
// ///
// /// \brief The CHLpeLia caching functions.
// ///

// //#include <cassert>
// //#include <atomic>

// #include "unihelpers.hpp"
// #include "DataCaching.hpp"
// //#include "backend_wrappers.hpp"


// Cache_p Global_Cache[32] = {NULL};
// int CBlock_ctr[32] = {0};
// int DevCache_ctr = 0;

// int globalock = 0;

// #if defined(FIFO) || defined(MRU) || defined(LRU)

// /**************************
//  ** LinkedList Functions **
//  **************************/

// LinkedList::LinkedList(Cache_p cache, std::string name){
// 	short lvl = 2;
// 	if(cache==NULL)
// 		error("LinkedList::LinkedList(): Creating cache that doesn't belong to a cache.\n");
// #ifdef CDEBUG
// 	lprintf(lvl-1, "|-----> [dev_id=%d] LinkedList::LinkedList(name=%s):\n", cache->dev_id, name.c_str());
// #endif
// 	lock_ll = 0;
// 	Parent = cache;
// 	Name = name;
// 	start = NULL;
// 	end = NULL;
// 	length = 0;
// 	iter = NULL;
// #ifdef CDEBUG
// 	lprintf(lvl-1, "<-----| [dev_id=%d] LinkedList::LinkedList(name=%s)\n", Parent->dev_id, Name.c_str());
// #endif
// }

// LinkedList::~LinkedList(){
// 	short lvl = 2;
// #ifdef CDEBUG
// lprintf(lvl-1, "|-----> [dev_id=%d] LinkedList::~LinkedList(name=%s)\n", Parent->dev_id, Name.c_str());
// #endif
// 	lock();
// 	Node_LL_p tmp;
// 	tmp = start;
// 	for(int i=0; i<length && tmp!=end; i++){
// 		tmp = start;
// 		start = tmp->next;
// 		delete tmp;
// 	}
// 	unlock();
// #ifdef CDEBUG
// lprintf(lvl-1, "<-----| [dev_id=%d] LinkedList::~LinkedList(name=%s)\n", Parent->dev_id, Name.c_str());
// #endif
// }

// void LinkedList::draw_queue(bool lockfree){
// 	short lvl = 1;
// 	if(!lockfree)
// 		lock();
// 	int count = 1;
// #if defined(FIFO)
// 	lprintf(lvl-1, " FIFO");
// #elif defined(MRU)
// 	lprintf(lvl-1, " MRU");
// #elif defined(LRU)
// 	lprintf(lvl-1, " LRU");
// #endif
// 	lprintf(lvl-1, " Queue:\
// 						\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
// 						\n||      Cache Id       | %d\
// 						\n|| - - - - - - - - - - - - - - - - - - - -\
// 						\n||      Device Id      | %d\
// 						\n|| - - - - - - - - - - - - - - - - - - - -\
// 						\n||        Name         | %s\
// 						\n|| - - - - - - - - - - - - - - - - - - - -\
// 						\n||       Length        | %d\
// 						\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n", Parent->id, Parent->dev_id, Name.c_str(), length);
// 	if(length>0){
// 	iter = start;
// 	lprintf(lvl-1, " Position: %d\
// 						\n_________________________________________\
// 						\n|       Idx       | %d\
// 						\n| - - - - - - - - - - - - - - - - - - - -\
// 						\n|      Valid      | %s\
// 						\n|________________________________________\
// 						\n", count, iter->idx, iter->valid ? "True" : "False");
// 		iter = iter->next;
// 		count++;
// 		while(iter!=NULL && count<=length){
// 			lprintf(lvl-1, " Position: %d\
// 								\n_________________________________________\
// 								\n|       Idx       | %d\
// 								\n| - - - - - - - - - - - - - - - - - - - -\
// 								\n|      Valid      | %s\
// 								\n|________________________________________\
// 								\n", count, iter->idx, iter->valid ? "True" : "False");
// 			iter = iter->next;
// 			count++;
// 		}
// 	}
// 	if(!lockfree)
// 		unlock();
// }

// void LinkedList::invalidate(Node_LL_p node, bool lockfree){
// 	short lvl = 2;
// #ifdef CDEBUG
// 	lprintf(lvl-1, "|-----> [dev_id=%d] LinkedList::invalidate(name=%s, node_id=%d)\n", Parent->dev_id, Name.c_str(), node->idx);
// #endif
// 	if(node == NULL)
// 		error("[dev_id=%d] LinkedList::invalidate: Node not found.\n", Parent->dev_id);
// 	else if(node->idx>=Parent->BlockNum)
// 		error("[dev_id=%d] LinkedList::invalidate: Node idx (%d) is larger than the BlockNum(%d).\n",
// 				Parent->dev_id,  node->idx, Parent->BlockNum);
// 	else{
// 		if(!lockfree)
// 			lock();
// 		node->valid = false;
// 		if(!lockfree)
// 			unlock();
// 		}
// #ifdef CDEBUG
// 	lprintf(lvl-1, "<-----| [dev_id=%d] LinkedList::invalidate(name=%s, node_id=%d)\n", Parent->dev_id, Name.c_str(), node->idx);
// #endif
// }

// void LinkedList::push_back(int idx, bool lockfree){
// 	// Pushes an element in the back of the queue.
// 	short lvl = 2;
// #ifdef CDEBUG
// 	lprintf(lvl-1, "|-----> [dev_id=%d] LinkedList::push_back(name=%s, idx=%d)\n", Parent->dev_id, Name.c_str(), idx);
// #endif
// 	if(!lockfree)
// 		lock();
// 	if(idx >= Parent->BlockNum)
// 		error("[dev_id=%d] LinkedList::push_back(name=%s): Index given(%d) is larger than the number of blocks(%d).\n", Parent->dev_id, Name.c_str(), idx, Parent->BlockNum);
// 	else if(idx < 0)
// 		error("[dev_id=%d] LinkedList::push_back(name=%s): Index given(%d) is not valid.\n", Parent->dev_id, Name.c_str(), idx);
// 	else if(length > Parent->BlockNum)
// 		error("[dev_id=%d] LinkedList::push_back(name=%s): Called to put another element but max length is reached with BlockNum=%d and length=%d\n", Parent->dev_id, Name.c_str(), Parent->BlockNum, length);
// 	Node_LL_p tmp = new Node_LL();
// 	tmp->next = NULL;
// 	tmp->previous = NULL;
// 	tmp->idx = idx;
// 	if(start == NULL){ // Queue is empty.
// 		start = tmp;
// 		end = tmp;
// 	}
// 	else{ // Queue not empty.
// 		end->next = tmp;
// 		tmp->previous = end;
// 		end = tmp;
// 	}
// 	length++;
// 	if(!lockfree)
// 		unlock();
// #ifdef CDEBUG
// 	lprintf(lvl-1, "<-----| [dev_id=%d] LinkedList::push_back(name=%s, idx=%d)\n", Parent->dev_id, Name.c_str(), idx);
// #endif

// }


// Node_LL_p LinkedList::start_iterration(){
// 	// Returns the first valid element without removing it.
// 	short lvl = 2;
// #ifdef CDEBUG
// 	lprintf(lvl-1, "|-----> [dev_id=%d] LinkedList::start_iterration(name=%s)\n", Parent->dev_id, Name.c_str());
// #endif
// 	iter = start;
// 	if(iter == NULL){
// 		Node_LL_p tmp_node = new Node_LL();
// 		tmp_node->idx = -1;
// 	#ifdef CDEBUG
// 		lprintf(lvl-1, "<-----| [dev_id=%d] LinkedList::start_iterration(name=%s, node_idx=%d)\n", Parent->dev_id, Name.c_str(), tmp_node->idx);
// 	#endif
// 		return tmp_node;
// 	}
// 	else{
// 	#ifdef CDEBUG
// 		lprintf(lvl-1, "<-----| [dev_id=%d] LinkedList::start_iterration(name=%s, node_idx=%d)\n", Parent->dev_id, Name.c_str(), iter->idx);
// 	#endif
// 		return iter;
// 	}
// }

// Node_LL_p LinkedList::next_in_line(){
// 	// Returns next element in iterration.
// 	short lvl = 2;
// #ifdef CDEBUG
// 	lprintf(lvl-1, "|-----> [dev_id=%d] LinkedList::next_in_line(name=%s)\n", Parent->dev_id, Name.c_str());
// #endif
// 	if(iter == NULL)
// 		error("[dev_id=%d] LinkedList::next_in_line(name=%s): Iterration not started. Call check_first.\n", Parent->dev_id, Name.c_str());
// 	else{
// 		iter = iter->next;
// 		if(iter == NULL){
// 			Node_LL_p tmp_node = new Node_LL();
// 			tmp_node->idx = -1;
// 		#ifdef CDEBUG
// 			lprintf(lvl-1, "<-----| [dev_id=%d] LinkedList::next_in_line(name=%s, node_idx=%d)\n", Parent->dev_id, Name.c_str(), tmp_node->idx);
// 		#endif
// 			return tmp_node;
// 		}
// 		else{
// 		#ifdef CDEBUG
// 			lprintf(lvl-1, "<-----| [dev_id=%d] LinkedList::next_in_line(name=%s, node_idx=%d)\n", Parent->dev_id, Name.c_str(), iter->idx);
// 		#endif
// 			return iter;
// 		}
// 	}
// }

// Node_LL_p LinkedList::remove(Node_LL_p node, bool lockfree){
// 	// Removes first element.
// 	short lvl = 2;

// #ifdef CDEBUG
// 	if(node!=NULL)
// 		lprintf(lvl-1, "|-----> [dev_id=%d] LinkedList::remove(name=%s, idx=%d)\n", Parent->dev_id, Name.c_str(), node->idx);
// #endif
// 	if(node == NULL) // Node not found.
// 		error("[dev_id=%d] LinkedList::remove(name=%s): Input node not found. Can't pop element.\n", Parent->dev_id, Name.c_str());
// 	else if(node->idx >= Parent->BlockNum)
// 		error("[dev_id=%d] LinkedList::remove(name=%s): Index of given node(%d) is larger than the number of blocks(%d).\n", Parent->dev_id, Name.c_str(), node->idx, Parent->BlockNum);
// 	else if(node->idx < 0)
// 		error("[dev_id=%d] LinkedList::remove(name=%s): Index of given node(%d) is negative.\n", Parent->dev_id, Name.c_str(), node->idx);
// 	else{
// 		if(!lockfree)
// 			lock();
// 		if(length == 0)
// 			error("[dev_id=%d] LinkedList::remove(name=%s): Queue empty.\n", Parent->dev_id, Name.c_str());
// 		else if(length == 1){
// 			if(start == node){
// 				start = NULL;
// 				end = NULL;
// 			}
// 			else
// 				error("[dev_id=%d] LinkedList::remove(name=%s): Node to be removed does not belong to the list.\n", Parent->dev_id, Name.c_str());
// 		}
// 		else if(node == start){
// 			start = node->next;
// 			if(start != NULL)
// 				start->previous = NULL;
// 		}
// 		else if(node == end){
// 			end = node->previous;
// 			if(end != NULL)
// 				end->next = NULL;
// 		}
// 		else{
// 			(node->previous)->next = node->next;
// 			(node->next)->previous = node->previous;
// 		}
// 		node->next = NULL;
// 		node->previous = NULL;
// 		length--;
// 		if(!lockfree)
// 			unlock();
// #ifdef CDEBUG
// lprintf(lvl-1, "<-----| [dev_id=%d] LinkedList::remove(name=%s, idx=%d)\n", Parent->dev_id, Name.c_str(), node->idx);
// #endif
// 		return node;
// 	}
// }

// void LinkedList::put_first(Node_LL* node, bool lockfree){
// 	// Puts the element first on the list.
// 	short lvl = 2;

// #ifdef CDEBUG
// lprintf(lvl-1, "|-----> [dev_id=%d] LinkedList::put_first(name=%s, idx=%d)\n", Parent->dev_id, Name.c_str(), node->idx);
// #endif

// 	if(node == NULL)
// 		error("[dev_id=%d] LinkedList::put_first(name=%s): Called to put a node that doesn't exist.\n", Parent->dev_id, Name.c_str());
// 	else if(node->idx < 0)
// 		error("[dev_id=%d] LinkedList::put_first(name=%s): Called to put an invalid node .\n", Parent->dev_id, Name.c_str());
// 	else{
// 		if(!lockfree){
// 			lock();
// 		}
// 		// Add it to the new queue
// 		if(length == 0){
// 			end = node;
// 			start = node;
// 			node->next = NULL;
// 			node->previous = NULL;
// 		}
// 		else{
// 			node->next = start;
// 			start->previous = node;
// 			start = node;
// 			node->previous = NULL;
// 		}
// 		length++;
// 		if(!lockfree){
// 			unlock();
// 		}
// 	}
// #ifdef CDEBUG
// lprintf(lvl-1, "<-----| [dev_id=%d] LinkedList::put_first(name=%s, idx=%d)\n", Parent->dev_id, Name.c_str(), node->idx);
// #endif
// }

// void LinkedList::put_last(Node_LL_p node, bool lockfree){
// 	// Takes a node and puts it in the end of the list.
// 	short lvl = 2;

// #ifdef CDEBUG
// lprintf(lvl-1, "|-----> [dev_id=%d] LinkedList::put_last(name=%s, idx=%d)\n", Parent->dev_id, Name.c_str(), node->idx);
// #endif

// 	if(node == NULL)
// 		error("[dev_id=%d] LinkedList::put_last(name=%s): Called to put a node that doesn't exist.\n", Parent->dev_id, Name.c_str());
// 	else if(node->idx < 0)
// 		error("[dev_id=%d] LinkedList::put_last(name=%s): Called to put an invalid node .\n", Parent->dev_id, Name.c_str());
// 	else{
// 		if(!lockfree){
// 			lock();
// 		}
// 		// Add it to the new queue
// 		if(length == 0){
// 			end = node;
// 			start = node;
// 			node->next = NULL;
// 			node->previous = NULL;
// 		}
// 		else{
// 			node->previous = end;
// 			end->next = node;
// 			end = node;
// 			node->next = NULL;
// 		}
// 		length++;
// 		if(!lockfree){
// 			unlock();
// 		}
// 	}
// #ifdef CDEBUG
// lprintf(lvl-1, "<-----| [dev_id=%d] LinkedList::put_last(name=%s, idx=%d)\n", Parent->dev_id, Name.c_str(), node->idx);
// #endif
// }


// void LinkedList::lock(){
// 	// lprintf(0, " |------> [dev_id=%d] Locking %s\n", Parent->dev_id, Name.c_str());
// 	// lprintf(0, "|>>>>>> [dev_id=%d] list \n", Parent->dev_id);
// 	while(__sync_lock_test_and_set(&lock_ll, 1));
// 	// lprintf(0, "<<<<<<| [dev_id=%d] list \n", Parent->dev_id);
// 	// lprintf(0, " <------| [dev_id=%d] Locked %s\n", Parent->dev_id, Name.c_str());
// 	// Lock++;
// 	// Lock.lock();
// }

// void LinkedList::unlock(){
// 	// lprintf(0, " |------> [dev_id=%d] Unlocking %s\n", Parent->dev_id, Name.c_str());
//    __sync_lock_release(&lock_ll);
// 	// lprintf(0, " <------| [dev_id=%d] Unlocked %s\n", Parent->dev_id, Name.c_str());
// 	// Lock--;
// }

// bool LinkedList::is_locked(){
// 	if(lock_ll==0)
// 		return false;
// 	return true;
// }

// #endif

// const char* print_state(state in_state){
// 	switch(in_state){
// 		case(INVALID):
// 			return "INVALID";
// 		case(NATIVE):
// 			return "NATIVE";
// 		case(EXCLUSIVE):
// 			return "EXCLUSIVE";
// 		case(SHARABLE):
// 			return "SHARABLE";
// 		case(AVAILABLE):
// 			return "AVAILABLE";
// 		default:
// 			error("print_state: Unknown state\n");
// 	}
// 	return "";
// }

// /***************************
//  ** Cache Block Functions **
//  ***************************/

// CacheBlock::CacheBlock(int block_id, Cache_p block_parent, long long block_size){
// 	// Constructor for Blocks.
// 	// Args:
// 	// - block_id: Id of the block.
// 	// - block_parent: Cache that the block belongs to.
// 	// - block_size: Size of usable memory in the block.

// 	short lvl = 2;

// 	if(block_parent!=NULL && block_id>=0 && block_id<block_parent->BlockNum && block_size>0){
// 	#ifdef CDEBUG
// 		lprintf(lvl-1, "|-----> [dev_id=%d] CacheBlock::CacheBlock(id=%d, Parent_id=%d, Size=%llu)\n", block_parent->dev_id, block_id, block_parent->id, block_size);
// 	#endif
// 		Lock = 0;
// 		id = block_id;
// 		// TODO: Name?
// 		Owner_p = NULL;
// 		WritebackData_p = NULL;
// 		Parent = block_parent;
// 		Size = block_size;
// 		PendingReaders = 0;
// 		PendingWriters = 0;

// 		Adrs = NULL; // Will be set by cache allocate.
// 		State = INVALID; 	//	Technically no data in, so INVALID?
// 												//	But AVAILABLE correct for scheduling out as well...
// 		Available = new Event(Parent->dev_id);
// 	#ifdef CDEBUG
// 		lprintf(lvl-1, "<-----| [dev_id=%d] CacheBlock::CacheBlock()\n", Parent->dev_id);
// 	#endif
// 	}
// 	else{
// 		if(block_parent==NULL)
// 			error("CacheBlock::CacheBlock(): Constructor called with no cache to belong.\n");
// 		else if(block_id<0 && block_id>=block_parent->BlockNum)
// 			error("[dev_id=%d] CacheBlock::CacheBlock(): Constructor called with invalid id=%d.\n", block_parent->dev_id, block_id);
// 		else
// 			error("[dev_id=%d] CacheBlock::CacheBlock(): Constructor called with invalid mem size=%llu.\n", block_parent->dev_id, block_size);
// 	}
// }

// CacheBlock::~CacheBlock(){
// 	short lvl = 2;
// 	// Destructor of the block.
// #ifdef CDEBUG
// lprintf(lvl-1, "|-----> [dev_id=%d] CacheBlock::~CacheBlock()\n", Parent->dev_id);
// #endif
// 	lock();
// 	//reset(true, true);
// 	if(Owner_p){
// 		*Owner_p = NULL;
// 		Owner_p = NULL;
// 	}
// 	free(WritebackData_p);
// 	if(State != NATIVE){
// 		CHLFree(Adrs, Parent->dev_id);
// 		delete Available;
// #ifdef CDEBUG
// 		lprintf(lvl-1, "------- [dev_id=%d] CacheBlock::~CacheBlock(): Deleting non-NATIVE block id =%d\n",
// 			Parent->dev_id, id);
// #endif
// 	}
// 	else{;
// #ifdef CDEBUG
// 		lprintf(lvl-1, "------- [dev_id=%d] CacheBlock::~CacheBlock(): Refrain from deleting NATIVE block id =%d\n",
// 			Parent->dev_id, id);
// #endif
// 	}
// 	unlock();
// #ifdef CDEBUG
// lprintf(lvl-1, "<-----| [dev_id=%d] CacheBlock::~CacheBlock()\n", Parent->dev_id);
// #endif
// }

// void CacheBlock::draw_block(bool lockfree){
// 	// Draws the block for debugging purposes.
// 	short lvl=0;

// 	if(!lockfree)
// 		lock();
// 	lprintf(lvl-1, " Block:   \
// 						\n_________________________________________\
// 						\n|       Id        | %d\
// 						\n| - - - - - - - - - - - - - - - - - - - -\
// 						\n|      Name       | %s\
// 						\n| - - - - - - - - - - - - - - - - - - - -\
// 						\n|      Size       | %llu\
// 						\n| - - - - - - - - - - - - - - - - - - - -\
// 						\n|      State      | %s \
// 						\n| - - - - - - - - - - - - - - - - - - - -\
// 						\n| Pending Readers | %d\
// 						\n| - - - - - - - - - - - - - - - - - - - -\
// 						\n| Pending Writers | %d\
// 						\n|________________________________________\
// 						\n", id, Name.c_str(), Size, print_state(State), PendingReaders.load(), PendingWriters.load());
// 	if(!lockfree)
// 		unlock();
// }

// void CacheBlock::add_reader(bool lockfree){
// 	short lvl = 2;
// #ifdef CDEBUG
// 	lprintf(lvl-1, "|-----> [dev_id=%d] CacheBlock::add_reader(block_id=%d)\n", Parent->dev_id, id);
// #endif
// 	if(!lockfree){
// #if defined(FIFO) || defined(MRU) || defined(LRU)
// 		Parent->InvalidQueue->lock();
// 		Parent->Queue->lock();
// #endif
// 		lock();
// 	}
// 	if(State==SHARABLE || State==EXCLUSIVE || State == NATIVE)
// 		PendingReaders++;
// 	else if(State==AVAILABLE){
// 		PendingReaders++;
// 		set_state(SHARABLE, true);
// 	}
// 	else
// 		error("[dev_id=%d] CacheBlock::add_reader(): Can't add reader. Block has State=%s\n", Parent->dev_id, print_state(State));
// 	if(!lockfree){
// 		unlock();
// #if defined(FIFO) || defined(MRU) || defined(LRU)
// 		Parent->InvalidQueue->unlock();
// 		Parent->Queue->unlock();
// #endif
// 	}
// #ifdef CDEBUG
// 	lprintf(lvl-1, "<-----| [dev_id=%d] CacheBlock::add_reader(block_id=%d)\n", Parent->dev_id, id);
// #endif
// }

// void CacheBlock::add_writer(bool lockfree){
// 	short lvl = 2;
// 	if(!lockfree){
// #if defined(FIFO) || defined(MRU) || defined(LRU)
// 		Parent->InvalidQueue->lock();
// 		Parent->Queue->lock();
// #endif
// 		lock();
// 	}
// #ifdef CDEBUG
// 	lprintf(lvl-1, "|-----> [dev_id=%d] CacheBlock::add_writer(block_id=%d)\n", Parent->dev_id, id);
// #endif
// 	if(State==EXCLUSIVE || State==NATIVE)
// 		PendingWriters++;
// 	else if(State==AVAILABLE || State==SHARABLE){
// 		lprintf(lvl-1, "------- [dev_id=%d] CacheBlock::add_writer(block_id=%d): Making exclusive....\n", Parent->dev_id, id);
// 		set_state(EXCLUSIVE, true);
// 		PendingWriters++;
// 	}
// 	else
// 		error("[dev_id=%d] CacheBlock::add_reader(): Can't add reader. Block has State=%s\n", Parent->dev_id, print_state(State));
// 	if(!lockfree){
// 		unlock();
// #if defined(FIFO) || defined(MRU) || defined(LRU)
// 		Parent->InvalidQueue->unlock();
// 		Parent->Queue->unlock();
// #endif
// 	}
// #ifdef CDEBUG
// 	lprintf(lvl-1, "<-----| [dev_id=%d] CacheBlock::add_writer(block_id=%d)\n", Parent->dev_id, id);
// #endif
// }

// void CacheBlock::remove_reader(bool lockfree){
// 	short lvl = 2;
// #ifdef CDEBUG
// lprintf(lvl-1, "|-----> [dev_id=%d] CacheBlock::remove_reader(block_id=%d)\n", Parent->dev_id, id);
// #endif
// 	if(!lockfree){
// #if defined(FIFO) || defined(MRU) || defined(LRU)
// 		Parent->InvalidQueue->lock();
// 		Parent->Queue->lock();
// #endif
// 		lock();
// 	}

// 	if(PendingReaders.load()>0)
// 		PendingReaders--;
// 	else
// 		error("[dev_id=%d] CacheBlock::remove_reader(): Can't remove reader. There are none.\n", Parent->dev_id);
// 	update_state(true);
// #if defined(MRU)
// 	Node_LL_p node = Parent->Queue->remove(Parent->Hash[id], true);
// 	Parent->Queue->put_first(node, true);
// #elif defined(LRU)
// 	Node_LL_p node = Parent->Queue->remove(Parent->Hash[id], true);
// 	Parent->Queue->put_last(node, true);
// #endif
// 	if(!lockfree){
// 		unlock();
// #if defined(FIFO) || defined(MRU) || defined(LRU)
// 		Parent->Queue->unlock();
// 		Parent->InvalidQueue->unlock();
// #endif
// 	}
// #ifdef CDEBUG
// lprintf(lvl-1, "<-----| [dev_id=%d] CacheBlock::remove_reader(block_id=%d)\n", Parent->dev_id, id);
// #endif
// }

// void CacheBlock::remove_writer(bool lockfree){
// 	short lvl = 2;
// #ifdef CDEBUG
// 	lprintf(lvl-1, "|-----> [dev_id=%d] CacheBlock::remove_writer(block_id=%d)\n", Parent->dev_id, id);
// #endif
// 	if(!lockfree){
// #if defined(FIFO) || defined(MRU) || defined(LRU)
// 		Parent->InvalidQueue->lock();
// 		Parent->Queue->lock();
// #endif
// 		lock();
// 	}
// 	if(PendingWriters.load()>0)
// 		PendingWriters--;
// 	else
// 		error("[dev_id=%d] CacheBlock::remove_writer(block_id=%d): Can't remove writer. There are none.\n", Parent->dev_id, id);
// 	update_state(true);
// #if defined(MRU)
// 	Node_LL_p node = Parent->Queue->remove(Parent->Hash[id], true);
// 	Parent->Queue->put_first(node, true);
// #elif defined(LRU)
// 	Node_LL_p node = Parent->Queue->remove(Parent->Hash[id], true);
// 	Parent->Queue->put_last(node, true);
// #endif
// 	if(!lockfree){
// 		unlock();
// #if defined(FIFO) || defined(MRU) || defined(LRU)
// 		Parent->Queue->unlock();
// 		Parent->InvalidQueue->unlock();
// #endif
// 	}
// #ifdef CDEBUG
// 	lprintf(lvl-1, "<-----| [dev_id=%d] CacheBlock::remove_writer(block_id=%d)\n", Parent->dev_id, id);
// #endif
// }

// void* CBlock_RR_wrap(void* CBlock_wraped){
// 	CBlock_wrap_p CBlock_unwraped = (CBlock_wrap_p) CBlock_wraped;
// 	CBlock_unwraped->CBlock->remove_reader(CBlock_unwraped->lockfree);
// 	return NULL;
// }

// void* CBlock_RW_wrap(void* CBlock_wraped){
// 	CBlock_wrap_p CBlock_unwraped = (CBlock_wrap_p) CBlock_wraped;
// 	//if(!CBlock_unwraped->CBlock->PendingWriters.load())
// 	//	printf("|-----> CBlock_RW_wrap: suspicious PendingWriters = %d\n", CBlock_unwraped->CBlock->PendingWriters.load());
// 	CBlock_unwraped->CBlock->remove_writer(CBlock_unwraped->lockfree);
// 	return NULL;
// }

// void* CBlock_INV_wrap(void* CBlock_wraped){
// 	CBlock_wrap_p CBlock_unwraped = (CBlock_wrap_p) CBlock_wraped;
// 	CBlock_unwraped->CBlock->Available->soft_reset();
// 	CBlock_unwraped->CBlock->set_state(INVALID, CBlock_unwraped->lockfree);
// 	return NULL;
// }

// void* CBlock_RR_INV_wrap(void* CBlock_wraped){
// 	CBlock_wrap_p CBlock_unwraped = (CBlock_wrap_p) CBlock_wraped;
// 	if(!CBlock_unwraped->lockfree) CBlock_unwraped->CBlock->lock();
// 	CBlock_unwraped->CBlock->remove_reader(true);
// 	CBlock_unwraped->CBlock->Available->soft_reset();
// 	CBlock_unwraped->CBlock->set_state(INVALID, true);
// 	if(!CBlock_unwraped->lockfree) CBlock_unwraped->CBlock->unlock();
// 	return NULL;
// }

// void* CBlock_RW_INV_wrap(void* CBlock_wraped){
// 	CBlock_wrap_p CBlock_unwraped = (CBlock_wrap_p) CBlock_wraped;
// 	if(!CBlock_unwraped->lockfree) CBlock_unwraped->CBlock->lock();
// 	//if(!CBlock_unwraped->CBlock->PendingWriters.load())
// 	//	printf("|-----> CBlock_RW_INV_wrap: suspicious PendingWriters = %d\n", CBlock_unwraped->CBlock->PendingWriters.load());
// 	CBlock_unwraped->CBlock->remove_writer(true);
// 	CBlock_unwraped->CBlock->Available->soft_reset();
// 	CBlock_unwraped->CBlock->set_state(INVALID, true);
// 	if(!CBlock_unwraped->lockfree) CBlock_unwraped->CBlock->unlock();
// 	return NULL;
// }

// void CacheBlock::set_owner(void** owner_adrs, bool lockfree){
// 	short lvl = 2;
// 	#ifdef CDEBUG
// 		lprintf(lvl-1, "|-----> CacheBlock::set_owner(owner_adrs=%p)\n", owner_adrs);
// 	#endif

// 	if(!lockfree)
// 		lock();
// 	Owner_p = owner_adrs;
// 	if(!lockfree)
// 		unlock();
// #ifdef CDEBUG
// 	lprintf(lvl-1, "<-----| CacheBlock::set_owner()\n");
// #endif
// }

// void CacheBlock::reset(bool lockfree, bool forceReset){
// 	// Resets block attibutes if it's AVAILABLE to be used again.
// 	short lvl = 2;
// #ifdef CDEBUG
// 	lprintf(lvl-1, "|-----> [dev_id=%d] CacheBlock::reset(block_id=%d)\n", Parent->dev_id, id);
// #endif
// 	if(State==INVALID || State==AVAILABLE || forceReset){
// 		if(!lockfree){
// 		#if defined(FIFO) || defined(MRU) || defined(LRU)
// 			Parent->InvalidQueue->lock();
// 			Parent->Queue->lock();
// 		#endif
// 			lock();
// 		}
// 		PendingReaders = 0;
// 		PendingWriters = 0;
// 		free(WritebackData_p);
// 		WritebackData_p = NULL;
// 		Available->reset();

// 		if(forceReset && State==NATIVE){
// 			Adrs = NULL;
// 			State = INVALID;
// 		#if defined(FIFO) || defined(MRU) || defined(LRU)
// 			Node_LL_p node = Parent->Queue->remove(Parent->Hash[id], true);
// 			Parent->InvalidQueue->put_last(node, true);
// 			node->valid=false;
// 		#endif
// 		}
// 		else{
// 			set_state(INVALID, true);
// 		}

// 		if(Owner_p){
// 			*Owner_p = NULL;
// 			Owner_p = NULL;
// 		}
// 		if(!lockfree){
// 			unlock();
// 		#if defined(FIFO) || defined(MRU) || defined(LRU)
// 			Parent->InvalidQueue->unlock();
// 			Parent->Queue->unlock();
// 		#endif
// 		}
// 	#ifdef CDEBUG
// 		if(forceReset)
// 			lprintf(lvl-1, "------- [dev_id=%d] CacheBlock::reset(): Block with id=%d forced to be reseted.\n", Parent->dev_id, id);
// 		else
// 			lprintf(lvl-1, "------- [dev_id=%d] CacheBlock::reset(): Block with id=%d reseted.\n", Parent->dev_id, id);
// 	#endif
// 	}
// 	else
// 		error("[dev_id=%d] CacheBlock::reset(): Reset was called on a %s block.\n", Parent->dev_id, print_state(State));
// #ifdef CDEBUG
// 	lprintf(lvl-1, "<-----| [dev_id=%d] CacheBlock::reset(block_id=%d)\n", Parent->dev_id, id);
// #endif
// }

// void CacheBlock::init_writeback_info(CBlock_p WB_block, int* RW_master_p,
// 	int dim1, int dim2, int ldim, int ldim_wb, int dtype_sz, CQueue_p wb_queue, bool lockfree){
// 	short lvl = 2;
// #ifdef CDEBUG
// 	lprintf(lvl-1, "|-----> [dev_id=%d] CacheBlock::init_writeback_info(block_id=%d)\n", Parent->dev_id, id);
// #endif

// 	if(WritebackData_p!=NULL) error("[dev_id=%d] CacheBlock::init_writeback_info(block_id=%d):\
// 		Called with a previously allocated WritebackData_p\n", Parent->dev_id, id);

// 	if(!lockfree){
// #if defined(FIFO) || defined(MRU) || defined(LRU)
// 		Parent->InvalidQueue->lock();
// 		Parent->Queue->lock();
// #endif
// 		lock();
// 	}

// 	WritebackData_p = (writeback_info_p) malloc(sizeof(struct writeback_info));
// 	WritebackData_p->Native_block = WB_block;
// 	WritebackData_p->WB_master_p = RW_master_p;
// 	WritebackData_p->dim1 = dim1;
// 	WritebackData_p->dim2 = dim2;
// 	WritebackData_p->ldim = ldim;
// 	WritebackData_p->ldim_wb = ldim_wb;
// 	WritebackData_p->dtype_sz = dtype_sz;
// 	WritebackData_p->wb_queue = wb_queue;
// 	if(!lockfree){
// #if defined(FIFO) || defined(MRU) || defined(LRU)
// 		Parent->InvalidQueue->unlock();
// 		Parent->Queue->unlock();
// #endif
// 		unlock();
// 	}
// 	#ifdef CDEBUG
// 	lprintf(lvl-1, "<-----| [dev_id=%d] CacheBlock::init_writeback_info(block_id=%d)\n", Parent->dev_id, id);
// 	#endif
// }

// void CacheBlock::write_back(bool lockfree){
// 	short lvl = 2;

// 	#ifdef CDEBUG
// 	lprintf(lvl-1, "|-----> [dev_id=%d] CacheBlock::write_back(block_id=%d)\n", Parent->dev_id, id);
// 	#endif

// 	if(WritebackData_p->Native_block==NULL)
// 	error("[dev_id=%d] CacheBlock::write_back(block_id=%d): Can't write back. Native block is NULL.\n", Parent->dev_id, id);
// 	else{
// 		if(!lockfree){
// 			#if defined(FIFO) || defined(MRU) || defined(LRU)
// 			Parent->InvalidQueue->lock();
// 			Parent->Queue->lock();
// 			#endif
// 			lock();
// 			WritebackData_p->Native_block->lock();
// 		}
// 		CBlock_p Write_back_Native_block = WritebackData_p->Native_block;
// 		/// We always lock for now, since we can't lock externally (since reset also resets WritebackData_p)
// 		CHLMemcpy2DAsync(WritebackData_p->Native_block->Adrs, WritebackData_p->ldim_wb, Adrs, WritebackData_p->ldim,
// 			WritebackData_p->dim1, WritebackData_p->dim2, WritebackData_p->dtype_sz,
// 			WritebackData_p->Native_block->Parent->dev_id, Parent->dev_id, WritebackData_p->wb_queue);
// 		*(WritebackData_p->WB_master_p) = WritebackData_p->Native_block->Parent->dev_id;
// 		if(!lockfree){
// 			add_reader(true);
// 			Write_back_Native_block->add_writer(true);
// 			Write_back_Native_block->unlock();
// 			unlock();
// 			#if defined(FIFO) || defined(MRU) || defined(LRU)
// 			Parent->Queue->unlock();
// 			Parent->InvalidQueue->unlock();
// 			#endif
// 		}
// 		WritebackData_p->wb_queue->sync_barrier();
// #ifdef CDEBUG
// 		lprintf(lvl-1, "------- [dev_id=%d] CacheBlock::write_back(block_id=%d): Writeback complete\n", Parent->dev_id, id);
// #endif
// 		if(!lockfree){
// 			Write_back_Native_block->remove_writer();
// 		}
// 		// One command, no need to lock here.
// 		// reset(lockfree, true);
// #ifdef CDEBUG
// 		lprintf(lvl-1, "------- [dev_id=%d] CacheBlock::write_back(block_id=%d): Reset block complete\n", Parent->dev_id, id);
// #endif
// 	}

// 	#ifdef CDEBUG
// 	lprintf(lvl-1, "<-----| [dev_id=%d] CacheBlock::write_back(block_id=%d)\n", Parent->dev_id, id);
// 	#endif
// }

// void CacheBlock::allocate(bool lockfree){
// 	// Allocates a cache block if not already pointing to some memory (not null!)
// 	short lvl = 2;
// #ifdef CDEBUG
// 	lprintf(lvl, "|-----> [dev_id=%d] CacheBlock::allocate(block_id=%d)\n", Parent->dev_id, id);
// #endif
// 	if(Adrs == NULL){
// 		Adrs = CHLMalloc(Size, Parent->dev_id);
// #ifdef CDEBUG
// 		lprintf(lvl, "------- [dev_id=%d] CacheBlock::allocate(block_id=%d): Allocated Adrs = %p\n", Parent->dev_id, id, Adrs);
// #endif
// 	}
// 	else{
// 		#ifdef CDEBUG
// 			lprintf(lvl, "------- [dev_id=%d] CacheBlock::allocate(block_id=%d) -> Supposedly already allocated block...", Parent->dev_id, id);
// 		#endif
// 	}
// #ifdef CDEBUG
// 	lprintf(lvl, "<-----| [dev_id=%d] CacheBlock::allocate(block_id=%d)\n", Parent->dev_id, id);
// #endif
// }

// state CacheBlock::get_state(){
// 	if(id < 0 || id >= Parent->BlockNum)
// 		error("[dev_id=%d] CacheBlock::get_state(): Invalid block id=%d\n", Parent->dev_id, id);
// 	return State;
// }

// state CacheBlock::set_state(state new_state, bool lockfree){
// 	// Forces a new state.
// 	short lvl = 2;
// #ifdef CDEBUG
// 	lprintf(lvl-1, "|-----> [dev_id=%d] CacheBlock::set_state(block_id=%d, prior_state=%s)\n", Parent->dev_id, id, print_state(State));
// #endif
// 	if(id < 0 || id >= Parent->BlockNum)
// 		error("[dev_id=%d] CacheBlock::set_state(%s): Invalid block id=%d\n", Parent->dev_id, print_state(new_state), id);
// 	if(!lockfree){
// 	#if defined(FIFO) || defined(MRU) || defined(LRU)
// 		Parent->InvalidQueue->lock();
// 		Parent->Queue->lock();
// 	#endif
// 		lock();
// 	}
// 	state old_state = State;
// 	if(old_state == NATIVE){;
// #ifdef CDEBUG
// 	lprintf(lvl-1, "<-----| [dev_id=%d] CacheBlock::set_state(block_id=%d, new_state=%s):\
// 	Tried to set state of NATIVE block, ignoring...\n", Parent->dev_id, id, print_state(new_state));
// #endif
// 	}
// 	else{
// 		State = new_state;
// 	#if defined(FIFO) || defined(MRU) || defined(LRU)
// 		if(State == INVALID && old_state != INVALID){
// 			Node_LL_p node = Parent->Queue->remove(Parent->Hash[id], true);
// 			Parent->InvalidQueue->put_last(node, true);
// 			node->valid=false;
// 		}
// 	#endif
// 	}
// 	if(!lockfree){
// 		unlock();
// 	#if defined(FIFO) || defined(MRU) || defined(LRU)
// 		Parent->Queue->unlock();
// 		Parent->InvalidQueue->unlock();
// 	#endif
// 	}

// #ifdef CDEBUG
// 	lprintf(lvl-1, "<-----| [dev_id=%d] CacheBlock::set_state(block_id=%d, new_state=%s)\n", Parent->dev_id, id, print_state(State));
// #endif
// 	return old_state;
// }

// int CacheBlock::update_state(bool lockfree){
// 	// Updates the state of the block. It cannot raise the state but only lower.
// 	short lvl = 2;
// #ifdef CDEBUG
// lprintf(lvl-1, "|-----> [dev_id=%d] CacheBlock::update_state(block_id=%d)\n", Parent->dev_id, id);
// #endif
// 	int ret = 0;
// 	if(!lockfree){
// 	#if defined(FIFO) || defined(MRU) || defined(LRU)
// 		Parent->InvalidQueue->lock();
// 		Parent->Queue->lock();
// 	#endif
// 		lock();
// 	}
// 	state prev_state = State;
// 	if(PendingWriters > 0){
// 		if(State!=EXCLUSIVE && State!=NATIVE)
// 			error("[dev_id=%d] CacheBlock::update_state(): Block has writers but state was %s.\n", Parent->dev_id, print_state(State));
// 	}
// 	else if(PendingReaders > 0){
// 		if(State==EXCLUSIVE){
// 			; // Do nothing - Not allowed to schedule out EXCLUSIVE blocks, unless we implement a writeback-to-native mechanism
// 		}
// 		else if(State!=SHARABLE && State!=NATIVE)
// 			error("[dev_id=%d] CacheBlock::update_state(): Block has readers but state was %s.\n", Parent->dev_id, print_state(State));
// 	}
// 	else if(State == SHARABLE){
// 		set_state(AVAILABLE, true);
// 		ret = 1;
// 	}
// 	else if(State == EXCLUSIVE){
// 		; // Do nothing - Not allowed to schedule out EXLUSIVE blocks, unless we implement a writeback-to-native mechanism
// 	}
// #ifdef CDEBUG
// 	if(ret==1)
// 		lprintf(lvl-1, "------- [dev_id=%d] CacheBlock::update_state(block_id=%d): Block state was changed %s -> %s \n", Parent->dev_id, id, print_state(prev_state), print_state(State));
// 	else
// 		lprintf(lvl-1, "------- [dev_id=%d] CacheBlock::update_state(block_id=%d): Block state is still %s \n", Parent->dev_id, id, print_state(State));
// #endif
// 	if(!lockfree){
// 		unlock();
// 	#if defined(FIFO) || defined(MRU) || defined(LRU)
// 		Parent->Queue->unlock();
// 		Parent->InvalidQueue->unlock();
// 	#endif
// 	}
// #ifdef CDEBUG
// lprintf(lvl-1, "<-----| [dev_id=%d] CacheBlock::update_state(block_id=%d)\n", Parent->dev_id, id);
// #endif
// 	return ret;
// }

// void CacheBlock::lock(){
// 	// lprintf(0, "|>>>>>> [dev_id=%d] cacheblock \n", Parent->dev_id);
// 	while(__sync_lock_test_and_set(&Lock, 1));
// 	// lprintf(0, "<<<<<<| [dev_id=%d] cacheblock \n", Parent->dev_id);
// 	// Lock++;
// 	// Lock.lock();
// }

// void CacheBlock::unlock(){
// 	__sync_lock_release(&Lock);
// 	// Lock--;
// }

// bool CacheBlock::is_locked(){
// 	if(Lock==0)
// 		return false;
// 	return true;
// }

// /*********************
//  ** Cache Functions **
//  *********************/

// Cache::Cache(int dev_id_in, long long block_num, long long block_size){
// 	// Constructor for caches
// 	short lvl = 2;
// #ifdef CDEBUG
// 	lprintf(lvl-1, "|-----> [dev_id=%d] Cache::Cache(block_num = %lld, block_size = %lld)\n", dev_id_in, block_num, block_size);
// #endif
//   Lock = 0;
// 	id = DevCache_ctr++;
// 	dev_id = dev_id_in;
// 	BlockSize = block_size;
// 	SerialCtr = 0;
// 	BlockNum = block_num;
// 	Size = BlockSize*BlockNum;
// 	Blocks =  (CBlock_p*) malloc (BlockNum * sizeof(CBlock_p));
// 	for (int idx = 0; idx < BlockNum; idx++) Blocks[idx] = new CacheBlock(idx, this, BlockSize); // Or NULL here and initialize when requested? not sure
// 	cont_buf_head = NULL;

// #if defined(FIFO) || defined(MRU) || defined(LRU)
// 	Hash = (Node_LL_p*) malloc(BlockNum * sizeof(Node_LL_p));
// 	InvalidQueue = new LinkedList(this, "InvalidQueue");
// 	for(int idx = 0; idx < BlockNum; idx++){
// 		InvalidQueue->push_back(idx, true);
// 		Hash[idx] = InvalidQueue->end;
// 		Hash[idx]->valid=false;
// 	}
// 	Queue = new LinkedList(this, "Queue");
// #endif
// #ifdef CDEBUG
// 	lprintf(lvl-1, "<-----| [dev_id=%d] Cache::Cache()\n", dev_id_in);
// #endif
// 	return ;
// }

// Cache::~Cache(){
// 	short lvl = 2;
// #ifdef CDEBUG
// lprintf(lvl-1, "|-----> [dev_id=%d] Cache::~Cache()\n", dev_id);
// #endif
// 	lock();
// 	DevCache_ctr--;
// #ifdef ENABLE_CACHE_CONTINUOUS_ALLOC
// 	for (int idx = 0; idx < BlockNum; idx++) if(Blocks[idx]!=NULL) Blocks[idx]->Adrs = NULL;
// 	//if(cont_buf_head)
// 	CHLFree(cont_buf_head, dev_id);
// 	cont_buf_head = NULL;
// #endif
// 	for (int idx = 0; idx < BlockNum; idx++) delete Blocks[idx];
// 	free(Blocks);
// #if defined(FIFO)
// 	free(Hash);
// 	delete InvalidQueue;
// 	delete Queue;
// #endif
// 	unlock();
// #ifdef CDEBUG
// lprintf(lvl-1, "<-----| [dev_id=%d] Cache::~Cache()\n", dev_id);
// #endif
// 	return ;
// }

// void Cache::reset(bool lockfree, bool forceReset){
// 	short lvl = 2;
// #ifdef CDEBUG
// lprintf(lvl-1, "|-----> [dev_id=%d] Cache::reset()\n", dev_id);
// #endif
// 	if(!lockfree) lock();
// 	// Block->reset makes all nodes INVALID and they end up in InvalidQueue.
// 	for (int idx = 0; idx < BlockNum; idx++) Blocks[idx]->reset(lockfree, forceReset);
// #ifdef STEST
// 	timer = 0; // Keeps total time spend in cache operations-code
// #endif
// 	SerialCtr = 0;
// // #if defined(FIFO) || defined(MRU) || defined(LRU)
// // 	if(!lockfree){
// // 		InvalidQueue->lock();
// // 		Queue->lock();
// // 	}
// // 	Node_LL_p node = Queue->start_iterration();
// // 	int i = 0;
// // 	while(i < Queue->length && node->idx >= 0){
// // 		node->valid=false;
// // 		node = Queue->next_in_line();
// // 		i++;
// // 	}
// // 	if(InvalidQueue->length>0){
// // 		if(Queue->length>0){
// // 			InvalidQueue->end->next = Queue->start;
// // 			Queue->start->previous = InvalidQueue->end;
// // 			InvalidQueue->end = Queue->end;
// // 			InvalidQueue->length += Queue->length;
// // 		}
// // 	}
// // 	else{
// // 		InvalidQueue->start = Queue->start;
// // 		InvalidQueue->end = Queue->end;
// // 		InvalidQueue->length = Queue->length;
// // 	}
// // 	Queue->start = NULL;
// // 	Queue->end = NULL;
// // 	Queue->length = 0;
// // 	if(!lockfree){
// // 		InvalidQueue->unlock();
// // 		Queue->unlock();
// // 	}
// // #endif
// 	if(!lockfree) unlock();
// 	// draw_cache(true, true, true);
// #ifdef CDEBUG
// lprintf(lvl-1, "<-----| [dev_id=%d] Cache::reset()\n", dev_id);
// #endif
// 	return ;
// }

// void Cache::draw_cache(bool print_blocks, bool print_queue, bool lockfree){
// 	short lvl = 0;

// 	if(!lockfree)
// 		lock();
// 	lprintf(lvl-1, " Cache:\
// 						\n==========================================\
// 						\n||      Cache Id       | %d\
// 						\n|| - - - - - - - - - - - - - - - - - - - -\
// 						\n||      Device Id      | %d\
// 						\n|| - - - - - - - - - - - - - - - - - - - -\
// 						\n||        Name         | %s\
// 						\n|| - - - - - - - - - - - - - - - - - - - -\
// 						\n||        Size         | %llu\
// 						\n|| - - - - - - - - - - - - - - - - - - - -\
// 						\n||  Number of blocks   | %d\
// 						\n|| - - - - - - - - - - - - - - - - - - - -\
// 						\n||   Size of blocks    | %llu\
// 						\n|| - - - - - - - - - - - - - - - - - - - -", id, dev_id, Name.c_str(), Size, BlockNum, BlockSize);
// #if defined(NAIVE)
// 	lprintf(lvl-1,"\n||  Scheduling policy  | NAIVE");
// #elif defined(FIFO)
// 	lprintf(lvl-1,"\n||  Scheduling policy  | FIFO");
// #elif defined(MRU)
// 	lprintf(lvl-1,"\n||  Scheduling policy  | MRU");
// #elif defined(LRU)
// 	lprintf(lvl-1,"\n||  Scheduling policy  | LRU");
// #endif

// 	lprintf(lvl-1, "\n==========================================");

// 	if(print_blocks){
// 		lprintf(lvl-1, "======================================\
// 							\n|| Start of Blocks in Cache ||\
// 							\n==============================\n");
// 		for(int i=0; i<BlockNum; i++)
// 			if(Blocks[i]!=NULL)
// 				Blocks[i]->draw_block(lockfree);
// 		lprintf(lvl-1, "============================\
// 							\n|| End of Blocks in Cache ||\
// 							\n================================================================================\n");

// 	}
// 	if(print_queue){
// 		#if defined(NAIVE)
// 		lprintf(lvl-1, "There is no Queue.\n");
// 		#elif defined(FIFO) || defined(MRU) || defined(LRU)
// 		lprintf(lvl-1, "|| Start of Queue with Invalid Blocks ||\
// 							\n=======================================\n");
// 		InvalidQueue->draw_queue(lockfree);
// 		lprintf(lvl-1, "============================\
// 							\n|| End of Queue with Invalid Blocks ||\
// 							\n================================================================================\n");
// 		lprintf(lvl-1, "|| Start of Queue with Valid Blocks ||\
// 							\n=======================================\n");
// 		Queue->draw_queue(lockfree);
// 		lprintf(lvl-1, "============================\
// 							\n|| End of Queue with Valid Blocks ||\
// 							\n================================================================================\n");
// 		#endif
// 	}
// 	lprintf(lvl-1, "\n");
// 	if(!lockfree)
// 		unlock();
// }

// #ifdef ENABLE_CACHE_CONTINUOUS_ALLOC
// /// Allocates all cacheblocks in cache, but in a single continuous piece of memory.
// void Cache::allocate(bool lockfree){
// 	short lvl = 2;
// #ifdef CDEBUG
// 	lprintf(lvl-1, "|-----> [dev_id=%d] Cache::allocate-continuous()\n", dev_id);
// #endif
// 	long long total_sz = 0, total_offset = 0;
// 	for(int i=0; i<BlockNum; i++) if(Blocks[i]!=NULL && Blocks[i]->Adrs==NULL) total_sz+= Blocks[i]->Size;
// 	if(!cont_buf_head) cont_buf_head = CHLMalloc(total_sz, dev_id);
// 	for(int i=0; i<BlockNum; i++)
// 		if(Blocks[i]!=NULL){
// 			if(Blocks[i]->Adrs==NULL){
// 				Blocks[i]->Adrs = cont_buf_head + total_offset;
// 				total_offset+=Blocks[i]->Size;
// 			}
// 		}
// 		else error("[dev_id=%d] Cache::allocate-continuous() -> Blocks[%d] was NULL\n", dev_id, i);
// #ifdef CDEBUG
// 	lprintf(lvl-1, "<-----| [dev_id=%d] Cache::allocate-continuous()\n", dev_id);
// #endif
// }

// #else

// void Cache::allocate(bool lockfree){
// 	// Allocates all cacheblocks in cache
// 	short lvl = 2;
// #ifdef CDEBUG
// 	lprintf(lvl-1, "|-----> [dev_id=%d] Cache::allocate()\n", dev_id);
// #endif
// 	for(int i=0; i<BlockNum; i++)
// 		if(Blocks[i]!=NULL) Blocks[i]->allocate(lockfree);
// 		else error("[dev_id=%d] Cache::allocate() -> Blocks[%d] was NULL\n", dev_id, i);
// #ifdef CDEBUG
// 	lprintf(lvl-1, "<-----| [dev_id=%d] Cache::allocate()\n", dev_id);
// #endif
// }
// #endif



// CBlock_p Cache::assign_Cblock(state start_state, bool lockfree){
// 	// Assigns a block from cache to be used for memory.
// 	// State options are:
// 	// - INVALID: Raise error.
// 	// - NATIVE:
// 	// - EXCLUSIVE: Will add one writer.
// 	// - SHARABLE: Will add one reader.
// 	// - AVAILABLE: Will be initialized as AVAILABLE and someone might take it.

// 	short lvl = 2;
// #ifdef CDEBUG
// lprintf(lvl-1, "|-----> [dev_id=%d] Cache::assign_Cblock()\n", dev_id);
// #endif

// 	CBlock_p result = NULL;
// 	if(!lockfree) lock(); // Lock cache
// #if defined(NAIVE)
// 	if (SerialCtr >= BlockNum){
// #endif
// 		int remove_block_idx = -42;
// 	#if defined(NAIVE)
// 		remove_block_idx = CacheSelectBlockToRemove_naive(this, lockfree);
// 	#elif defined(FIFO) || defined(MRU) || defined(LRU)
// 		Node_LL_p remove_block;
// 		remove_block = CacheSelectBlockToRemove_fifo_mru_lru(this, lockfree);
// 		remove_block_idx = remove_block->idx;
// 	#endif
// 		if(remove_block_idx < 0){ // Check again
// 		#if defined(NAIVE)
// 			remove_block_idx = CacheSelectBlockToRemove_naive(this, lockfree);
// 		#elif defined(FIFO) || defined(MRU) || defined(LRU)
// 			remove_block = CacheSelectBlockToRemove_fifo_mru_lru(this, lockfree);
// 			remove_block_idx = remove_block->idx;
// 		#endif
// 			if(remove_block_idx < 0){ // Check for exclusive
// 			#if defined(NAIVE)
// 				remove_block_idx = CacheSelectExclusiveBlockToRemove_naive(this, lockfree);
// 			#elif defined(FIFO) || defined(MRU) || defined(LRU)
// 				remove_block = CacheSelectExclusiveBlockToRemove_fifo_mru_lru(this, lockfree);
// 				remove_block_idx = remove_block->idx;
// 			#endif
// 			}
// 		}
// 		if(remove_block_idx >= 0){
// 			result = Blocks[remove_block_idx];
// 			if(!lockfree){
// 			#if defined(FIFO) || defined(MRU) || defined(LRU)
// 			 	InvalidQueue->lock();
// 				Queue->lock();
// 			#endif
// 				result->lock();
// 			}
// 			// Reset doesn't do anything to the node_ll if block is already INVALID.
// 			result->reset(true,false);

// 	#if defined(FIFO)
// 			Queue->put_last(remove_block, true);
// 	#elif defined(MRU)
// 			Queue->put_first(remove_block, true);
// 	#elif defined(LRU)
// 			Queue->put_last(remove_block, true);
// 	#endif

// 		#ifdef CDUBUG
// 			lprintf(lvl-1,"------ [dev_id=%d] Cache::assign_Cblock(): Block with id=%d reseted.\n", dev_id, remove_block_idx);
// 		#endif
// 			// Set state
// 			if(start_state==INVALID)
// 				error("[dev_id=%d] Cache::assign_Cblock(): New block called to be initialized as invalid\n", dev_id);
// 			else if(start_state==NATIVE)
// 				result->set_state(NATIVE, true);
// 			else if(start_state==EXCLUSIVE){
// 				result->set_state(EXCLUSIVE, true);
// 				result->add_writer(true);
// 			}
// 			else if(start_state==SHARABLE){
// 				result->set_state(SHARABLE, true);
// 				result->add_reader(true);
// 			}
// 			else if(start_state==AVAILABLE){
// 				result->set_state(AVAILABLE, true);
// 			}
// 			else
// 				error("[dev_id=%d] Cache::assign_Cblock(): Uknown state(%s)\n", dev_id, print_state(start_state));
// 		#if defined(FIFO) || defined(MRU) || defined(LRU)
// 			Hash[result->id]->valid = true;
// 		#endif
// 			if(!lockfree){
// 				result->unlock();
// 			#if defined(FIFO) || defined(MRU) || defined(LRU)
// 				Queue->unlock();
// 				InvalidQueue->unlock();
// 			#endif
// 				unlock(); // Unlock cache
// 			}
// 		}
// 		else{
// 	#if defined(FIFO) || defined(MRU) || defined(LRU)
// 			delete remove_block;
// 	#endif

// 			if(!lockfree)  unlock(); // Unlock cache
// 			result = assign_Cblock(start_state, lockfree);
// 		}
// #if defined(NAIVE)
// 	}
// 	else{
// 		result = Blocks[SerialCtr];
// 		if(!lockfree)
// 			result->lock();

// 		result->reset(true,false);
// 		SerialCtr++;
// 		// Set state
// 		if(start_state==INVALID)
// 			error("[dev_id=%d] Cache::assign_Cblock(): New block called to be initialized as invalid\n", dev_id);
// 		else if(start_state==NATIVE)
// 			result->set_state(NATIVE, true);
// 		else if(start_state==EXCLUSIVE){
// 			result->set_state(EXCLUSIVE, true);
// 			result->add_writer(true);
// 		}
// 		else if(start_state==SHARABLE){
// 			result->set_state(SHARABLE, true);
// 			result->add_reader(true);
// 		}
// 		else if(start_state==AVAILABLE){
// 			result->set_state(AVAILABLE, true);
// 		}
// 		else
// 			error("[dev_id=%d] Cache::assign_Cblock(): Uknown state(%s)\n", dev_id, print_state(start_state));

// 		if(!lockfree){
// 			result->unlock();
// 			unlock(); // Unlock cache
// 		}
// 	}
// #endif
// #ifdef CDEBUG
// lprintf(lvl-1, "<-----| [dev_id=%d] Cache::assign_Cblock()\n", dev_id);
// #endif
//   return result;
// }

// void Cache::lock(){
// 	// lprintf(0, "|>>>>>> [dev_id=%d] cache \n", dev_id);
// 	while(__sync_lock_test_and_set(&Lock, 1));
// 	// lprintf(0, "<<<<<<| [dev_id=%d] cache \n", dev_id);
// 	// Lock++;
// 	// Lock.lock();
// }

// void Cache::unlock(){
// 	__sync_lock_release(&Lock);
// 	// Lock--;
// }

// bool Cache::is_locked(){
// 	if(Lock==0)
// 		return false;
// 	return true;
// }

// /*********************
//  ** Other Functions **
//  *********************/

// #if defined(NAIVE)
// int CacheSelectBlockToRemove_naive(Cache_p cache, bool lockfree){
// 	short lvl = 2;

// 	if (cache == NULL)
// 		error("CacheSelectBlockToRemove_naive(): Called on empty buffer\n");

// #ifdef CDEBUG
// 	lprintf(lvl-1, "|-----> [dev_id=%d] CacheSelectBlockToRemove_naive()\n",cache->dev_id);
// #endif

// 	int result_idx = -1;

// 	for (int idx = 0; idx < cache->BlockNum; idx++){ // Iterate through cache serially.
// 		if(!lockfree)
// 			cache->Blocks[idx]->lock();
// 		cache->Blocks[idx]->update_state(true);
// 		state tmp_state = cache->Blocks[idx]->get_state(); // Update all events etc for idx.
// 		if(tmp_state == INVALID || tmp_state == AVAILABLE){ // Indx can be removed if there are no pending events.
// 			result_idx = idx;
// 			cache->Blocks[idx]->set_state(INVALID, true);
// 			if(!lockfree)
// 				cache->Blocks[idx]->unlock();
// 		#ifdef CDEBUG
// 			lprintf(lvl-1, "------- [dev_id=%d] CacheSelectBlockToRemove_naive(): Found available block. Invalidated.\n",cache->dev_id);
// 		#endif
// 			break;
// 		}
// 		if(!lockfree)
// 			cache->Blocks[idx]->unlock();
// 	}
// 	#ifdef CDEBUG
// 		lprintf(lvl-1, "<-----| [dev_id=%d] CacheSelectBlockToRemove_naive()\n",cache->dev_id);
// 	#endif
// 	return result_idx;
// }

// int CacheSelectExclusiveBlockToRemove_naive(Cache_p cache, bool lockfree){
// 	short lvl = 2;

// 	if (cache == NULL)
// 		error("CacheSelectExclusiveBlockToRemove_naive(): Called on empty buffer\n");

// #ifdef CDEBUG
// 	lprintf(lvl-1, "|-----> [dev_id=%d] CacheSelectExclusiveBlockToRemove_naive()\n",cache->dev_id);
// #endif

// 	int result_idx = -1;
// 	CBlock_p native_block;

// 	for (int idx = 0; idx < cache->BlockNum; idx++){ // Iterate through cache serially.
// 		if(!lockfree)
// 			cache->Blocks[idx]->lock();
// 		// cache->Blocks[idx]->update_state(true);
// 		state tmp_state = cache->Blocks[idx]->get_state(); // Update all events etc for idx.
// 		if(tmp_state == EXCLUSIVE){
// 			native_block = cache->Blocks[idx]->WritebackData_p->Native_block;
// 			if(!lockfree)
// 				native_block->lock();
// 			if(cache->Blocks[idx]->PendingReaders==0 && cache->Blocks[idx]->PendingWriters==0){
// 				result_idx = idx;
// 				native_block->add_writer(true);
// 				cache->Blocks[idx]->add_reader(true);
// 				if(!lockfree){
// 					native_block->unlock();
// 					cache->Blocks[idx]->unlock();
// 				}
// 				cache->Blocks[idx]->write_back(true);
// 				if(!lockfree){
// 					cache->Blocks[idx]->lock();
// 					native_block->lock();
// 				}
// 				cache->Blocks[idx]->reset(true, true);
// 				native_block->remove_writer(true);
// 				cache->Blocks[idx]->set_state(INVALID, true);

// 			#ifdef CDEBUG
// 				lprintf(lvl-1, "------- [dev_id=%d] CacheSelectExclusiveBlockToRemove_naive(): Found exclusive block with no pernding operations on it. Invalidated.\n",cache->dev_id);
// 			#endif
// 				if(!lockfree){
// 					native_block->unlock();
// 					cache->Blocks[idx]->unlock();
// 				}
// 				break;
// 			}
// 			if(!lockfree)
// 				native_block->unlock();
// 		}
// 		if(!lockfree)
// 			cache->Blocks[idx]->unlock();
// 	}
// 	#ifdef CDEBUG
// 		lprintf(lvl-1, "<-----| [dev_id=%d] CacheSelectExclusiveBlockToRemove_naive()\n",cache->dev_id);
// 	#endif
// 	return result_idx;
// }

// #elif defined(FIFO) || defined(MRU) || defined(LRU)
// Node_LL_p CacheSelectBlockToRemove_fifo_mru_lru(Cache_p cache, bool lockfree){
// 	short lvl = 2;
// 	if (cache == NULL)
// 		error("CacheSelectBlockToRemove_fifo_mru_lru(): Called on empty buffer\n");

// #ifdef CDEBUG
// lprintf(lvl-1, "|-----> [dev_id=%d] CacheSelectBlockToRemove_fifo_mru_lru()\n", cache->dev_id);
// #endif
// 	Node_LL_p result_node;
// 	if(cache->InvalidQueue->length > 0){
// 		result_node = cache->InvalidQueue->remove(cache->InvalidQueue->start, lockfree);
// 		cache->Blocks[result_node->idx]->set_state(INVALID, lockfree);
// 	}
// 	else{
// 		result_node = new Node_LL();
// 		result_node->idx = -1;
// 		state tmp_state = INVALID;
// 		if(!lockfree){
// 			cache->InvalidQueue->lock();
// 			cache->Queue->lock();
// 		}
// 		Node_LL_p node = cache->Queue->start_iterration();
// 		int i=0;
// 		if(node->idx >= 0){
// 			if(!lockfree)
// 				cache->Blocks[node->idx]->lock();
// 			tmp_state = cache->Blocks[node->idx]->get_state();
// 			while(tmp_state != AVAILABLE){
// 				if(!lockfree)
// 					cache->Blocks[node->idx]->unlock();
// 				node = cache->Queue->next_in_line();
// 				i++;
// 				if(node->idx >= 0 && i < cache->Queue->length){
// 					if(!lockfree)
// 						cache->Blocks[node->idx]->lock();
// 					tmp_state = cache->Blocks[node->idx]->get_state();
// 				}
// 				else
// 					break;
// 			}
// 		}
// 		if(node->idx >=0 && i < cache->Queue->length){
// 			if(tmp_state == AVAILABLE){
// 				delete(result_node);
// 				cache->Blocks[node->idx]->set_state(INVALID, true);
// 				result_node = cache->InvalidQueue->remove(node, true);
// 			#ifdef CDEBUG
// 				lprintf(lvl-1, "------- [dev_id=%d] CacheSelectBlockToRemove_fifo_mru_lru(): Found available block. Invalidated.\n",cache->dev_id);
// 			#endif
// 			}
// 			if(!lockfree)
// 				cache->Blocks[result_node->idx]->unlock();
// 		}
// 		if(!lockfree){
// 			cache->Queue->unlock();
// 			cache->InvalidQueue->unlock();
// 		}
// 	}

// #ifdef CDEBUG
// lprintf(lvl-1, "<-----| [dev_id=%d] CacheSelectBlockToRemove_fifo_mru_lru()\n",cache->dev_id);
// #endif
// 	return result_node;
// }

// Node_LL_p CacheSelectExclusiveBlockToRemove_fifo_mru_lru(Cache_p cache, bool lockfree){
// 	short lvl = 2;
// 	if (cache == NULL)
// 		error("CacheSelectExclusiveBlockToRemove_fifo_mru_lru(): Called on empty buffer\n");

// #ifdef CDEBUG
// lprintf(lvl-1, "|-----> [dev_id=%d] CacheSelectExclusiveBlockToRemove_fifo_mru_lru()\n", cache->dev_id);
// #endif
// 	Node_LL_p result_node;
// 	if(cache->InvalidQueue->length > 0){
// 		result_node = cache->InvalidQueue->remove(cache->InvalidQueue->start, lockfree);
// 		cache->Blocks[result_node->idx]->set_state(INVALID, lockfree);
// 	}
// 	else{

// 		result_node = new Node_LL();
// 		result_node->idx = -1;
// 		state tmp_state = INVALID;
// 		if(!lockfree){
// 			cache->InvalidQueue->lock();
// 			cache->Queue->lock();
// 		}
// 		Node_LL_p node = cache->Queue->start_iterration();
// 		int i=0;
// 		if(node->idx >= 0){
// 			if(!lockfree)
// 				cache->Blocks[node->idx]->lock();
// 			tmp_state = cache->Blocks[node->idx]->get_state(); // Update all events etc for idx.
// 			while(tmp_state != EXCLUSIVE || cache->Blocks[node->idx]->PendingReaders>0 || cache->Blocks[node->idx]->PendingWriters>0){
// 				if(!lockfree)
// 					cache->Blocks[node->idx]->unlock();
// 				node = cache->Queue->next_in_line();
// 				i++;
// 				if(node->idx >= 0 && i < cache->Queue->length){
// 					if(!lockfree)
// 						cache->Blocks[node->idx]->lock();
// 					tmp_state = cache->Blocks[node->idx]->get_state(); // Update all events etc for idx.
// 				}
// 				else
// 					break;
// 			}
// 		}
// 		if(node->idx >=0 && i < cache->Queue->length){
// 			if(tmp_state == EXCLUSIVE){
// 				CBlock_p native_block = cache->Blocks[node->idx]->WritebackData_p->Native_block;
// 				if(!lockfree)
// 					native_block->lock();
// 				if(cache->Blocks[node->idx]->PendingReaders==0 && cache->Blocks[node->idx]->PendingWriters==0 && cache->Hash[node->idx]->valid){
// 					delete result_node;
// 					cache->Hash[node->idx]->valid = false;
// 					// native_block->add_writer(true);
// 					// cache->Blocks[node->idx]->add_reader(true);
// 					if(!lockfree){
// 						native_block->unlock();
// 						cache->Blocks[node->idx]->unlock();
// 						cache->Queue->unlock();
// 						cache->InvalidQueue->unlock();
// 					}
// 					cache->Blocks[node->idx]->write_back();
// 					if(!lockfree){
// 						cache->InvalidQueue->lock();
// 						cache->Queue->lock();
// 						cache->Blocks[node->idx]->lock();
// 						native_block->lock();
// 					}
// 					cache->Blocks[node->idx]->reset(true, true);
// 					// native_block->remove_writer(true);
// 					// cache->Blocks[node->idx]->set_state(INVALID, true);
// 					result_node = cache->InvalidQueue->remove(node, true);
// 				#ifdef CDEBUG
// 					lprintf(lvl-1, "------- [dev_id=%d] CacheSelectExclusiveBlockToRemove_fifo_mru_lru(): Found exclusive block with no pernding operations on it. Invalidated.\n",cache->dev_id);
// 				#endif
// 				}
// 				if(!lockfree)
// 					native_block->unlock();
// 			}
// 			if(!lockfree)
// 				cache->Blocks[result_node->idx]->unlock();
// 		}
// 		if(!lockfree){
// 			cache->Queue->unlock();
// 			cache->InvalidQueue->unlock();
// 		}
// 	}
// #ifdef CDEBUG
// lprintf(lvl-1, "<-----| [dev_id=%d] CacheSelectExclusiveBlockToRemove_fifo_mru_lru()\n",cache->dev_id);
// #endif
// 	return result_node;
// }
// #endif
