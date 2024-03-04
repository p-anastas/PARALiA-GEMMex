///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief "DataTile" class function implementations.
///

#include "DataTile.hpp"
#include "Resource_manager.hpp"
#include "backend_wrappers.hpp"

int Tile2D_num = 0;

//----------------------------------------------General class-----------------------------------------//

int Tile2D::get_dtype_size() {
    if (dtype == DOUBLE) return sizeof(double);
    else if (dtype == FLOAT) return sizeof(float);
    else error("dtypesize: Unknown type");
    return -1;
}

int Tile2D::size() { return get_dtype_size()*dim1*dim2;}

Tile2D::Tile2D(void *in_addr, int in_dim1, int in_dim2,
               int in_ldim, int inGrid1, int inGrid2, dtype_enum dtype_in, CBlock_p init_loc_block_p){
#ifdef DDEBUG
	fprintf(stderr, "|-----> Tile2D(%d)::Tile2D(in_addr(%d) = %p,%d,%d,%d, %d, %d)\n",
			Tile2D_num, CHLGetPtrLoc(in_addr), in_addr, in_dim1, in_dim2, in_ldim, inGrid1, inGrid2);
#endif
	dtype = dtype_in;
	dim1 = in_dim1;
	dim2 = in_dim2;
	GridId1 = inGrid1;
	GridId2 = inGrid2;
	id = Tile2D_num;
	W_op_dev_id = W_op_queue_ctr = W_op_num = -42;
	W_op_fired = 0; 
	W_op_complete = W_ready = NULL;
	W_op_params = NULL; 
	W_op_dep_num = 0; 
	Tile2D_num++;
	short init_loc = CHLGetPtrLoc(in_addr);
	for (int iloc = 0; iloc < CHL_MEMLOCS; iloc++){
		if (iloc == init_loc){
			W_init_loc = iloc;
			//loc_map[iloc] = 0;
			block_ETA[iloc] = 0; 
			StoreBlock[iloc] = init_loc_block_p;
			StoreBlock[iloc]->Adrs = in_addr;
			//StoreBlock[iloc]->set_owner((void **)&StoreBlock[iloc], false);
			ldim[iloc] = in_ldim;
			StoreBlock[iloc]->Available->record_to_queue(NULL);
		}
		else{
			StoreBlock[iloc] = NULL;
			ldim[iloc] = in_dim1;
			//loc_map[iloc] = -42;
			block_ETA[iloc] = -42; 
		} 
	}
#ifdef DDEBUG
  fprintf(stderr, "<-----|\n");
#endif
}

Tile2D::~Tile2D()
{
  delete W_op_complete; 
  delete W_ready; 
  Tile2D_num--;
  if(!Tile2D_num) warning("Tile2D::~Tile2D destructor incomplete, TBC\n");
}

//----------------------------------------------Tile caching------------------------------------------//

/*void Tile2D::set_loc_idx(int loc_idx, int val){
    loc_map[loc_idx] = val; 
}

void Tile2D::try_set_loc_idx(int loc_idx, int val){
    if (loc_map[loc_idx] == -42) loc_map[loc_idx] = val; 
}*/

void Tile2D::fetch(LinkRoute_p in_route)
{
#ifdef DEBUG
  fprintf(stderr, "|-----> Tile2D(%d:[%d,%d])::fetch(%s)\n", id, GridId1, GridId2,
  	printlist(in_route->hop_uid_list, in_route->hop_num));
#endif
	CBlock_p block_ptr[in_route->hop_num] = {NULL};
	block_ptr[0] = StoreBlock[(in_route->hop_uid_list[0])]; 
	in_route->hop_buf_list[0] = block_ptr[0]->Adrs;
	in_route->hop_ldim_list[0] = ldim[in_route->hop_uid_list[0]];
	for(int inter_hop = 1 ; inter_hop < in_route->hop_num; inter_hop++){
		in_route->hop_ldim_list[inter_hop] = ldim[in_route->hop_uid_list[inter_hop]];
#ifndef PRODUCTION
		if(WRP == RONLY && StoreBlock[(in_route->hop_uid_list[inter_hop])] != NULL) 
			error("Tile2D(%d:[%d,%d])::fetch(%s) -> StoreBlock[%d] is already defined\n", id, GridId1, GridId2, 
				printlist(in_route->hop_uid_list, in_route->hop_num), in_route->hop_uid_list[inter_hop]);
#endif 
		if(WRP == RONLY) block_ptr[inter_hop] = StoreBlock[(in_route->hop_uid_list[inter_hop])] = 
			current_SAB[(in_route->hop_uid_list[inter_hop])]->assign_Cblock(SHARABLE,false);
		else{
			block_ptr[inter_hop] = current_SAB[(in_route->hop_uid_list[inter_hop])]->assign_Cblock(EXCLUSIVE,false);
			if(inter_hop == in_route->hop_num - 1) StoreBlock[(in_route->hop_uid_list[inter_hop])] = block_ptr[inter_hop];
		}
		in_route->hop_buf_list[inter_hop] = block_ptr[inter_hop]->Adrs;
		in_route->hop_event_list[inter_hop-1] = block_ptr[inter_hop]->Available;
		in_route->hop_cqueue_list[inter_hop-1] = 
			recv_queues[(in_route->hop_uid_list[inter_hop])][(in_route->hop_uid_list[inter_hop-1])];
	}
	// Wait until the source of the transfer has the data available (might be obsolete based on implementation)
	//in_route->print();
	in_route->hop_cqueue_list[0]->wait_for_event(block_ptr[0]->Available);

	// Fire the full transfer chain described by in_route
	FasTCHLMemcpy2DAsync(in_route, dim1, dim2, get_dtype_size());

	CQueue_p used_queue = in_route->hop_cqueue_list[in_route->hop_num-2];

	// TODO: Is an updated loc_map needed for PARALiA 3.0?
	//if(WRP == WR || WRP == W_REDUCE || WRP == WR_LAZY) loc_map[(in_route->hop_uid_list[in_route->hop_num-1])] = 42;
	//else for(int inter_hop = 1 ; inter_hop < in_route->hop_num; inter_hop++) loc_map[(in_route->hop_uid_list[inter_hop])] = 42; 

#ifdef DEBUG
	fprintf(stderr, "<-----|\n");
#endif
  return;
}

//--------------------------------------------Tile properties-----------------------------------------//

WR_properties Tile2D::get_WRP(){
    return WRP;
}

void Tile2D::set_WRP(WR_properties inprop){
    this->WRP = inprop;
}

const char* Tile2D::get_WRP_string(){
    switch (WRP)
    {
    case RONLY: 
      return "RONLY";
    case WONLY:  
      return "WONLY";
    case WR:  
      return "WR";
    case WR_LAZY:  
      return "WR_LAZY";   
    case W_REDUCE:  
      return "W_REDUCE";
    default:
        error("Tile2D::get_WRP_string(): Unknown WRP\n");
    }
    return "UNREACHABLE";
}

//--------------------------------------------WTile properties-----------------------------------------//

inline int get_next_queue_ctr(int dev_id){
	exec_queue_ctr[(dev_id)]++;
	if (exec_queue_ctr[(dev_id)] == MAX_BACKEND_L) exec_queue_ctr[(dev_id)] = 0; 
	return exec_queue_ctr[(dev_id)];
}

void Tile2D::run_operation(int W_op_id, LinkRoute_p lazy_route)
{
#ifdef DEBUG
	fprintf(stderr, "|-----> Tile2D(%d:[%d,%d])::run_operation(W_op_dev_id=%d)\n", id, GridId1, GridId2, W_op_dev_id);
#endif
#ifndef PRODUCTION
	if(W_op_dev_id == -42) error("Tile2D(%d:[%d,%d])::run_operation(W_op_dev_id=%d)\n Uninitialized W_op_dev_id", 
    	id, GridId1, GridId2, W_op_dev_id);
#endif 
	CHLSelectDevice(W_op_dev_id);
	CQueue_p assigned_exec_queue = NULL;
	if(W_op_queue_ctr == -42)
		W_op_queue_ctr = get_next_queue_ctr(W_op_dev_id);
	assigned_exec_queue = exec_queue[W_op_dev_id][W_op_queue_ctr];

	if (!W_op_fired && (WRP == WR_LAZY || WRP == W_REDUCE)){
		StoreBlock[W_op_dev_id] = current_SAB[W_op_dev_id]->assign_Cblock(EXCLUSIVE,false);
		assigned_exec_queue->record_event(StoreBlock[W_op_dev_id]->Available);
	}
#ifndef PRODUCTION
	if (StoreBlock[W_op_dev_id] == NULL)
		error("Tile2D(%d:[%d,%d])::run_operation(W_op_dev_id=%d): Write storeblock is NULL\n",
			id, GridId1, GridId2, W_op_dev_id);
#endif
	// Extract the Cblock buffers from the (now defined) Tiles. 
	gemm_backend_in<double>*  ptr_ker_translate = (gemm_backend_in<double>*) W_op_params[W_op_id];
	ptr_ker_translate->A = &((Tile2D_p)ptr_ker_translate->A_tile_v)->StoreBlock[W_op_dev_id]->Adrs;
	ptr_ker_translate->B = &((Tile2D_p)ptr_ker_translate->B_tile_v)->StoreBlock[W_op_dev_id]->Adrs;
	ptr_ker_translate->C = &/*((Tile2D_p)ptr_ker_translate->C_tile_v)->*/StoreBlock[W_op_dev_id]->Adrs;
	// Make execution queue block until all input dependencies are met. 
	assigned_exec_queue->wait_for_event(((Tile2D_p)ptr_ker_translate->A_tile_v)->StoreBlock[W_op_dev_id]->Available);
	assigned_exec_queue->wait_for_event(((Tile2D_p)ptr_ker_translate->B_tile_v)->StoreBlock[W_op_dev_id]->Available);
	assigned_exec_queue->wait_for_event(/*((Tile2D_p)ptr_ker_translate->C_tile_v)->*/StoreBlock[W_op_dev_id]->Available);

	// Fire the W_op_id - ith operation to the queue
	assigned_exec_queue->run_operation(W_op_params[W_op_id], W_op_name, W_op_dev_id);
	W_op_fired++;
	if (W_op_fired == W_op_num){
		if (WRP == WR_LAZY) WR_lazy_combine(lazy_route);
		else if (WRP == W_REDUCE) WReduce_backup_C();
		W_op_complete->record_to_queue(assigned_exec_queue);
//#ifdef ENABLE_SEND_RECV_OVERLAP
//		writeback();
//#endif
	}
#ifndef ASYNC_ENABLE
	CHLSyncCheckErr();
#endif
#ifdef DEBUG
	fprintf(stderr, "<-----|\n");
#endif
}

void Tile2D::writeback(LinkRoute_p out_route){
#ifdef DEBUG
	fprintf(stderr, "|-----> Tile2D(%d:[%d,%d])::writeback(%s)\n", 
    id, GridId1, GridId2, printlist(out_route->hop_uid_list, out_route->hop_num));
#endif
	// If operation location is also the output location, do not perform any writeback. 
	if (W_op_dev_id == W_init_loc) return;
	CBlock_p WB_block = StoreBlock[W_init_loc];
#ifndef PRODUCTION
	if (!(WRP == WR || WRP == WR_LAZY || WRP == W_REDUCE))
    error("Tile2D::writeback -> Tile(%d.[%d,%d]) has WRP = %s\n",
			id, GridId1, GridId2, get_WRP_string());
	if (StoreBlock[W_op_dev_id] == NULL || StoreBlock[W_op_dev_id]->State == INVALID)
		error("Tile2D::writeback -> Tile(%d.[%d,%d]) Storeblock[%d] is NULL\n",
			id, GridId1, GridId2, W_op_dev_id);
	if (WB_block == NULL)
		error("Tile2D::writeback -> Tile(%d.[%d,%d]) WB_block at %d is NULL\n",
      id, GridId1, GridId2, W_init_loc);
	if(W_op_dev_id != out_route->hop_uid_list[0])
		error("Tile2D::writeback error -> W_op_dev_id [%d] != (out_route->hop_uid_list[0]) [%d]", 
		W_op_dev_id, out_route->hop_uid_list[0]);
#endif
    out_route->hop_ldim_list[0] = ldim[W_op_dev_id];
    out_route->hop_buf_list[0] = StoreBlock[W_op_dev_id]->Adrs;
    CBlock_p block_ptr[out_route->hop_num] = {NULL};
    block_ptr[0] = StoreBlock[W_op_dev_id]; 

    for(int inter_hop = 1 ; inter_hop < out_route->hop_num; inter_hop++){
		out_route->hop_ldim_list[inter_hop] = ldim[out_route->hop_uid_list[inter_hop]];
		if(inter_hop < out_route->hop_num - 1){
			block_ptr[inter_hop] = current_SAB[(out_route->hop_uid_list[inter_hop])]->assign_Cblock(EXCLUSIVE,false);
			out_route->hop_event_list[inter_hop-1] = block_ptr[inter_hop]->Available;
		}
		else{
			block_ptr[inter_hop] = WB_block;
			out_route->hop_event_list[inter_hop-1] = W_wb_complete;
		}
		out_route->hop_buf_list[inter_hop] = block_ptr[inter_hop]->Adrs;
    	out_route->hop_cqueue_list[inter_hop-1] = wb_queues[(out_route->hop_uid_list[inter_hop])][(out_route->hop_uid_list[inter_hop-1])];
    }
	out_route->print();
	/// Wait for all operations on the tile to be complete 
    out_route->hop_cqueue_list[0]->wait_for_event(W_op_complete);

#ifdef DEBUG
    fprintf(stderr, "Tile2D::writeback WRP = %s, Road = %s \n", get_WRP_string() ,
    	printlist(out_route->hop_uid_list, out_route->hop_num));
#endif
    FasTCHLMemcpy2DAsync(out_route, dim1, dim2, get_dtype_size());

	if (WRP == W_REDUCE) WReduce_combine();
	else W_ready->record_to_queue(out_route->hop_cqueue_list[out_route->hop_num - 2]);

#ifndef ASYNC_ENABLE
	CHLSyncCheckErr();
#endif
  return; 
}

void Tile2D::WR_lazy_combine(LinkRoute_p lazy_route){
	backup_C = StoreBlock[W_op_dev_id];
    StoreBlock[W_op_dev_id] = current_SAB[W_op_dev_id]->assign_Cblock(EXCLUSIVE,false);
    fetch(lazy_route);
	/// Swap backup_C back to StoreBlock 
	CBlock_p temp_block = StoreBlock[W_op_dev_id]; 
	StoreBlock[W_op_dev_id] = backup_C; 
	backup_C = temp_block; 
	axpy_backend_in<double>* backend_axpy_wrapper = (axpy_backend_in<double>*) malloc(sizeof(struct axpy_backend_in<double>));
    backend_axpy_wrapper->N = dim1*dim2;
    backend_axpy_wrapper->incx = backend_axpy_wrapper->incy = 1;
    backend_axpy_wrapper->alpha = reduce_mult;
    backend_axpy_wrapper->dev_id = W_op_dev_id;
    backend_axpy_wrapper->x = (void**) &(temp_block->Adrs);
    backend_axpy_wrapper->y = (void**) &(StoreBlock[W_op_dev_id]->Adrs);
	/// Wait for WR tile fetch to be complete
    exec_queue[W_op_dev_id][W_op_queue_ctr]->wait_for_event(StoreBlock[W_op_dev_id]->Available);
	/// Perform  C = reduce_mult * C' + C (axpy) at the compute location for this tile (W_op_dev_id)
    exec_queue[W_op_dev_id][W_op_queue_ctr]->run_operation(backend_axpy_wrapper, "Daxpy", W_op_dev_id);
    W_op_complete->record_to_queue(exec_queue[W_op_dev_id][W_op_queue_ctr]);
}

void Tile2D::WReduce_backup_C(){
    //if(W_init_loc >= CHL_WORKERS) W_init_loc = CHL_WORKER_CLOSE_TO_MEMLOC[W_master];
	backup_C = StoreBlock[W_init_loc];
	StoreBlock[W_init_loc] = current_SAB[W_init_loc]->assign_Cblock(EXCLUSIVE,false);
	backup_C_ldim = ldim[W_init_loc];
	ldim[W_init_loc] = dim2; 
	//loc_map[W_op_dev_id] = 42; TODO: P3 - is this needed somewhere?
}

void Tile2D::WReduce_combine(){
    //if(W_init_loc >= CHL_WORKERS) W_init_loc = CHL_WORKER_CLOSE_TO_MEMLOC[W_master];
    if (reduce_queue_ctr[W_init_loc] == REDUCE_WORKERS_PERDEV - 1) reduce_queue_ctr[W_init_loc] = 0; 
    else reduce_queue_ctr[W_init_loc]++;
    CQueue_p WB_exec_queue = reduce_queue[W_init_loc][reduce_queue_ctr[W_init_loc]];
	CBlock_p temp_block = StoreBlock[W_init_loc]; 
	StoreBlock[W_init_loc] = backup_C;
	ldim[W_init_loc] = backup_C_ldim; 

	slaxpby_backend_in<double>* backend_slaxpby_wrapper = (slaxpby_backend_in<double>*) malloc(sizeof(struct slaxpby_backend_in<double>));
    backend_slaxpby_wrapper->N = dim1;
    backend_slaxpby_wrapper->incx = 1;
    backend_slaxpby_wrapper->incy = 1;
    backend_slaxpby_wrapper->alpha = 1.0;
    backend_slaxpby_wrapper->beta = reduce_mult;
    backend_slaxpby_wrapper->dev_id = W_init_loc;
    backend_slaxpby_wrapper->x = (void**) &(temp_block->Adrs);
    backend_slaxpby_wrapper->y = (void**) &(StoreBlock[W_init_loc]->Adrs);
    backend_slaxpby_wrapper->slide_x = dim2;
    backend_slaxpby_wrapper->slide_y = ldim[W_init_loc];

	/// Wait for the writeback of C' to be complete
    WB_exec_queue->wait_for_event(W_wb_complete);

	/// Perform  C = 1.0 * C' + reduce_mult * C (axpby) at the initial data location for this tile (W_init_loc)
	WB_exec_queue->run_operation(backend_slaxpby_wrapper, "Dslaxpby", W_init_loc);

	W_ready->record_to_queue(WB_exec_queue);

}



/*
/// PARALia 2.0 - timed queues and blocks

void Tile2D::ETA_add_task(long double task_duration, int dev_id){
	block_ETA[(dev_id)] += task_duration;
}

void Tile2D::ETA_set(long double new_workload_t, int dev_id){
	block_ETA[(dev_id)] = new_workload_t; 
}

long double Tile2D::ETA_get(int dev_id){
	return block_ETA[(dev_id)];
}

long double Tile2D::ETA_fetch_estimate(int target_id){
  long double result = 0; 
  if(loc_map[(target_id)]){
    int temp_val = loc_map[(target_id)];
    set_loc_idx((target_id), 2);

    LinkRoute_p best_route = new LinkRoute();
    best_route->starting_hop = 0;
    result = best_route->optimize(this, 0);

    set_loc_idx((target_id), temp_val);
  }
  return result; 
}*/

void Tile2D::reset(void* new_adrr, int new_init_chunk, CBlock_p new_init_loc_block_p){
#ifdef DDEBUG
	fprintf(stderr, "|-----> Tile2D()::reset(%p, %d)\n", new_adrr, new_init_chunk);
#endif
	W_op_dev_id = W_op_queue_ctr = W_op_num = -42;
	W_op_fired = 0; 
	W_op_complete = W_ready = NULL;
	W_op_params = NULL; 
	W_op_dep_num = 0; 
	short init_loc = CHLGetPtrLoc(new_adrr);
	short init_loc_idx = (init_loc);
	for (int iloc = 0; iloc < CHL_MEMLOCS; iloc++)
	{
		if (iloc == init_loc_idx){
			StoreBlock[iloc] = new_init_loc_block_p;
			StoreBlock[iloc]->Adrs = new_adrr;
			ldim[iloc] = new_init_chunk;
			StoreBlock[iloc]->Available->record_to_queue(NULL);
		}
		else{
			StoreBlock[iloc] = NULL;
			//loc_map[iloc] = -42;
		} 
	}
#ifdef DDEBUG
	fprintf(stderr, "<-----|\n");
#endif
}