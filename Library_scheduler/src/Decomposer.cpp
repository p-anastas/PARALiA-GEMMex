///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief The "Asset" related function implementations.
///

#include "Decomposer.hpp"
#include "chl_smart_wrappers.hpp"

Decom2D::Decom2D( void* in_adr, int in_dim1, int in_dim2, int in_ldim, char in_transpose, dtype_enum dtype_in){
#ifdef DEBUG
  	fprintf(stderr, "|-----> Decom2D::Decom2D(%p, %d,%d, %d, trans, dtype)\n", in_adr, in_dim1, in_dim2, in_ldim);
#endif
  ldim = in_ldim;
  dim1 = in_dim1;
  dim2 = in_dim2;
  adrs = in_adr;
  loc = CHLGetPtrLoc(in_adr);
  transpose = in_transpose;
  dtype = dtype_in;
}

Decom1D::Decom1D(void* in_adr, int in_dim, int in_inc, dtype_enum dtype_in){
  dim1 = in_dim;
  dim2 = 1; 
  GridSz2 = 1; 
  adrs = in_adr;
  loc = CHLGetPtrLoc(in_adr);
  inc = in_inc;
  dtype = dtype_in;
}

void Decomposer::Reset(void* new_adrs, int T1, int T2, long new_chunk_size, Buffer_p* init_loc_cache_p){
    error("Must never be called for parent Decomposer class\n");
    return;
}

void Decom2D::Reset(void* new_adrs, int new_T1, int new_T2, long new_chunk_size, Buffer_p* init_loc_cache_p){
  #ifdef DEBUG
  	fprintf(stderr, "|-----> Decom1D::Reset(Buffer_p=%p, loc = %d)\n", init_loc_cache_p, loc);
  #endif
  adrs = new_adrs;
  ldim = new_chunk_size; 
  int current_ctr;
  void* tile_addr = NULL;
  int loc_idx = loc;
  for (int itt1 = 0; itt1 < GridSz1; itt1++){
		for (int itt2 = 0 ; itt2 < GridSz2; itt2++){
      current_ctr = itt1*GridSz2 + itt2;
      if (transpose == 'N'){
         if (dtype == DOUBLE) tile_addr = ((double*)adrs) + (itt1*new_T1 + itt2*new_T2*ldim);
         else if(dtype == FLOAT) tile_addr = ((float*)adrs) + (itt1*new_T1 + itt2*new_T2*ldim);
         else error("Decom2D::Reset: dtype not implemented");
         Tile_map[current_ctr]->reset(tile_addr, ldim, init_loc_cache_p[loc_idx]->assign_Cblock(NATIVE, true));
       }
      else if (transpose == 'T'){
        if (dtype == DOUBLE) tile_addr = ((double*)adrs) + (itt1*new_T1*ldim + itt2*new_T2);
         else if(dtype == FLOAT)  tile_addr = ((float*)adrs) + (itt1*new_T1*ldim + itt2*new_T2);
        else error("Decom2D::Reset: dtype not implemented");
        Tile_map[current_ctr]->reset(tile_addr, ldim, init_loc_cache_p[loc_idx]->assign_Cblock(NATIVE, true));
      }
      else error("Decom2D::Reset: Unknown transpose type\n");

     }
   }
   #ifdef DEBUG
   	fprintf(stderr, "<-----|\n");
   #endif
}

void Decom1D::Reset(void* new_adrs, int new_T1, int new_T2, long new_chunk_size, Buffer_p* init_loc_cache_p){
  #ifdef DEBUG
  	fprintf(stderr, "|-----> Decom1D::Reset(Buffer_p=%p, loc = %d)\n", init_loc_cache_p, loc);
  #endif
  adrs = new_adrs;
  inc = new_chunk_size;
  int current_ctr;
  void* tile_addr = NULL;
  int loc_idx = loc;
  for (int itt = 0; itt < GridSz1; itt++){
    current_ctr = itt;
    if (dtype == DOUBLE) tile_addr = ((double*)adrs) + itt*new_T1*inc;
    else if(dtype == FLOAT) tile_addr = ((float*)adrs) + itt*new_T1*inc;
    else error("Decom1D::Reset: dtype not implemented");
        Tile_map[current_ctr]->reset(tile_addr, inc, init_loc_cache_p[loc_idx]->assign_Cblock(NATIVE, true));
  }
   #ifdef DEBUG
   	fprintf(stderr, "<-----|\n");
   #endif
}

void Decomposer::InitTileMap(int T1, int T2, Buffer_p* init_loc_cache_p){
    error("Must never be called for parent Decomposer class\n");
    return;
}

void Decom2D::InitTileMap(int T1, int T2, Buffer_p* init_loc_cache_p){
  #ifdef DEBUG
  	fprintf(stderr, "|-----> Decom2D::InitTileMap(%d,%d)\n", T1, T2);
  #endif

  GridSz1 = dim1/T1;
	GridSz2 = dim2/T2;
  int T1Last = dim1%T1, T2Last = dim2%T2;
  // TODO: Padding instead of resize so all data fit in buffer without complex mechanism.
  // Can degrade performance for small div sizes.
  if (T1Last > 0) GridSz1++;
  else T1Last=T1;
  if (T2Last > 0) GridSz2++;
  else T2Last=T2;
	//if (T1Last > T1/4) GridSz1++;
	//else T1Last+=T1;
  //if (T2Last > T2/4) GridSz2++;
  //else T2Last+=T2;

  Tile_map = (DataTile_p*) malloc(sizeof(DataTile_p)*GridSz1*GridSz2);

  int current_ctr, T1tmp, T2tmp;
  void* tile_addr = NULL;
  int loc_idx = CHLGetPtrLoc(adrs);
  for (int itt1 = 0; itt1 < GridSz1; itt1++){
    if ( itt1 == GridSz1 - 1) T1tmp = T1Last;
    else  T1tmp = T1;
		for (int itt2 = 0 ; itt2 < GridSz2; itt2++){
      if ( itt2 == GridSz2 - 1) T2tmp = T2Last;
      else  T2tmp = T2;
      current_ctr = itt1*GridSz2 + itt2;
      /// For column major format assumed with T1tmp = rows and T2tmp = cols
      if (transpose == 'N'){
         if (dtype == DOUBLE) tile_addr = ((double*)adrs) + (itt1*T1 + itt2*T2*ldim);
         else if(dtype == FLOAT) tile_addr = ((float*)adrs) + (itt1*T1 + itt2*T2*ldim);
         else error("Decom2D::InitTileMap: dtype not implemented");
         Tile_map[current_ctr] = new Tile2D(tile_addr, T1tmp, T2tmp, ldim, itt1, itt2, dtype, init_loc_cache_p[loc_idx]->assign_Cblock(NATIVE, true));
       }
      else if (transpose == 'T'){
        if (dtype == DOUBLE) tile_addr = ((double*)adrs) + (itt1*T1*ldim + itt2*T2);
         else if(dtype == FLOAT)  tile_addr = ((float*)adrs) + (itt1*T1*ldim + itt2*T2);
        else error("Decom2D::InitTileMap: dtype not implemented");
        Tile_map[current_ctr] = new Tile2D(tile_addr, T2tmp, T1tmp, ldim, itt2, itt1, dtype, init_loc_cache_p[loc_idx]->assign_Cblock(NATIVE, true));
      }
      else error("Decom2D::InitTileMap: Unknown transpose type\n");

     }
   }
   #ifdef DEBUG
   	fprintf(stderr, "<-----|\n");
   #endif
}

void Decom1D::InitTileMap(int T, int dummy, Buffer_p* init_loc_cache_p){
  #ifdef DEBUG
  	fprintf(stderr, "|-----> Decom1D::InitTileMap(%d)\n", T);
  #endif

  GridSz1 = dim1/T;
  int TLast = dim1%T;
  // TODO: Padding instead of resize so all data fit in buffer without complex mechanism.
  // Can degrade performance for small div sizes.
  if (TLast > 0) GridSz1++;
  else TLast=T;

  Tile_map = (DataTile_p*) malloc(sizeof(DataTile_p)*GridSz1);

  int current_ctr, Ttmp;
  void* tile_addr = NULL;
  for (int itt = 0; itt < GridSz1; itt++){
    if ( itt == GridSz1 - 1) Ttmp = TLast;
    else  Ttmp = T;
    current_ctr = itt;
    if (dtype == DOUBLE) tile_addr = ((double*)adrs) + itt*T*inc;
    else if(dtype == FLOAT) tile_addr = ((float*)adrs) + itt*T*inc;
    else error("Decom1D::InitTileMap: dtype not implemented");
    Tile_map[current_ctr] = new Tile1D(tile_addr, Ttmp, inc, itt, dtype, init_loc_cache_p[CHLGetPtrLoc(adrs)]->assign_Cblock(NATIVE, true));
  }
  #ifdef DEBUG
  	fprintf(stderr, "<-----|\n");
  #endif
}

DataTile_p Decomposer::getTile(int iloc1, int iloc2){
  if(iloc1 >= GridSz1) error("Decomposer::getTile : iloc1 >= GridSz1 (%d vs %d)\n", iloc1, GridSz1);
  else if(iloc2 >= GridSz2) error("Decomposer::getTile : iloc2 >= GridSz2 (%d vs %d)\n", iloc2, GridSz2);
  return Tile_map[iloc1*GridSz2 + iloc2];
}

void Decomposer::DestroyTileMap(){
  int current_ctr;
  for (int itt1 = 0; itt1 < GridSz1; itt1++)
    for (int itt2 = 0 ; itt2 < GridSz2; itt2++){
      current_ctr = itt1*GridSz2 + itt2;
      delete Tile_map[current_ctr];
    }
  free(Tile_map);
}

void Decomposer::WBTileMap(){
  if (Tile_map[0]->WRP == W_REDUCE) error("Decomposer::WBTileMap not implemented for ALGO_WREDUCE\n");
  for (int itt1 = 0; itt1 < GridSz1; itt1++)
    for (int itt2 = 0 ; itt2 < GridSz2; itt2++)
      // TODO: this does not reuse pathing (for !SK_FIRE_WHEN_READY)
      Tile_map[itt1*GridSz2 + itt2]->writeback(NULL, NULL);
}

void Decomposer::SyncTileMap(){
  for (int itt1 = 0; itt1 < GridSz1; itt1++)
    for (int itt2 = 0 ; itt2 < GridSz2; itt2++){
      Tile_map[itt1*GridSz2 + itt2]->W_complete->sync_barrier();
      if(!strcmp(OUTPUT_ALGO_MODE,"ALGO_WREDUCE")) Tile_map[itt1*GridSz2 + itt2]->
        W_reduce->sync_barrier();
      //if(!strcmp(OUTPUT_ALGO_MODE,"ALGO_WREDUCE")) for (int devi = 0; devi < CHL_MEMLOCS; devi++)
      //for (int ctri = 0; ctri < REDUCE_WORKERS_PERDEV; ctri++) 
      //if(reduce_queue[devi] && reduce_queue[devi][ctri]) reduce_queue[devi][ctri]->sync_barrier();
      //Tile_map[itt1*GridSz2 + itt2]->StoreBlock[Tile_map[itt1*GridSz2 + itt2]->
      //get_initial_location()]->Available->sync_barrier();
    }
}


void Decomposer::DrawTileMap(){
  fprintf(stderr, " Decomposer representation: \
                 \n ______________________ \
                 \n| id[GridId1, GridId2] |\
                 \n| - - - - - - - - - - -|\
                 \n|    (dim1 X dim2)     |\
                 \n| - - - - - - - - - - -|\
                 \n| - - WR_properties - -|\
                 \n| - - - - - - - - - - -|\
                 \n| - - - loc_list  - - -|\
                 \n|______________________|\n\n");

  for (int itt1 = 0; itt1 < GridSz1; itt1++){
    for (int itt2 = 0 ; itt2 < GridSz2; itt2++)
      fprintf(stderr, " ______________________ ");
    fprintf(stderr, "\n");
    for (int itt2 = 0 ; itt2 < GridSz2; itt2++)
      fprintf(stderr, "|%4d[%6d,%6d]   |",
      Tile_map[itt1*GridSz2 + itt2]->id,
      Tile_map[itt1*GridSz2 + itt2]->GridId1,
      Tile_map[itt1*GridSz2 + itt2]->GridId2);
    fprintf(stderr, "\n");
    for (int itt2 = 0 ; itt2 < GridSz2; itt2++)
      fprintf(stderr, "| - - - - - - - - - - -|");
    fprintf(stderr, "\n");
    for (int itt2 = 0 ; itt2 < GridSz2; itt2++)
      fprintf(stderr, "|  (%6d X %6d)   |",
      Tile_map[itt1*GridSz2 + itt2]->dim1,
      Tile_map[itt1*GridSz2 + itt2]->dim2);
    fprintf(stderr, "\n");
    for (int itt2 = 0 ; itt2 < GridSz2; itt2++)
      fprintf(stderr, "| - - - - - - - - - - -|");
    fprintf(stderr, "\n");
    for (int itt2 = 0 ; itt2 < GridSz2; itt2++)
      fprintf(stderr, "%s", Tile_map[itt1*GridSz2 + itt2]->get_WRP_string());
    fprintf(stderr, "\n");
    for(int loctr = 0; loctr < CHL_MEMLOCS; loctr++){
      for (int itt2 = 0 ; itt2 < GridSz2; itt2++)
        fprintf(stderr, "| - - - - - - - - - - -|");
      fprintf(stderr, "\n");
      for (int itt2 = 0 ; itt2 < GridSz2; itt2++)
        fprintf(stderr, "%s", printlist<int>(Tile_map[itt1*GridSz2 + itt2]->loc_map, CHL_MEMLOCS));
      fprintf(stderr, "\n");
    }
    for (int itt2 = 0 ; itt2 < GridSz2; itt2++)
     fprintf(stderr, "|______________________|");
    fprintf(stderr, "\n\n");
  }
}


long Decomposer::get_chunk_size(){
    error("Must never be called for parent Decomposer class\n");
    return -42;
}


long Decom2D::get_chunk_size(){
    return ldim;
}

long Decom1D::get_chunk_size(){
    return inc;
}

long Decomposer::get_mem_size(){
    error("Must never be called for parent Decomposer class\n");
    return -42;
}

long Decom2D::get_mem_size(){
    return dim2*ldim*dtypesize();
}


long Decom1D::get_mem_size(){
    return dim1*inc*dtypesize();
}

void Decomposer::set_chunk_size(long value){
    error("Must never be called for parent Decomposer class\n");
    return;
}

void Decom2D::set_chunk_size(long value){
    ldim = value;
}

void Decom1D::set_chunk_size(long value){
    inc = value;
}