///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief A transfer microbenchmark aiming to capture how link overlap effects their bandwidth
///

#include <unistd.h>
#include <cassert>

#include <numa.h>

#include "smart_wrappers.hpp"
#include "grid_amalgamation.hpp"
#include "microbenchmarks.hpp"

int main(const int argc, const char *argv[]) {

	int ctr = 1, use_square_tile = 0, loc = CHL_MEMLOCS -1, elem_size = 8;
	long long N = 1;

	switch (argc) {
	case (3):
		N = atol(argv[ctr++]);
		use_square_tile = atoi(argv[ctr++]);
		break;
	default:
		error("Incorrect input arguments. Usage: ./correct_run N use_square_tile:\n"
		"N : A test buff dimension\n"
		"use_square_tile: if 0, allocate 1D vector of size N, else a NxN 2D matrix\n");
  	}

	void* test_buff[CHL_MEMLOCS];
	void* test_buff_2[CHL_MEMLOCS];
	for(loc = 0; loc < CHL_MEMLOCS; loc++){
		if(use_square_tile) test_buff[loc] = CHLMalloc(N*N*elem_size, loc, 1);
		else test_buff[loc] = CHLMalloc(N*elem_size, loc, 1);
		int test_loc = CHLGetPtrLoc(test_buff[loc]);
		if(test_loc!= loc) error("chl_test_mememd: Buffer that should be at %d returned loc =%d\n", loc, test_loc);
		else fprintf(stderr, "chl_test_mememd: Allocation and GetPtrLoc at correct location %d\n", loc);
		void* test_buff_offset = (char*)test_buff[loc] + (N*elem_size)/2;
		test_loc = CHLGetPtrLoc(test_buff_offset);
		if(test_loc!= loc) error("chl_test_mememd: Buffer that should be at %d returned loc =%d\n", loc, test_loc);
		else fprintf(stderr, "chl_test_mememd: GetPtrLoc + offset at correct location %d\n", loc);
		int loc2 = (loc == CHL_MEMLOCS -1) ? loc - 1 : loc + 1;
		test_buff_2[loc] = CHLMalloc(N*elem_size, loc2, 1);
		int test_loc2 = CHLGetPtrLoc(test_buff_2[loc]);
		if(test_loc2!= loc2) error("chl_test_mememd: Buffer that should be at %d returned loc =%d\n", loc2, test_loc2);
		else fprintf(stderr, "chl_test_mememd: 2nd allocation and GetPtrLoc at correct location %d\n", loc2);
		test_buff_offset = (char*)test_buff_2[loc] + N*elem_size/2;
		test_loc2 = CHLGetPtrLoc(test_buff_offset);
		if(test_loc2!= loc2) error("chl_test_mememd: Buffer that should be at %d returned loc =%d\n", loc2, test_loc2);
		else fprintf(stderr, "chl_test_mememd: GetPtrLoc + offset at correct location %d\n", loc2);

		test_buff_offset = (char*)test_buff[loc] + 1;
		test_loc2 = CHLGetPtrLoc(test_buff_offset);
		if(test_loc2!= loc) error("chl_test_mememd: Buffer that should be at %d returned loc =%d\n", loc, test_loc2);
		else fprintf(stderr, "chl_test_mememd: GetPtrLoc in first buffer at %d still correct\n", loc);

		//for (long long itter = 0; itter < allocated_data_num; itter++)
		//	allocated_data[itter]->print_slices(); 		
	}
  	return 0;
}
