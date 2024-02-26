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
	int dim = 1;

	switch (argc) {
	case (3):
		dim = atol(argv[ctr++]);
		loc = atoi(argv[ctr++]);
		break;
	default:
		error("Incorrect input arguments. Usage: ./correct_run N loc use_square_tile:\n"
		"dim : A test buff dimension\n"
		"loc: the initial allocation memory. If = -1 will check smart numa-split malloc.\n");
  	}

	int maxDim = std::min(MAX_DIM_TRANS, (int) CHLGetMaxDimSqAsset2D(2*(CHL_WORKERS), sizeof(double), STEP_TRANS, -1));

	for(short dev_id_idx = 0 ; dev_id_idx < CHL_WORKERS; dev_id_idx++){
		maxDim = std::min(maxDim, (int)CHLGetMaxDimSqAsset2D(2, sizeof(double), STEP_TRANS, (dev_id_idx)));
	}
	int ldhost = std::min((int) (2*MAX_DIM_TRANS), (int) CHLGetMaxDimSqAsset2D(2*2*(CHL_WORKERS), sizeof(double), STEP_TRANS, -1)),
		lddev = maxDim;
	fprintf(stderr,"\nchl_test_mememd: \nSystem = %s\nmaxDim = %d, ldhost = %d, lddev = %d\n", 
		TESTBED, maxDim, ldhost, lddev);
	fprintf(stderr,"------------------------------------------------------------------------------------------------------------------------------------------------------\n");

	//maxDim/=8;
	if (dim <= 0) error("chl_2D_throughput_test: Invalid dim = %d \n", dim);
	if (dim > maxDim) error("chl_2D_throughput_test: Input dim = %d too large for system memories\n", dim);
	int active_unit_num = CHL_WORKERS, active_unit_id_list[CHL_WORKERS], rev_active_unit_num = CHL_WORKERS, rev_active_unit_id_list[CHL_WORKERS];
	double temp_bw, temp_rev_bw, temp_dummy, temp_bid_sim, temp_bid, temp_bid_rev;
	for(short dev_id_idx = 0 ; dev_id_idx < CHL_WORKERS; dev_id_idx++) 
		active_unit_id_list[dev_id_idx] = rev_active_unit_id_list[dev_id_idx] = dev_id_idx;

	double timer = csecond();
	void* loc_buffs[CHL_WORKERS*2], *worker_buffs[CHL_WORKERS*2];
	long long chunk_offset[CHL_WORKERS];
	if(loc == -1){
		long long chunk = PAGE_sz, len_in_pages = ldhost*((long long) ldhost*elemSize)/chunk;
		for(short dev_id_idx = 0 ; dev_id_idx < CHL_WORKERS; dev_id_idx++){
			void* temp_buff1 = CHLMallocHostTouchSmart(ldhost, ldhost, elemSize, 'N');
			void* temp_buff2 = CHLMallocHostTouchSmart(ldhost, ldhost, elemSize, 'N');
			loc_buffs[2*dev_id_idx] = search_sub_addrs_at_memloc(temp_buff1, CHL_WORKER_CLOSE_TO_MEMLOC[dev_id_idx], 
				chunk, ldhost*((long long) ldhost*elemSize), &chunk_offset[dev_id_idx]);
			loc_buffs[2*dev_id_idx + 1] = search_sub_addrs_at_memloc(temp_buff2, CHL_WORKER_CLOSE_TO_MEMLOC[dev_id_idx], 
				chunk, ldhost*((long long) ldhost*elemSize), &chunk_offset[dev_id_idx]);
		}
		timer = csecond() - timer;
		fprintf(stderr, "Allocation loc_buffs size = (%d x %d) x %d complete:\t alloc_timer=%lf ms\n", ldhost, ldhost, elemSize, timer  * 1000);
		ldhost = dim;
		long long pages_per_node[get_hw_numa_num()] = {0};
		int page_node[len_in_pages] = {-1};
		for(short dev_id_idx = 0 ; dev_id_idx < CHL_WORKERS; dev_id_idx++){
			long long remaining_length = len_in_pages - chunk_offset[dev_id_idx];
			get_hw_numa_list(loc_buffs[2*dev_id_idx], chunk, remaining_length, pages_per_node, page_node);
			if (remaining_length < 1024) fprintf(stderr, "Sample loc_buffs[%d][0] in loc %d has %lld pages on numa nodes:\n->pages_per_node = %s\n->page_node = %s\n", 
				dev_id_idx, loc,  remaining_length, printlist(pages_per_node, get_hw_numa_num()), printlist(page_node, remaining_length));
			else fprintf(stderr, "Sample loc_buffs[%d][0] in loc %d has %lld pages on numa nodes:\n->pages_per_node = %s\n->page_node = too_big_to_print\n", 
				dev_id_idx, loc, remaining_length, printlist(pages_per_node, get_hw_numa_num()));
		}
		loc = CHL_MEMLOCS - 1;
	}
	else{
		for(short dev_id_idx = 0 ; dev_id_idx < CHL_WORKERS; dev_id_idx++){
			loc_buffs[2*dev_id_idx] = CHLMalloc(ldhost*((long long) ldhost*elemSize), loc, 1);
			loc_buffs[2*dev_id_idx + 1] = CHLMalloc(ldhost*((long long) ldhost*elemSize), loc, 1);
		}
	}
	timer = csecond();
	for(short dev_id_idx = 0 ; dev_id_idx < CHL_WORKERS; dev_id_idx++){
		worker_buffs[2*dev_id_idx] = CHLMalloc(lddev*lddev*elemSize, (dev_id_idx), 1);
		worker_buffs[2*dev_id_idx + 1] = CHLMalloc(lddev*lddev*elemSize, (dev_id_idx), 1);
	}
	timer = csecond() - timer;
	fprintf(stderr, "Allocation worker_buffs size = (%d x %d) x %d complete:\t alloc_timer=%lf ms\n", lddev, lddev, elemSize, timer  * 1000);
				fprintf(stderr,"------------------------------------------------------------------------------------------------------------------------------------------------------\n");
	//return 1;

	for(short dev_id_idx = 0 ; dev_id_idx < CHL_WORKERS; dev_id_idx++){
		//printf("dev_id = %d, dev_id_idx = %d, dev_id_idy = %d, CHL_WORKERS = %d\n", dev_id, dev_id_idx, dev_id_idy, CHL_WORKERS);
		short queue_id = (loc >= CHL_WORKERS || loc == -1)? dev_id_idx : loc, rev_queue_id = dev_id_idx;
		h2d_queue_list[dev_id_idx] = new CommandQueue(queue_id, COMMUNICATION);
		d2h_queue_list[dev_id_idx] = new CommandQueue(rev_queue_id, COMMUNICATION);
	}

	fprintf(stderr, "Warming up");
	/// Warmup.
	for (int it = 0; it < 3; it++){
		fprintf(stderr, ".");
		for(short dev_id_idx = 0 ; dev_id_idx < CHL_WORKERS; dev_id_idx++){
		if(loc!=dev_id_idx) h2d_queue_list[dev_id_idx]->memcpy2DAsync(worker_buffs[2*dev_id_idx], lddev,
			loc_buffs[2*dev_id_idx], ldhost,
			dim, dim, elemSize,
			(dev_id_idx), loc, 1);
			h2d_queue_list[dev_id_idx]->sync_barrier();
		}
		for(short dev_id_idx = 0 ; dev_id_idx < CHL_WORKERS; dev_id_idx++){
			if(loc!=dev_id_idx) d2h_queue_list[dev_id_idx]->memcpy2DAsync(loc_buffs[2*dev_id_idx + 1], ldhost,
				worker_buffs[2*dev_id_idx + 1], lddev,
				dim, dim, elemSize,
				loc, (dev_id_idx), 1);
			d2h_queue_list[dev_id_idx]->sync_barrier();
		}
	}
	fprintf(stderr, " complete.\n");
	CHLSyncCheckErr();


	perform_microbenchmark_2D(loc_buffs, worker_buffs, loc, (((long long) maxDim*maxDim)*elemSize*active_unit_num), lddev, ldhost, 
		dim, active_unit_id_list, active_unit_num, NULL, -1, NULL, 0, NULL, &temp_bw, &temp_dummy, &temp_dummy);
	perform_microbenchmark_2D(loc_buffs, worker_buffs, loc, (((long long) maxDim*maxDim)*elemSize*active_unit_num), lddev, ldhost, 
		-1, NULL, 0, NULL, dim, rev_active_unit_id_list, rev_active_unit_num, NULL, &temp_dummy, &temp_rev_bw, &temp_dummy);
	perform_microbenchmark_2D(loc_buffs, worker_buffs, loc, (((long long) maxDim*maxDim)*elemSize*active_unit_num), lddev, ldhost, 
		dim, active_unit_id_list, active_unit_num, NULL, dim, rev_active_unit_id_list, rev_active_unit_num, NULL, &temp_bid, &temp_bid_rev, &temp_bid_sim);

	fprintf(stderr, "chl_test_mememd: BW microbench loc --> WORKERS ( %d --> %s ): %.2lf GB/s\n", loc, printlist<int>(active_unit_id_list, active_unit_num), temp_bw); 
	fprintf(stderr, "chl_test_mememd: BW microbench loc <-- WORKERS ( %d <-- %s ): %.2lf GB/s\n", loc, printlist<int>(active_unit_id_list, active_unit_num), temp_rev_bw); 
	fprintf(stderr, "chl_test_mememd: BW microbench loc <-> WORKERS ( %d <-> %s ): %.2lf GB/s [ -> %.2lf GB/s, <- %.2lf GB/s ]\n", 
		loc, printlist<int>(active_unit_id_list, active_unit_num), loc, temp_bid_sim, temp_bid, temp_bid_rev); 

  	return 0;
}
