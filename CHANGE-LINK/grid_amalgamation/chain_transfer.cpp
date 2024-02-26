///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief Some CUDA function calls with added error-checking
///

#include "smart_wrappers.hpp"
#include "grid_amalgamation.hpp"

#ifdef TTEST /// C programmers hate him
int fast_trans_ctr = 0;
long long bytes[100000] = {0};
int inter_hop_locs[100000][5];
double inter_hop_timers[100000][4][3];
int timer_ctr[64][64] = {{0}};
double link_gbytes_s[64][64] = {{0}};
int hop_log_lock = 0; /// This might slow down things, but it is needed. 

void reseTTEST(){
	for(int k = 0; k < fast_trans_ctr; k++){
				bytes[k] = 0;
				for(int l = 0; l < 5; l++){
					inter_hop_locs[k][l] = -42;
					if(l < 4) for(int m = 0; m < 3; m++) inter_hop_timers[k][l][m] = 0;
				} 
	}
	fast_trans_ctr = 0;
	for (int d1 = 0; d1 < CHL_MEMLOCS; d1++)
		for (int d2 = 0; d2 < CHL_MEMLOCS; d2++){
			timer_ctr[d1][d2] = 0; 
			link_gbytes_s[d1][d2] = 0; 
		}
}
#endif

//#define SPLIT_2D_ROWISE
void FasTCHLMemcpy2DAsync(LinkRoute_p roadMap, long int rows, long int cols, short elemSize){
	if(roadMap->hop_num - roadMap->starting_hop < 2) error("FasTCHLMemcpy2DAsync: Cannot copy with less than 2 locations\n");
#ifdef TTEST
	while(__sync_lock_test_and_set(&hop_log_lock, 1));
	if (roadMap->hop_num > 5) error("FasTCHLMemcpy2DAsync(dest = %d, src = %d) exeeded 5 hops in TTEST Mode\n",
			roadMap->hop_uid_list[roadMap->starting_hop], roadMap->hop_uid_list[roadMap->hop_num]);
	if (fast_trans_ctr > 100000) error("FasTCHLMemcpy2DAsync(dest = %d, src = %d) exeeded 100000 transfers in TTEST Mode\n",
			roadMap->hop_uid_list[roadMap->starting_hop], roadMap->hop_uid_list[roadMap->hop_num]);
	if(!fast_trans_ctr) reseTTEST();
	bytes[fast_trans_ctr] = rows*cols*elemSize;
#endif
	if (roadMap->hop_num - roadMap->starting_hop == 2){
#ifdef TTEST
		inter_hop_locs[fast_trans_ctr][0] = roadMap->hop_uid_list[roadMap->starting_hop];
		CHLSetTimerAsync(&(inter_hop_timers[fast_trans_ctr][0][0]));
		roadMap->hop_cqueue_list[roadMap->starting_hop]->add_host_func((void*)&CHLSetTimerAsync, 
			(void*) &(inter_hop_timers[fast_trans_ctr][0][1]));
#endif
		roadMap->hop_cqueue_list[roadMap->starting_hop]->memcpy2DAsync(roadMap->hop_buf_list[roadMap->starting_hop+1], roadMap->hop_ldim_list[roadMap->starting_hop+1],
										roadMap->hop_buf_list[roadMap->starting_hop], roadMap->hop_ldim_list[roadMap->starting_hop],
										rows, cols, elemSize, 
										roadMap->hop_uid_list[roadMap->starting_hop+1], 
										roadMap->hop_uid_list[roadMap->starting_hop], 
										0);
		if(roadMap->hop_event_list[roadMap->starting_hop]) roadMap->hop_event_list[roadMap->starting_hop]->record_to_queue(roadMap->hop_cqueue_list[roadMap->starting_hop]);
#ifdef TTEST
		roadMap->hop_cqueue_list[roadMap->starting_hop]->add_host_func((void*)&CHLSetTimerAsync, 
			(void*) &(inter_hop_timers[fast_trans_ctr][0][2]));
#endif
	}
	else{
		if (rows/roadMap->streaming_workers < 100) roadMap->streaming_workers = (rows/100) + 1 ;
		if (cols/roadMap->streaming_workers < 100) roadMap->streaming_workers = (cols/100) + 1 ;
		Event_p step_events[roadMap->hop_num][roadMap->streaming_workers];
		for(int uid_ctr = roadMap->starting_hop; uid_ctr < roadMap->hop_num - 1; uid_ctr++){
#ifdef SPLIT_2D_ROWISE
			long int local_rows = rows/roadMap->streaming_workers;
#else
			long int local_cols = cols/roadMap->streaming_workers;
#endif
			//roadMap->hop_cqueue_list[uid_ctr]->request_parallel_backend();
			for(int steps = 0; steps < roadMap->streaming_workers; steps++){
				if(uid_ctr < roadMap->hop_num - 1) step_events[uid_ctr][steps] = new Event(roadMap->hop_uid_list[uid_ctr+1]);
#ifdef SPLIT_2D_ROWISE
				long buff_offset_dest = steps* elemSize * local_rows,
				buff_offset_src = steps * elemSize * local_rows;
				if(steps == roadMap->streaming_workers -1) local_rows+= rows%roadMap->streaming_workers;
#else
				long buff_offset_dest = steps* elemSize * local_cols * roadMap->hop_ldim_list[uid_ctr + 1],
				buff_offset_src = steps * elemSize * local_cols * roadMap->hop_ldim_list[uid_ctr];
				if(steps == roadMap->streaming_workers -1) local_cols+= cols%roadMap->streaming_workers;
#endif
				if(uid_ctr > 0) roadMap->hop_cqueue_list[uid_ctr]->wait_for_event(step_events[uid_ctr-1][steps]);
#ifdef TTEST
				if(!steps){
					inter_hop_locs[fast_trans_ctr][uid_ctr - roadMap->starting_hop] = roadMap->hop_uid_list[uid_ctr];
					CHLSetTimerAsync(&(inter_hop_timers[fast_trans_ctr][uid_ctr - roadMap->starting_hop][0]));
					roadMap->hop_cqueue_list[uid_ctr]->add_host_func((void*)&CHLSetTimerAsync, 
						(void*) &(inter_hop_timers[fast_trans_ctr][uid_ctr - roadMap->starting_hop][1]));
				}
#endif
				roadMap->hop_cqueue_list[uid_ctr]->memcpy2DAsync(roadMap->hop_buf_list[uid_ctr + 1] + buff_offset_dest, roadMap->hop_ldim_list[uid_ctr + 1],
											roadMap->hop_buf_list[uid_ctr] + buff_offset_src, roadMap->hop_ldim_list[uid_ctr],
#ifdef SPLIT_2D_ROWISE
											local_rows, cols, elemSize,
#else
											rows, local_cols, elemSize,
#endif
											roadMap->hop_uid_list[uid_ctr + 1], roadMap->hop_uid_list[uid_ctr], 0);
				if(uid_ctr < roadMap->hop_num - 1) step_events[uid_ctr][steps]->record_to_queue(roadMap->hop_cqueue_list[uid_ctr]);
			}
			if(roadMap->hop_event_list[uid_ctr]) roadMap->hop_event_list[uid_ctr]->record_to_queue(roadMap->hop_cqueue_list[uid_ctr]);
#ifdef TTEST
		roadMap->hop_cqueue_list[uid_ctr]->add_host_func((void*)&CHLSetTimerAsync, 
			(void*) &(inter_hop_timers[fast_trans_ctr][uid_ctr - roadMap->starting_hop][2]));
#endif
		}
	}
#ifdef TTEST
	inter_hop_locs[fast_trans_ctr][roadMap->hop_num - 1] = roadMap->hop_uid_list[roadMap->hop_num-1];
	fast_trans_ctr++;
	__sync_lock_release(&hop_log_lock);
#endif	
}

#ifdef TTEST
void HopMemcpyPrint(){
	lprintf(0,"\n Hop Tranfers Full:\n");
	FILE* fp = fopen("temp_hop_trans.log", "w+");
	for(int k = 0; k < fast_trans_ctr; k++){
		int src = inter_hop_locs[k][0], dest = inter_hop_locs[k][1], iloc = 1;
		for(int l = 2; l < 5; l++) if(inter_hop_locs[k][l]!= -42){
			dest = inter_hop_locs[k][l];
			iloc = l; 
		}
		int dest_sh, src_sh;
		dest_sh = dest; 
		src_sh = src;
		
		timer_ctr[(dest_sh)][(src_sh)]++;
		double time = (inter_hop_timers[k][iloc-1][2] - inter_hop_timers[k][0][1]), pipe_time = (inter_hop_timers[k][iloc-1][2] - inter_hop_timers[k][0][0]);
		link_gbytes_s[(dest_sh)][(src_sh)]+=Gval_per_s(bytes[k], time);
		//lprintf(0, "Hop Trasfer %d->%d -> road: %s total_t = %lf ms ( %.3lf Gb/s ), pipelined_t = %lf ms ( %.3lf Gb/s )\n", 
		//	inter_hop_locs[k][0], inter_hop_locs[k][iloc], printlist(inter_hop_locs[k], iloc+1),
		//1000*time, Gval_per_s(bytes[k], time), 1000*pipe_time, Gval_per_s(bytes[k], pipe_time));
		fprintf(fp, "%d,%d,%s,%ld,%lf,%lf,%lf\n", inter_hop_locs[k][0], inter_hop_locs[k][iloc], printlist(inter_hop_locs[k], iloc+1), bytes[k], 
			inter_hop_timers[k][0][0], inter_hop_timers[k][0][1], inter_hop_timers[k][iloc-1][2]);
		/*for (int inter_transfers = 0; inter_transfers < iloc ; inter_transfers++){
			double time = (inter_hop_timers[k][inter_transfers][2] - inter_hop_timers[k][inter_transfers][1]), 
			pipe_time = (inter_hop_timers[k][inter_transfers][2] - inter_hop_timers[k][inter_transfers][0]);
			lprintf(1, "link trasfer %d->%d : total_t = %lf ms ( %.3lf Gb/s ), pipelined_t = %lf ms ( %.3lf Gb/s )\n", 
				inter_hop_locs[k][inter_transfers], inter_hop_locs[k][inter_transfers+1], 1000*time, Gval_per_s(bytes[k], time), 1000*pipe_time, Gval_per_s(bytes[k], pipe_time));
		}*/
	}
		
	lprintf(0,"\n Hop Tranfer Map (Full chain):\n      |");
	for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
		lprintf(0, "  %s  |", mem_name(d2));
	lprintf(0, "\n      |");
	for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
		lprintf(0, "----------");
	lprintf(0, "\n");
	for (int d1 = 0; d1 < CHL_MEMLOCS; d1++){
		lprintf(0, "%s | ", mem_name(d1));
		for (int d2 = 0; d2 < CHL_MEMLOCS; d2++){
			lprintf(0, "%6d  | ", timer_ctr[d1][d2]);
		}
		lprintf(0, "\n");
	}

	lprintf(0,"\n Hop Tranfer Map (Full chain) Achieved Bandwidths (GB/s):\n      |");
	for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
		lprintf(0, "  %s   |", mem_name(d2));
	lprintf(0, "\n      |");
	for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
		lprintf(0, "-----------");
	lprintf(0, "\n");
	for (int d1 = 0; d1 < CHL_MEMLOCS; d1++){
		lprintf(0, "%s | ", mem_name(d1));
		for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
			if (timer_ctr[d1][d2]) lprintf(0, "% 6.2lf  | ", link_gbytes_s[d1][d2]/timer_ctr[d1][d2]);
			else lprintf(0, "    -    | ");
		lprintf(0, "\n");
	}

	for (int d1 = 0; d1 < CHL_MEMLOCS; d1++)
		for (int d2 = 0; d2 < CHL_MEMLOCS; d2++){
			timer_ctr[d1][d2] = 0; 
			link_gbytes_s[d1][d2] = 0; 
		}

	for(int k = 0; k < fast_trans_ctr; k++){
		for(int l = 1; l < 5; l++) if(inter_hop_locs[k][l]!= -42){
			int dest_sh, src_sh;
			if (best_grid_edge_bws[(inter_hop_locs[k][l])][(inter_hop_locs[k][l-1])] == -1.0 && 
				best_grid_edge_replaced[(inter_hop_locs[k][l])][(inter_hop_locs[k][l-1])][0] != -1){
				dest_sh = best_grid_edge_replaced[(inter_hop_locs[k][l])][(inter_hop_locs[k][l-1])][0];
				src_sh = best_grid_edge_replaced[(inter_hop_locs[k][l])][(inter_hop_locs[k][l-1])][1];
			}
			else{
				dest_sh = inter_hop_locs[k][l]; 
				src_sh = inter_hop_locs[k][l-1];
			}
			
			timer_ctr[(dest_sh)][(src_sh)]++;
			double time = (inter_hop_timers[k][l-1][2] - inter_hop_timers[k][l-1][1]), 
			pipe_time = (inter_hop_timers[k][l-1][2] - inter_hop_timers[k][l-1][0]);
			link_gbytes_s[(dest_sh)][(src_sh)]+=Gval_per_s(bytes[k], time);
		}
	}

	lprintf(0,"\n Hop Tranfer Map (All hops):\n      |");
	for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
		lprintf(0, "  %s  |", mem_name(d2));
	lprintf(0, "\n      |");
	for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
		lprintf(0, "----------");
	lprintf(0, "\n");
	for (int d1 = 0; d1 < CHL_MEMLOCS; d1++){
		lprintf(0, "%s | ", mem_name(d1));
		for (int d2 = 0; d2 < CHL_MEMLOCS; d2++){
			lprintf(0, "%6d  | ", timer_ctr[d1][d2]);
		}
		lprintf(0, "\n");
	}

	lprintf(0,"\n Hop Tranfer Map (All hops) Achieved Bandwidths (GB/s):\n      |");
	for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
		lprintf(0, "  %s   |", mem_name(d2));
	lprintf(0, "\n      |");
	for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
		lprintf(0, "-----------");
	lprintf(0, "\n");
	for (int d1 = 0; d1 < CHL_MEMLOCS; d1++){
		lprintf(0, "%s | ", mem_name(d1));
		for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
			if (timer_ctr[d1][d2]) lprintf(0, "% 6.2lf   | ", link_gbytes_s[d1][d2]/timer_ctr[d1][d2]);
			else lprintf(0, "   -     | ");
		lprintf(0, "\n");
	}
	fclose(fp);
	reseTTEST();
}
#endif
