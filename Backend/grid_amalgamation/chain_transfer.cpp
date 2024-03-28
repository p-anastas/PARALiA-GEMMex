///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief Some CUDA function calls with added error-checking
///

#include "chl_smart_wrappers.hpp"
#include "chl_grid_amalgamation.hpp"

#ifdef TTEST /// C programmers hate him
int fast_trans_ctr = 0;
long long bytes[100000] = {0};
int transfer_link[100000][2];
int timer_ctr[64][64] = {{0}};
double link_gbytes_s[64][64] = {{0}};
Event_timer_p event_time[100000];

//int hop_log_lock = 0; /// This might slow down things, but it is needed. 

void reseTTEST(){
	for(int k = 0; k < fast_trans_ctr; k++){
		bytes[k] = 0;
		transfer_link[k][0] = transfer_link[k][1] = -42;
		delete event_time[k];
	} 
	fast_trans_ctr = 0;
	for (int d1 = 0; d1 < CHL_MEMLOCS; d1++)
		for (int d2 = 0; d2 < CHL_MEMLOCS; d2++){
			timer_ctr[d1][d2] = 0; 
			link_gbytes_s[d1][d2] = 0; 
		}
}
#endif

void LinkRoute::print(){
	fprintf(stderr,  "LinkRoute(%p)::print (hop_num = %d)\n", this, hop_num);
	fprintf(stderr,  "hop_uid_list = %s\n", printlist<int>(hop_uid_list, hop_num));
	fprintf(stderr,  "hop_ldim_list = %s\n", printlist<int>(hop_ldim_list, hop_num));
	fprintf(stderr,  "hop_buf_list = %s\n", printlist<void*>(hop_buf_list, hop_num));
	fprintf(stderr,  "hop_cqueue_list = %s\n", printlist<void*>((void**)hop_cqueue_list, hop_num-1));
	fprintf(stderr,  "hop_event_list = %s\n", printlist<void*>((void**)hop_event_list, hop_num-1));
}

//#define SPLIT_2D_ROWISE
void FasTCHLMemcpy2DAsync(LinkRoute_p roadMap, long int rows, long int cols, short elemSize){
	if(roadMap->hop_num - roadMap->starting_hop < 2) error("FasTCHLMemcpy2DAsync: Cannot copy with less than 2 locations\n");
#ifdef DDEBUG
	fprintf(stderr, "FasTCHLMemcpy2DAsync(%p, %ld, %ld, %d)\n",
		roadMap, rows, cols, elemSize);
#endif	
#ifdef TTEST
	//while(__sync_lock_test_and_set(&hop_log_lock, 1));
	if (roadMap->hop_num > 5) error("FasTCHLMemcpy2DAsync(dest = %d, src = %d) exeeded 5 hops in TTEST Mode\n",
			roadMap->hop_uid_list[roadMap->starting_hop], roadMap->hop_uid_list[roadMap->hop_num]);
	if (fast_trans_ctr > 100000) error("FasTCHLMemcpy2DAsync(dest = %d, src = %d) exeeded 100000 transfers in TTEST Mode\n",
			roadMap->hop_uid_list[roadMap->starting_hop], roadMap->hop_uid_list[roadMap->hop_num]);
	if(!fast_trans_ctr) reseTTEST();
#endif
	if (roadMap->hop_num - roadMap->starting_hop == 2){
#ifdef TTEST
		transfer_link[fast_trans_ctr][1] = roadMap->hop_uid_list[roadMap->starting_hop];
		transfer_link[fast_trans_ctr][0] = roadMap->hop_uid_list[roadMap->starting_hop + 1];

		//CHLSetTimerAsync(&(transfer_time[fast_trans_ctr][0]));
		event_time[fast_trans_ctr] = new Event_timer(roadMap->hop_cqueue_list[roadMap->starting_hop]->dev_id);
		event_time[fast_trans_ctr]->start_point(roadMap->hop_cqueue_list[roadMap->starting_hop]);
		//roadMap->hop_cqueue_list[roadMap->starting_hop]->add_host_func((void*)&CHLSetTimerAsync, 
		//	(void*) &(transfer_time[fast_trans_ctr][1]));
#endif
		roadMap->hop_cqueue_list[roadMap->starting_hop]->memcpy2DAsync(roadMap->hop_buf_list[roadMap->starting_hop+1], roadMap->hop_ldim_list[roadMap->starting_hop+1],
										roadMap->hop_buf_list[roadMap->starting_hop], roadMap->hop_ldim_list[roadMap->starting_hop],
										rows, cols, elemSize, 
										roadMap->hop_uid_list[roadMap->starting_hop+1], 
										roadMap->hop_uid_list[roadMap->starting_hop], 
										0);
		if(roadMap->hop_event_list[roadMap->starting_hop]) roadMap->hop_event_list[roadMap->starting_hop]->record_to_queue(roadMap->hop_cqueue_list[roadMap->starting_hop]);
#ifdef TTEST
		//roadMap->hop_cqueue_list[roadMap->starting_hop]->add_host_func((void*)&CHLSetTimerAsync, 
		//	(void*) &(transfer_time[fast_trans_ctr][2]));
		event_time[fast_trans_ctr]->stop_point(roadMap->hop_cqueue_list[roadMap->starting_hop]);
		bytes[fast_trans_ctr++] = rows*cols*elemSize;
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
				if(uid_ctr < roadMap->hop_num - 1) step_events[uid_ctr][steps] = new Event();
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
				//if(!steps){
					transfer_link[fast_trans_ctr][1] = roadMap->hop_uid_list[uid_ctr];
					transfer_link[fast_trans_ctr][0] = roadMap->hop_uid_list[uid_ctr+1];
					//CHLSetTimerAsync(&(transfer_time[fast_trans_ctr][0]));
					//roadMap->hop_cqueue_list[uid_ctr]->add_host_func((void*)&CHLSetTimerAsync, 
					//		(void*) &(transfer_time[fast_trans_ctr][1]));
					event_time[fast_trans_ctr] = new Event_timer(roadMap->hop_cqueue_list[uid_ctr]->dev_id);
					event_time[fast_trans_ctr]->start_point(roadMap->hop_cqueue_list[uid_ctr]);
				//}
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
#ifdef TTEST
			//roadMap->hop_cqueue_list[uid_ctr]->add_host_func((void*)&CHLSetTimerAsync, 
			//	(void*) &(transfer_time[fast_trans_ctr][2]));
				//if(!steps){
					event_time[fast_trans_ctr]->stop_point(roadMap->hop_cqueue_list[uid_ctr]);
					bytes[fast_trans_ctr++] = 
#ifdef SPLIT_2D_ROWISE
					local_rows * cols * elemSize;
#else
					rows * local_cols * elemSize;
#endif
				//}
#endif
			}
			if(roadMap->hop_event_list[uid_ctr]) roadMap->hop_event_list[uid_ctr]->record_to_queue(roadMap->hop_cqueue_list[uid_ctr]);
		}
	}
}

#ifdef TTEST
void HopMemcpyPrint(){
	fprintf(stderr,"\n Hop Tranfers Full:\n");
	//FILE* fp = fopen("temp_hop_trans.log", "w+");
	for(int k = 0; k < fast_trans_ctr; k++){
		int src, dest;
		if (best_grid_edge_bws[transfer_link[k][0]][transfer_link[k][1]] == -1.0 && 
				best_grid_edge_replaced[transfer_link[k][0]][transfer_link[k][1]][0] != -1){
				dest = best_grid_edge_replaced[transfer_link[k][0]][transfer_link[k][1]][0];
				src = best_grid_edge_replaced[transfer_link[k][0]][transfer_link[k][1]][1];
			}
		else{
			dest = transfer_link[k][0]; 
			src = transfer_link[k][1];
		}	
		timer_ctr[dest][src]++;
		double time = event_time[k]->sync_get_time()/1000;
		link_gbytes_s[dest][src]+=Gval_per_s(bytes[k], time);
		//fprintf(fp, "%d,%d,%ld,%lf,%lf,%lf\n", transfer_link[k][0], transfer_link[k][1], bytes[k], 
		//	transfer_time[k][0], transfer_time[k][1], transfer_time[k][2]);
	}
		
	fprintf(stderr,"\n Hop Tranfer Map:\n      |");
	for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
		fprintf(stderr, "  %s  |", mem_name(d2));
	fprintf(stderr, "\n      |");
	for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
		fprintf(stderr, "----------");
	fprintf(stderr, "\n");
	for (int d1 = 0; d1 < CHL_MEMLOCS; d1++){
		fprintf(stderr, "%s | ", mem_name(d1));
		for (int d2 = 0; d2 < CHL_MEMLOCS; d2++){
			fprintf(stderr, "%6d  | ", timer_ctr[d1][d2]);
		}
		fprintf(stderr, "\n");
	}

	fprintf(stderr,"\n Hop Tranfer Map Achieved Bandwidths (GB/s):\n      |");
	for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
		fprintf(stderr, "  %s   |", mem_name(d2));
	fprintf(stderr, "\n      |");
	for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
		fprintf(stderr, "-----------");
	fprintf(stderr, "\n");
	for (int d1 = 0; d1 < CHL_MEMLOCS; d1++){
		fprintf(stderr, "%s | ", mem_name(d1));
		for (int d2 = 0; d2 < CHL_MEMLOCS; d2++)
			if (timer_ctr[d1][d2]) fprintf(stderr, "% 6.2lf  | ", link_gbytes_s[d1][d2]/timer_ctr[d1][d2]);
			else fprintf(stderr, "    -    | ");
		fprintf(stderr, "\n");
	}
	//fclose(fp);
	reseTTEST();
}
#endif
