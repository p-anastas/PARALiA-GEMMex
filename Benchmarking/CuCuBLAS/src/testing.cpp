///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief Some testing helper functions for file output.
///

#include <cstdlib>
#include <cmath>
#include "PARALiA.hpp"

char* CoCoImplementationPrint(){
	char* string_out = (char*) malloc (2048*sizeof(char));
	char* string_helper = (char*) malloc (1024*sizeof(char));
#ifndef ASYNC_ENABLE
	sprintf(string_helper, "_NO-ASYNC");
	strcat(string_out, string_helper);
#endif
#ifndef BUFFER_REUSE_ENABLE
	sprintf(string_helper, "_NO-BUF-RE");
	strcat(string_out, string_helper);
#endif
#ifndef QUEUE_REUSE_ENABLE
	sprintf(string_helper, "_NO-UN-RE");
	strcat(string_out, string_helper);
#endif
	sprintf(string_helper, "_COMP_STREAMS_PERDEV-%d", MAX_BACKEND_L);
	strcat(string_out, string_helper);
	sprintf(string_helper, "_COMM_STREAM_BUFFERING-%d", STREAMING_BUFFER_OVERLAP);
	strcat(string_out, string_helper);
	if (!strcmp(OUTPUT_ALGO_MODE,"ALGO_WR"))
	sprintf(string_helper, "_ALGO-BASIC");
	else if (!strcmp(OUTPUT_ALGO_MODE,"ALGO_WR_LAZY"))
	sprintf(string_helper, "_ALGO-WR-LAZY");
	else if (!strcmp(OUTPUT_ALGO_MODE,"ALGO_WREDUCE"))
	sprintf(string_helper, "_ALGO-WREDUCE-%d", REDUCE_WORKERS_PERDEV);
	else if (!strcmp(OUTPUT_ALGO_MODE,"ALGO_AUTO"))
	sprintf(string_helper, "_ALGO_AUTO");
	strcat(string_out, string_helper);
#ifndef ENABLE_SEND_RECV_OVERLAP
	sprintf(string_helper, "_NO-SND-RCV-OVER");
	strcat(string_out, string_helper);
#endif
#ifdef SUBKERNELS_FIRE_LAZY
	sprintf(string_helper, "_SK-FIRE-LAZY");
	strcat(string_out, string_helper);
#endif
//#ifndef ENABLE_CPU_WORKLOAD
//	sprintf(string_helper, "_NO-CPU");
//	strcat(string_out, string_helper);
//#endif
/// Tile selection algorithm MD should go here, if any 
#ifdef APPLY_TILE_SL_TO_WORKLOAD_SPLIT
	sprintf(string_helper, "_MODEL-TILE-SL");
	strcat(string_out, string_helper);
#endif

	if (!strcmp(DISTRIBUTION, "2D-BLOCK-CYCLIC")) sprintf(string_helper, "_2D-BLOCK-CYCLIC-%s", ORDER_2DBC);
	else error("CoCoImplementationPrint(): Unknown Subkernel Distribution %s\n", DISTRIBUTION);
	strcat(string_out, string_helper);

	if (!strcmp(FETCH_ROUTING, "P2P_FETCH_FROM_INIT")) sprintf(string_helper, "_FETCH-P2P-INIT");
	else if (!strcmp(FETCH_ROUTING, "P2P_FETCH_FROM_GPU_SERIAL")) sprintf(string_helper, "_FETCH-P2P-SERIAL");
	else if (!strcmp(FETCH_ROUTING, "P2P_FETCH_FROM_GPU_DISTANCE")) sprintf(string_helper, "_FETCH-P2P-DISTANCE");
	else if (!strcmp(FETCH_ROUTING, "CHAIN_FETCH_SERIAL")) sprintf(string_helper, "_FETCH-CHAIN-SERIAL");
	else if (!strcmp(FETCH_ROUTING, "CHAIN_FETCH_RANDOM")) sprintf(string_helper, "_FETCH-CHAIN-RANDOM");
	else if (!strcmp(FETCH_ROUTING, "CHAIN_FETCH_TIME")) sprintf(string_helper, "_FETCH-CHAIN-TIME");
	else if (!strcmp(FETCH_ROUTING, "CHAIN_FETCH_QUEUE_WORKLOAD")) sprintf(string_helper, "_FETCH-CHAIN-QETA");
	else error("CoCoImplementationPrint(): Unknown Fetch routing %s\n", FETCH_ROUTING);
	strcat(string_out, string_helper);

	if (!strcmp(WB_ROUTING, "P2P_TO_INIT")) sprintf(string_helper, "_WB-P2P-INIT");
	else error("CoCoImplementationPrint(): Unknown WB routing %s\n", WB_ROUTING);
	strcat(string_out, string_helper);

	if (!strcmp(TASK_ORDER, "SERIAL")) sprintf(string_helper, "_TASK-ORDER-SERIAL");
	else if (!strcmp(TASK_ORDER, "FETCH_MINFETCH")) sprintf(string_helper, "_TASK-ORDER-MINFETCH");
	else if (!strcmp(TASK_ORDER, "FETCH_MINFETCH_THEN_MINPENDING")) sprintf(string_helper, "_TASK-ORDER-MINFETCH-THEN-MINPENDING");
	else if (!strcmp(TASK_ORDER, "FETCH_ETA")) sprintf(string_helper, "_TASK-ORDER-QETA");
	else if (!strcmp(TASK_ORDER, "FETCH_ETA_PLUS_MINPENDING")) sprintf(string_helper, "_TASK-ORDER-QETA-PLUS-MINPENDING");
	else error("CoCoImplementationPrint(): Unknown Task Order Algorithm %s\n", TASK_ORDER);
	strcat(string_out, string_helper);

#ifdef ENABLE_POWA
	if(strcmp(PREDICT_OPTIMIZE_TARGET,"PERF-PER-J")) sprintf(string_helper, "_PW-PRED-%s", PREDICT_OPTIMIZE_TARGET);
	else sprintf(string_helper, "_PW-PRED-%s-%.2lf", PREDICT_OPTIMIZE_TARGET, PERPER_LIMIT);
	strcat(string_out, string_helper);
#endif
	sprintf(string_helper, "_REPTILE-%d", REP_TILE);
	strcat(string_out, string_helper);
//#ifdef DDEBUG
	printf("%s\n", string_out);
//#endif
	free(string_helper);
	return string_out;

}

void ParseInputLvl3(const int argc, const char *argv[], ATC_p* predef_control_values,
		char* TransA, char* TransB, double* alpha, double* beta, long int* D1, long int* D2, long int* D3,
		int* loc1, int* loc2, int* loc3, int* outloc){
	if(argc != 16) error("Incorrect input arguments. Usage: ./correct_run\
	\n\tactive_unit_num(auto if <0)\n\tdev_ids(form example: 0101 for devices 0,2 - ignore if active_unit_num < 0)\n\tT(auto if <=0)\
	\n\tcache_max_size(auto if <0)\n\tTransA TransB alpha beta D1 D2 D3 loc1 loc2 loc3 outloc \n");

	int active_unit_num = atoi(argv[1]);
	long int dev_ids_token = atoi(argv[2]);
	int T = atoi(argv[3]);
	long long cache_limit = atof(argv[4]);

	ATC_p temp_p;

	//Tunning Parameters
	if (active_unit_num < 0 && cache_limit <= 0 && T <=0 ) temp_p = *predef_control_values = NULL;
	else{
		fprintf(stderr, "Using predefined control parameters from input\n");
		*predef_control_values = new ATC();
		temp_p = *predef_control_values;
		temp_p->cache_limit = cache_limit;
		temp_p->T = T;
		if (active_unit_num < 0){
			temp_p->active_unit_num = -1;
		}
		else{
			temp_p->active_unit_num = active_unit_num;
			int ctr = 0, itter = 0;
			do {
				if (dev_ids_token%10 == 1) temp_p->active_unit_id_list[ctr++] = (itter);
				itter++;
				dev_ids_token/=10;
			} while ( dev_ids_token > 0);
			if (ctr != active_unit_num) error("ParseInputLvl3: Read different device Ids in total (%d) than active_unit_num implied (%d)\n", ctr, active_unit_num);
		}
	}

	//Problem Parameters
	*TransA = argv[5][0];
	*TransB = argv[6][0];
	*alpha = atof(argv[7]);
	*beta = atof(argv[8]);
	*D1 = atoi(argv[9]);
	*D2 = atoi(argv[10]);
	*D3 = atoi(argv[11]);
	*loc1 = atoi(argv[12]);
	*loc2 = atoi(argv[13]);
	*loc3 = atoi(argv[14]);
	*outloc = atoi(argv[15]);

	fprintf(stderr, "ParseInputLvl3: ");
	if (*predef_control_values) (*predef_control_values)->print();
	fprintf(stderr, "Routine values:\n\tTransA: %c, TransB: %c\n\talpha: %lf, beta: %lf\n\tD1: %zu, D2: %zu, D3: %zu\n\tloc1: %d, loc2: %d, loc3: %d, outloc: %d\n",
	*TransA, *TransB, *alpha, *beta, *D1, *D2, *D3, *loc1, *loc2, *loc3, *outloc);

	return;
}

void CheckLogLvl3(char* filename, ATC_p predef_control_values, char TransA, char TransB, double alpha, double beta, long int D1, long int D2, long int D3, int loc1, int loc2, int loc3, int outloc){
	FILE* fp = fopen(filename,"r");
	if (!fp) {
		fp = fopen(filename,"w+");
		if (!fp) error("CheckLogLvl3: LogFile %s failed to open\n", filename);
		else warning("CheckLogLvl3: Generating Logfile %s...\n", filename);
	}
	char buffer[1024], search_string[1024];
	const char* control_str = (predef_control_values) ? predef_control_values->print_csv() : "-1,-1,-1,-1";
	sprintf(search_string, "%s, %c,%c,%.5lf,%.5lf,%zu,%zu,%zu,%d,%d,%d,%d", control_str, TransA, TransB, alpha, beta, D1, D2, D3, loc1, loc2, loc3, outloc);
	while (fgets(buffer, sizeof(buffer), fp) != NULL){
		if(strstr(buffer, search_string) != NULL){
   			fprintf(stderr,"CheckLogLvl3: entry %s, %c,%c,%.5lf,%.5lf,%zu,%zu,%zu,%d,%d,%d,%d found. Quiting...\n", control_str, TransA, TransB, alpha, beta, D1, D2, D3, loc1, loc2, loc3, outloc);
			fclose(fp);
			exit(1);
		}
	}

    	fclose(fp);
	return;
}

void StoreLogLvl3(char* filename, ATC_p predef_control_values, char TransA, char TransB, double alpha, double beta, long int D1, long int D2, long int D3, int loc1, int loc2, int loc3, int outloc, double timer, double pred_t, double pred_J){
	FILE* fp = fopen(filename,"a");
	if (!fp) error("report_results: LogFile failed to open");
	const char* control_str = (predef_control_values) ? predef_control_values->print_csv() : "-1,-1,-1,-1";
   	fprintf(fp,"%s, %c,%c,%.5lf,%.5lf,%zu,%zu,%zu,%d,%d,%d,%d, %e,%e,%e\n",  control_str, TransA, TransB, alpha, beta, D1, D2, D3, loc1, loc2, loc3, outloc, timer, pred_t, pred_J);

        fclose(fp);
	return;
}
