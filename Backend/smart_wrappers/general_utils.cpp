///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief Some convinient C/C++ utilities for CHL.
///

#include <float.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <omp.h>
#include <math.h>

#include "chl_smart_wrappers.hpp"

int MAX_BACKEND_L = MAX_BACKEND_L_IN;
int STREAMING_BUFFER_OVERLAP = STREAMING_BUFFER_OVERLAP_IN;
const char* OUTPUT_ALGO_MODE = OUTPUT_ALGO_MODE_IN;
const char* DISTRIBUTION = DISTRIBUTION_IN;
const char* ORDER_2DBC = ORDER_2DBC_IN;
const char* TASK_ORDER = TASK_ORDER_IN;
const char* FETCH_ROUTING = FETCH_ROUTING_IN;
const char* WB_ROUTING = WB_ROUTING_IN;
const char* PREDICT_OPTIMIZE_TARGET = PREDICT_OPTIMIZE_TARGET_IN;

double csecond(void) {
  struct timespec tms;

  if (clock_gettime(CLOCK_REALTIME, &tms)) {
    return (0.0);
  }
  /// seconds, multiplied with 1 million
  int64_t micros = tms.tv_sec * 1000000;
  /// Add full microseconds
  micros += tms.tv_nsec / 1000;
  /// round up if necessary
  if (tms.tv_nsec % 1000 >= 500) {
    ++micros;
  }
  return ((double)micros / 1000000.0);
}

int gcd (int a, int b, int c)
{
  for(int i = std::min(a, std::min(b, c)); i>1; i--) if(( a%i == 0 ) and ( b%i == 0 ) and ( c%i == 0 )) return i;
  return 1;
}

long factorial(const int n)
{
    long f = 1;
    for (int i=1; i<=n; ++i)
        f *= i;
    return f;
}

void tab_print(int lvl){
	for (int rep=0;rep<lvl;rep++) fprintf(stderr, "\t");
}

void _printf(const char *fmt, va_list ap)
{
    if (fmt) vfprintf(stderr, fmt, ap);
    //putc('\n', stderr);
}

void warning(const char *fmt, ...) {
//#ifdef DEBUG
	fprintf(stderr, "WARNING -> ");
	va_list ap;
	va_start(ap, fmt);
	_printf(fmt, ap);
	va_end(ap);
//#endif
}

void error(const char *fmt, ...) {
	fprintf(stderr, "ERROR ->");
	va_list ap;
	va_start(ap, fmt);
	_printf(fmt, ap);
	va_end(ap);
	exit(1);
}

void massert(bool condi, const char *fmt, ...) {
	if (!condi) {
		va_list ap;
		va_start(ap, fmt);
		_printf(fmt, ap);
		va_end(ap);
		exit(1);
  	}
}

void lprintf(short lvl, const char *fmt, ...){
	tab_print(lvl);
	va_list ap;
	va_start(ap, fmt);
	_printf(fmt, ap);
	va_end(ap);
}

const char *print_mem(mem_layout mem) {
  if (mem == ROW_MAJOR)
    return "Row major";
  else if (mem == COL_MAJOR)
    return "Col major";
  else
    return "ERROR";
}

template<typename VALUETYPE>
const char *printlist(VALUETYPE *list, int length)
{
	char* outstring = (char*) malloc(abs(length)*20*sizeof(char));
	std::string printfCmd(" ");
	sprintf(outstring, "[");
	if (std::is_same<VALUETYPE, short>::value) printfCmd += "%hd";
	else if (std::is_same<VALUETYPE, int>::value) printfCmd += "%d";
	else if (std::is_same<VALUETYPE, long>::value) printfCmd += "%ld";
	else if (std::is_same<VALUETYPE, long long>::value) printfCmd += "%lld";
	else if (std::is_same<VALUETYPE, float>::value) printfCmd += "%3.3f";
	else if (std::is_same<VALUETYPE, double>::value) printfCmd += "%3.3lf";
	else if (std::is_same<VALUETYPE, void*>::value) printfCmd += "%p";
	for (int i =0; i < length; i++) sprintf(outstring + strlen(outstring), printfCmd.c_str(), list[i]);
	sprintf(outstring + strlen(outstring), " ]");
	return outstring;
}

template const char *printlist<double>(double *list, int length);
template const char *printlist<float>(float *list, int length);
template const char *printlist<int>(int *list, int length);
template const char *printlist<short>(short *list, int length);
template const char *printlist<long>(long *list, int length);
template const char *printlist<long long>(long long *list, int length);
template const char *printlist<void*>(void** list, int length);



double dabs(double x){
	if (x < 0) return -x;
	else return x;
}

inline float Serror(float a, float b) {
  if (a == 0) return (float) dabs((float)(a - b));
  return (float) dabs(a - b)/a;
}

inline double Derror(double a, double b) {
  if (a == 0) return dabs(a - b);
  return dabs(a - b)/a;
}

long int Dvec_diff(double* a, double* b, long long size, double eps) {
	long int failed = 0;
	#pragma omp parallel for
	for (long long i = 0; i < size; i++)
		if (Derror(a[i], b[i]) > eps){
			#pragma omp atomic
			failed++;
		}
	return failed;
}

long int Svec_diff(float* a, float* b, long long size, float eps) {
  	long int failed = 0;
	//#pragma omp parallel for
  	for (long long i = 0; i < size; i++)
		if (Serror(a[i], b[i]) > eps){
			//#pragma omp atomic
			failed++;
		}
  	return failed;
}

short Stest_equality(float* C_comp, float* C, long long size) {
  long int acc = 4, failed;
  float eps = 1e-4;
  failed = Svec_diff(C_comp, C, size, eps);
  while (eps > FLT_MIN && !failed && acc < 30) {
    eps *= 0.1;
    acc++;
    failed = Svec_diff(C_comp, C, size, eps);
  }
  if (4==acc) {
  	fprintf(stderr, "Test failed %zu times\n", failed);
  	int ctr = 0, itt = 0;
  	while (ctr < 10 & itt < size){
  		if (Serror(C_comp[itt], C[itt]) > eps){
  			fprintf(stderr, "Baseline vs Tested: %.10f vs %.10f\n", C_comp[itt], C[itt]);
  			ctr++;
  		}
  		itt++;
  	}
    return 0;
  } else
    fprintf(stderr, "Test passed(Accuracy= %zu digits, %zu/%lld breaking for %zu)\n\n",
            acc, failed, size, acc + 1);
  return (short) acc;
}

short Dtest_equality(double* C_comp, double* C, long long size) {
  long int acc = 8, failed;
  double eps = 1e-8;
  failed = Dvec_diff(C_comp, C, size, eps);
  while (eps > DBL_MIN && !failed && acc < 30) {
    eps *= 0.1;
    acc++;
    failed = Dvec_diff(C_comp, C, size, eps);
  }
  if (8==acc) {
  	fprintf(stderr, "Test failed %zu times\n", failed);
  	int ctr = 0;
    long long itt = 0;
  	while (ctr < 10 & itt < size){
  		if (Derror(C_comp[itt], C[itt]) > eps){
  			fprintf(stderr, "Baseline vs Tested(adr = %p, itt = %lld): %.15lf vs %.15lf\n", &C[itt], itt, C_comp[itt], C[itt]);
  			ctr++;
  		}
  		itt++;
  	}
  return 0;
  }
  else
    fprintf(stderr, "Test passed(Accuracy= %zu digits, %zu/%lld breaking for %zu)\n\n",
            acc, failed, size, acc + 1);
  return (short) acc;
}

inline double Herror(__half a, __half b) {
  if (a == 0) return (__half) dabs(a - b);
  return (__half) dabs(a - b)/a;
}

long int Hvec_diff(__half* a, __half* b, long long size, __half eps) {
	long int failed = 0;
	#pragma omp parallel for
	for (long long i = 0; i < size; i++)
		if (Herror(a[i], b[i]) > eps){
			#pragma omp atomic
			failed++;
		}
	return failed;
}

short Htest_equality(__half* C_comp, __half* C, long long size) {
  long int acc = 2, failed;
  float eps = 1e-2;
  failed = Hvec_diff(C_comp, C, size, eps);
  while (eps > FLT_MIN && !failed && acc < 30) {
    eps *= 0.1;
    acc++;
    failed = Hvec_diff(C_comp, C, size, eps);
  }
  if (2==acc) {
  	fprintf(stderr, "Test failed %zu times\n", failed);
  	int ctr = 0, itt = 0;
  	while (ctr < 10 & itt < size){
  		if (Herror(C_comp[itt], C[itt]) > eps){
  			fprintf(stderr, "Baseline vs Tested: %lf vs %lf\n", (double) C_comp[itt], (double) C[itt]);
  			ctr++;
  		}
  		itt++;
  	}
    return 0;
  } else
    fprintf(stderr, "Test passed(Accuracy= %zu digits, %zu/%lld breaking for %zu)\n\n",
            acc, failed, size, acc + 1);
  return (short) acc;
}

long int count_lines(FILE* fp){
	if (!fp) error("count_lines: fp = 0 ");
	int ctr = 0;
	char chr = getc(fp);
	while (chr != EOF){
		//Count whenever new line is encountered
		if (chr == '\n') ctr++;
		//take next character from file.
		chr = getc(fp);
	}
	fseek(fp, 0, SEEK_SET);;
	return ctr;
}

int check_benchmark(char *filename){
	FILE* fp = fopen(filename,"r");
	if (!fp) {
		fp = fopen(filename,"w+");
		if (!fp) error("check_benchmark: LogFile %s failed to open\n", filename);
		else warning("Generating Logfile %s...\n", filename);
		fclose(fp);
	}
	else {
		fprintf(stderr,"Benchmark found: %s\n", filename);
		fclose(fp);
		return 1;
	}
	return 0;
}

double Gval_per_s(long long value, double time){
  return value / (time * 1e9);
}

long long gemm_ops(long int M, long int N, long int K){
	return (long long) M * N * (2 * K + 1);
}

long long gemm_mem_ops(long int M, long int N, long int K){
	return M * K + K * N + 2*(M * N);
}

long long gemm_mem_sz(long int M, long int N, long int K, int elemSize){
	return (M * K + K * N + M * N)*elemSize;
}

/*
long long gemv_flops(long int M, long int N){
	return (long long) M * (2 * N + 1);
}

long long gemv_memory(long int M, long int N, long int A_loc, long int x_loc, long int y_loc, short dsize){
	return (M * N * A_loc + N * x_loc + M * y_loc)*dsize;
}

long long axpy_flops(long int  N){
	return (long long) 2* N;
}

long long axpy_memory(long int  N, long int x_loc, long int y_loc, short dsize){
	return (long long) N*(x_loc + y_loc)*dsize;
}

long long dot_flops(long int  N){
	return (long long) 2* N;
}

long long dot_memory(long int  N, long int x_loc, long int y_loc, short dsize){
	return (long long) N*(x_loc + y_loc)*dsize;
}*/

void DECOM_2D(int val, int* valgrid_1_ptr, int* valgrid_2_ptr){
	// 2D Block cyclic device decomposition
	int D1_parts, D2_parts; 
	D1_parts = std::sqrt(val);
	D2_parts = D1_parts;
	if (D1_parts ==0) { D2_parts = val; D1_parts = 1; }
	else {
		// find the most square decomposition of val in D1_parts x D2_parts
		int g;
		for (g = D1_parts+1; g>0; --g)
		if (val % g == 0) break;
		if (g==0) { D1_parts = val; D2_parts = 1; }
		//if (g==0) { D1_parts = 1; D2_parts = val; }
		else { D1_parts = g; D2_parts = val/g; }
	}
	/// If ORDER_2DBC="D2_lesseq_D1", reverse layout. 
	if (!strcmp(ORDER_2DBC, "D2_lesseq_D1")){
		int tmp = D1_parts;
		D1_parts = D2_parts;
		D2_parts = tmp;
	}
	*valgrid_1_ptr = D1_parts;
	*valgrid_2_ptr = D2_parts;
}


int is_in_list(int elem, int* elem_list, int list_len){ 
	for (int idx = 0; idx < list_len; idx++)
		if(elem_list[idx] == elem) return 1; 
	return 0; 
}

int get_loc_in_list(int elem, int* elem_list, int list_len){ 
	for (int idx = 0; idx < list_len; idx++)
		if(elem_list[idx] == elem) return idx; 
	return -1; 
}


void translate_binary_to_unit_list(int case_id, int* active_unit_num_p, int* active_unit_id_list){
	int mask;
	*active_unit_num_p = 0;
	for (int mask_offset = 0; mask_offset < CHL_WORKERS; mask_offset++){
		mask =  1 << mask_offset;
		if (case_id & mask){
			active_unit_id_list[*active_unit_num_p] = mask_offset;
			(*active_unit_num_p)++;
		}
	}
}

int translate_unit_list_to_binary(int* active_unit_id_list, int active_unit_num){
	int case_id_out = 0; 
	for (int mask_offset = 0; mask_offset < CHL_WORKERS; mask_offset++)
		if (is_in_list(mask_offset, active_unit_id_list, active_unit_num)) 
			case_id_out+= pow(2, mask_offset);
	return case_id_out;
}

int binary_case_id_split(int case_id){
    int active_unit_id_list[CHL_WORKERS], active_unit_num;
    translate_binary_to_unit_list(case_id, &active_unit_num, active_unit_id_list);
	int out_active_unit_id_list[CHL_WORKERS], out_active_unit_num = 0;
	for(int idx = 0; idx < active_unit_num; idx+=2)
		out_active_unit_id_list[out_active_unit_num++] = active_unit_id_list[idx];
	if (!(out_active_unit_num*2 == active_unit_num)) 
		error("binary_case_id_split(%d): active_unit_num(=%d)%%2 is not 0\n", case_id, active_unit_num);
	//fprintf(stderr, "active_unit_id_list = %s, out_active_unit_id_list = %s\n", 
	//	printlist(active_unit_id_list, active_unit_num), printlist(out_active_unit_id_list, out_active_unit_num));
	return translate_unit_list_to_binary(out_active_unit_id_list, out_active_unit_num);
}

int is_subset(int case_id, int case_id_set){
	int mask;
	for (int mask_offset = 0; mask_offset < CHL_WORKERS; mask_offset++){
		mask =  1 << mask_offset;
		if (!(case_id_set & mask) && (case_id & mask)) return 0;
	}
	return 1;
}

int normal_equal(double n1, double n2){
	if(n1*(1+NORMALIZE_NEAR_SPLIT_LIMIT)>=n2 && n1*(1-NORMALIZE_NEAR_SPLIT_LIMIT)<=n2) return 1;
	return 0; 
}

int normal_less(double n1, double n2){
	if(n1*(1+NORMALIZE_NEAR_SPLIT_LIMIT)<n2) return 1;
	return 0; 
}

int normal_larger(double n1, double n2){
	if(n1*(1-NORMALIZE_NEAR_SPLIT_LIMIT)>n2) return 1;
	return 0; 
}
int normal_lessequal(double n1, double n2){
	return normal_equal(n1, n2) + normal_less(n1, n2); 
}
int normal_largerequal(double n1, double n2){
	return normal_equal(n1, n2) + normal_larger(n1, n2); 
}

void CHLSetTimerAsync(void* wrapped_timer_Ptr){
  double* timer = (double*) wrapped_timer_Ptr;
  *timer = csecond();
#ifdef DEBUG
  lprintf(6, "CHLSetTimerAsync(%p) ran succesfully.\n", wrapped_timer_Ptr);
#endif
}

template<typename VALUETYPE> void swap(VALUETYPE* elem1, VALUETYPE* elem2){
	VALUETYPE temp = *elem1; 
	*elem1 = *elem2;
	*elem2 = temp;
}

template void swap<int>(int* elem1, int* elem2);
template void swap<double>(double* elem1, double* elem2);

void bubble_sort_3_vals(int* case_ids, double vals[][3], double weight[3], int len){
	int i, j;
	for (i = 0; i < len - 1; i++) for (j = 0; j < len - i - 1; j++) {
		if ((weight[0]*vals[j][0]+ weight[1]*vals[j][1] + weight[2]*vals[j][2]) <
				(weight[0]*vals[j+1][0]+ weight[1]*vals[j+1][1] + weight[2]*vals[j+1][2])) {
			swap<int>(&(case_ids[j]), &(case_ids[j + 1]));
			swap<double>(&(vals[j][0]), &(vals[j+1][0]));
			swap<double>(&(vals[j][1]), &(vals[j+1][1]));
			swap<double>(&(vals[j][2]), &(vals[j+1][2]));
		}
	}
}

void bubble_sort_2_idx_3x3_vals(int* case_ids, int* case_ids_2, double vals[][3], double vals_2[][3], double vals_3[][3], double weight[3], int len){
	int i, j;
	for (i = 0; i < len - 1; i++) for (j = 0; j < len - i - 1; j++) {
		if ((weight[0]*vals[j][0]+ weight[1]*vals[j][1] + weight[2]*vals[j][2]) <
				(weight[0]*vals[j+1][0]+ weight[1]*vals[j+1][1] + weight[2]*vals[j+1][2])) {
			swap<int>(&(case_ids[j]), &(case_ids[j + 1]));
			swap<int>(&(case_ids_2[j]), &(case_ids_2[j + 1]));
			for (int idx = 0; idx < 3; idx++){
				swap<double>(&(vals[j][idx]), &(vals[j+1][idx]));
				swap<double>(&(vals_2[j][idx]), &(vals_2[j+1][idx]));
				swap<double>(&(vals_3[j][idx]), &(vals_3[j+1][idx]));
			}
		}
	}
}

void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}
