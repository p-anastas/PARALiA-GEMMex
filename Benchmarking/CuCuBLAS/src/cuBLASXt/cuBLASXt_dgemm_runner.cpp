///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief The start of Zawarudo
///

#include "backend_wrappers.hpp"
#include "PARALiA.hpp"
#include "BackenedLibsWrapped.hpp"
#include "Testing.hpp"

#define CBLASXT_MAX_SAFE_TILE 8192

int main(const int argc, const char *argv[]) {

	char TransA, TransB;
  	double alpha, beta;
	long int M, N, K;
	int A_loc, B_loc, C_loc, C_out_loc;

	ATC_p predef_control_values = NULL;
	ParseInputLvl3(argc, argv, &predef_control_values, &TransA, &TransB, &alpha, &beta, &M, &N, &K, &A_loc, &B_loc, &C_loc, &C_out_loc);

	char *filename = (char *) malloc(1024* sizeof(char));
	if (predef_control_values!= NULL){
		if(predef_control_values->T > 0) {
			if (predef_control_values->T > M || predef_control_values->T > N || predef_control_values->T > K)
				error("Given Tin=%ld bigger than problem dim\n", predef_control_values->T);
			else if (predef_control_values->T > CBLASXT_MAX_SAFE_TILE)
				error("Given Tin=%ld bigger than CBLASXT_MAX_SAFE_TILE\n", predef_control_values->T);
		}
		sprintf(filename, "%s/cuBLASXtDgemmRunner_predefined_vals_%s.log",
			TESTLIBDIR, VERSION);
	}
	else sprintf(filename, "%s/cuBLASXtDgemmRunner_%s.log", TESTLIBDIR, VERSION);
#ifdef CHECKLOG
	CheckLogLvl3(filename, predef_control_values, TransA, TransB, alpha, beta, M, N, K, A_loc, B_loc, C_loc, C_out_loc);
#endif
	long int cublasXt_tile;
	if (predef_control_values!= NULL && predef_control_values->T > 0) cublasXt_tile = predef_control_values->T;
	else cublasXt_tile = (long int) fmin(fmin(fmin(M,N),K)/2,CBLASXT_MAX_SAFE_TILE);
	double cache_limit;
	if (predef_control_values!= NULL && predef_control_values->cache_limit > 0) cache_limit = predef_control_values->cache_limit;
	else cache_limit = -1;
	int dev_num, *dev_ids;
	if (predef_control_values!= NULL && predef_control_values->active_unit_num > 0){
		dev_num = predef_control_values->active_unit_num;
		dev_ids = (int*) malloc(dev_num*sizeof(int));
		for(int idx =0; idx < predef_control_values->active_unit_num; idx++)
			dev_ids[idx] = predef_control_values->active_unit_id_list[idx];
	}
	else{
		dev_num = CHL_WORKERS;
		dev_ids = (int*) malloc(dev_num*sizeof(int));
		for (int i = 0; i < dev_num; i++) dev_ids[i] = i;
	}
	if(!predef_control_values) predef_control_values = new ATC();
#ifdef CHECKLOG
	CheckLogLvl3(filename, predef_control_values, TransA, TransB, alpha, beta, M, N, K, A_loc, B_loc, C_loc, C_out_loc);
#endif
	/// Matrix Layouts for CPU GEMM
	CBLAS_TRANSPOSE cpu_op_A, cpu_op_B;    // CblasNoTrans, CblasTrans
	cublasOperation_t gpu_op_A, gpu_op_B; // CUBLAS_OP_N, CUBLAS_OP_T

	long int ldA, ldB, ldC = M;
	TransposeTranslate(TransA, &cpu_op_A, &gpu_op_A, &ldA, M, K);
	TransposeTranslate(TransB, &cpu_op_B, &gpu_op_B, &ldB, K, N);

	/// Local Timers
	double cpu_timer = csecond();

	fprintf(stderr, "\nAllocating memory...");

	double *A, *B, *C;
	// allocate in device if loc = 0, otherwise allocate in pinned memory for benchmarks
	A = (double*) CHLMalloc(M * K*sizeof(double), A_loc, 0);
	B = (double*) CHLMalloc(N * K*sizeof(double), B_loc, 0);
	C = (double*) CHLMalloc(M * N*sizeof(double), C_loc, 1);
	//if(A_loc >= CHL_WORKERS) CHLTouche(A, M*K, sizeof(double));
 	//if(B_loc >= CHL_WORKERS) CHLTouche(B, N*K, sizeof(double));
 	//if(C_loc >= CHL_WORKERS) CHLTouche(C, M*N, sizeof(double));

	CHLSyncCheckErr();
	cpu_timer  = csecond() - cpu_timer;
	fprintf(stderr, "done.\nAlloc time:\t%lf ms\n\n",  cpu_timer  * 1000);
	cpu_timer = csecond();
	fprintf(stderr, "Initializing to random values (VALIDATE)...");
	CHLVecInit(A, K * M, 42, A_loc);
	CHLVecInit(B, K * N, 43, B_loc);
	CHLVecInit(C, M * N, 44, C_loc);
	CHLSyncCheckErr();
	cpu_timer  = csecond() - cpu_timer ;
	fprintf(stderr, "done.\nInit time:\t%lf ms\n\n",  cpu_timer  * 1000);

	// First call with T set to half the smaller problem dim (unless predefined or larger than CBLASXT_MAX_SAFE_TILE)
	cpu_timer = csecond();
	cuBLASXtDgemmWrap(TransA,  TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C, ldC,  cublasXt_tile, 0, dev_num, dev_ids);
	CHLSyncCheckErr();
	cpu_timer  = csecond() - cpu_timer;
	double best_standard_tile_t = cpu_timer;

	if (!(predef_control_values!= NULL && predef_control_values->T > 0)){
		// Second call with T set to smaller problem dim ( can be better for small problems with fat/thin matrices)
		long int cublasXt_min_dim = (long int) fmin(fmin(fmin(M,N),K),CBLASXT_MAX_SAFE_TILE);
		cpu_timer = csecond();
		cuBLASXtDgemmWrap(TransA,  TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C, ldC,  cublasXt_min_dim, 0, dev_num, dev_ids);
		CHLSyncCheckErr();
		cpu_timer  = csecond() - cpu_timer;
		if ( cpu_timer < best_standard_tile_t) {
			cublasXt_tile = cublasXt_min_dim;
			best_standard_tile_t = cpu_timer;
		}

	}
#ifdef CHECKLOG
	CheckLogLvl3(filename, predef_control_values, TransA, TransB, alpha, beta, M, N, K, A_loc, B_loc, C_loc, C_out_loc);
#endif
	double cublasXt_t = best_standard_tile_t;
	if (predef_control_values!= NULL && predef_control_values->T > 0){
		fprintf(stderr,"Running CUBLASXT DGEMM-> M = %zu, N = %zu, K = %zu, T = %zu\n", M, N, K, cublasXt_tile);
		cpu_timer  = csecond();
		cuBLASXtDgemmWrap(TransA, TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C, ldC,  cublasXt_tile, 0, dev_num, dev_ids);
		CHLSyncCheckErr();
		cpu_timer  = csecond() - cpu_timer;
		fprintf(stderr, "Total time:\t%lf ms\n", cpu_timer  * 1000);
		if (cublasXt_t > cpu_timer) cublasXt_t = cpu_timer;
	}
	else {
		for (long int T_trial = (((long int)fmax(fmin(fmin(M/8,N/8),K/8),1024))/1024)*1024; T_trial <= fmin(fmin(fmin(M,N),K),CBLASXT_MAX_SAFE_TILE); T_trial+=1024) if (M >= T_trial*1.5 || N >= T_trial*1.5 || K >= T_trial*1.5){
			fprintf(stderr,"Running CUBLASXT DGEMM-> M = %zu, N = %zu, K = %zu, T = %zu\n", M, N, K, T_trial);
			cpu_timer  = csecond();
			cuBLASXtDgemmWrap(TransA,  TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C, ldC,  T_trial, 0, dev_num, dev_ids);
			CHLSyncCheckErr();
			cpu_timer  = csecond() - cpu_timer;
			fprintf(stderr, "Total time:\t%lf ms\n", cpu_timer  * 1000);
			if (cpu_timer < cublasXt_t){
				cublasXt_t = cpu_timer;
				cublasXt_tile = T_trial;
			}
		}
		fprintf(stderr, "\nCUBLASXT DGEMM T_best = %zu : t = %lf ms ( %lf Gflops/s )\n\n", cublasXt_tile, cublasXt_t  * 1000, Gval_per_s(gemm_ops(M,N,K),cublasXt_t));
	}
	predef_control_values->T = cublasXt_tile;
	fprintf(stderr,"Running CUBLASXT DGEMM-> M = %zu, N = %zu, K = %zu T_best = %zu\n", M, N, K, cublasXt_tile);
#ifdef CHECKLOG
	CheckLogLvl3(filename, predef_control_values, TransA, TransB, alpha, beta, M, N, K, A_loc, B_loc, C_loc, C_out_loc);
#endif
	double min_t = cublasXt_t, max_t = 0, avg_t = 0;
	cpu_timer = csecond();
	int bench_it = 100;
	//TODO: bench if ( M >= 20000 && N >= 20000 && K >= 20000) bench_it = 20;
	bench_it = 10;
	for(int it = 0; it < bench_it; it++){
		cpu_timer = csecond();
		cuBLASXtDgemmWrap(TransA,  TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C, ldC,  cublasXt_tile, 0, dev_num, dev_ids);
		CHLSyncCheckErr();
		cpu_timer = csecond() - cpu_timer;
		StoreLogLvl3(filename, predef_control_values, TransA, TransB, alpha, beta, M, N, K, A_loc, B_loc, C_loc, C_out_loc, cpu_timer, -1, -1);
		if ( cpu_timer < min_t ) min_t = cpu_timer;
		if ( cpu_timer > max_t ) max_t = cpu_timer;
		avg_t += cpu_timer;
	}
	avg_t/=bench_it;
	fprintf(stderr, "cuBLASXt (%s):\n\tavg_t = %lf ms ( %lf Gflops/s )\n\tmin_t = %lf ms ( %lf Gflops/s )\n\tmax_t = %lf ms ( %lf Gflops/s )\n",
	predef_control_values->print_csv(),
	avg_t  * 1000, Gval_per_s(gemm_ops(M,N,K),avg_t),
	min_t  * 1000, Gval_per_s(gemm_ops(M,N,K),min_t),
	max_t  * 1000, Gval_per_s(gemm_ops(M,N,K),max_t));

	CHLSyncCheckErr();
	CHLFree(A, M * K* sizeof(double), A_loc);
	CHLFree(B, N * K* sizeof(double), B_loc);
	CHLFree(C, M * N* sizeof(double), C_loc);
	return 0;
}
