///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief The start of Zawarudo
///

#include "PARALiA.hpp"
#include "BackenedLibsWrapped.hpp"
#include "Testing.hpp"

#include "backend_wrappers.hpp"

#define CBLASXT_MAX_SAFE_TILE 8192

int main(const int argc, const char *argv[]) {
	char TransA, TransB;
  	double alpha_in, beta_in;
	long int M, N, K;
	int A_loc, B_loc, C_loc, C_out_loc;
	ATC_p predef_control_values = NULL, return_values = NULL;
	ParseInputLvl3(argc, argv, &predef_control_values, &TransA, &TransB, &alpha_in, &beta_in, &M, &N, &K, &A_loc, &B_loc, &C_loc, &C_out_loc);
  	float alpha = alpha_in, beta = beta_in;

	char *filename = (char *) malloc(2048* sizeof(char));
	if (predef_control_values!= NULL){
		if(predef_control_values->T > 0) {
			if (predef_control_values->T > M || predef_control_values->T > N || predef_control_values->T > K)
				error("Given Tin=%ld bigger than problem dim\n", predef_control_values->T);
			else if (predef_control_values->T > M/1.5 && predef_control_values->T > N/1.5 && predef_control_values->T > K/1.5)
				warning("Given Tin=%ld bigger than all problem dims/1.5\n", predef_control_values->T);
		}
		sprintf(filename, "%s/sgemm_runner_predefined_vals_%s_%s_strong_scaling_bench.log",
			TESTLIBDIR, CoCoImplementationPrint(), VERSION);
	}
	else sprintf(filename, "%s/sgemm_runner_%s_%s_strong_scaling_bench.log",
		TESTLIBDIR, CoCoImplementationPrint(), VERSION);
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

	fprintf(stderr, "\nAllocating memory...\n");
	//for(int d=0; d < CHL_MEMLOCS; d++) CHLEnableLinks(d, CHL_MEMLOCS);

	float *A, *B, *C;
	// allocate in device if loc = 0, otherwise allocate in pinned memory for benchmarks
	A = (float*) CHLMalloc(M * K*sizeof(float), A_loc, 0);
	B = (float*) CHLMalloc(N * K*sizeof(float), B_loc, 0);
	C = (float*) CHLMalloc(M * N*sizeof(float), C_loc, 1);
	if(A_loc == CHL_WORKERS) CHLTouche(A, M*K, sizeof(float));
 	if(B_loc == CHL_WORKERS) CHLTouche(B, N*K, sizeof(float));
	if(C_loc == CHL_WORKERS) CHLTouche(C, M*N, sizeof(float));
	CHLSyncCheckErr();
	cpu_timer  = csecond() - cpu_timer;
	fprintf(stderr, "Done: Alloc time:\t%lf ms\n\n",  cpu_timer  * 1000);
	cpu_timer = csecond();
	fprintf(stderr, "Initializing to random values (VALIDATE)...");
	CHLVecInit(A, K * M, 42, A_loc);
	CHLVecInit(B, K * N, 43, B_loc);
	CHLVecInit(C, M * N, 44, C_loc);
	CHLSyncCheckErr();
	cpu_timer  = csecond() - cpu_timer ;
	fprintf(stderr, "done.\nInit time:\t%lf ms\n\n",  cpu_timer  * 1000);

#ifdef RUNVALIDATION
	float *C_out, *C_out1, *C_buf;
	C_out  = (float*) CHLMalloc(M * N*sizeof(float), CHL_MEMLOCS-1, 1);
	C_out1  = (float*) CHLMalloc(M * N*sizeof(float), CHL_MEMLOCS-1, 1);
	C_buf  = (float*) CHLMalloc(M * N*sizeof(float), CHL_MEMLOCS-1, 1);

	CHLMemcpy(C_buf, C,  M * N *sizeof(float), CHL_MEMLOCS -1, C_loc);

	// Call for Validate start
	if (predef_control_values!= NULL) return_values = PARALiASgemmControled(TransA, TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C , ldC, predef_control_values);
	else return_values = PARALiASgemm(TransA, TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C , ldC);
	CHLSyncCheckErr();
	CHLMemcpy(C, C_buf,  M * N *sizeof(float), C_loc, CHL_MEMLOCS -1);

	// Call for Validate reuse
	if (predef_control_values!= NULL) return_values = PARALiASgemmControled(TransA, TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C , ldC, predef_control_values);
	else return_values = PARALiASgemm(TransA, TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C , ldC);
	CHLSyncCheckErr();
	for (int i = 0; i< CHL_MEMLOCS; i++) PARALiADevCacheFree(i);

 	CHLMemcpy(C_out, C,  M * N *sizeof(float), CHL_MEMLOCS -1, C_loc);
 	CHLMemcpy(C, C_buf,  M * N *sizeof(float), C_loc, CHL_MEMLOCS -1);

	// Validate with cuBLASXt (questionable but CPU validation can be slower by at least a factor)
	int dev_ids[CHL_WORKERS];
	for (int i = 0; i < CHL_WORKERS; i++) dev_ids[i] = i;
	cuBLASXtSgemmWrap(TransA,  TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C, ldC,  (long int) fmin(fmin(fmin(M,N),K)/2,CBLASXT_MAX_SAFE_TILE), 0, CHL_WORKERS, dev_ids);
	CHLMemcpy(C_out1, C,  M * N *sizeof(float), CHL_MEMLOCS -1, C_loc);
 	if(Stest_equality(C_out1, C_out, M * N) < 4) error("Insufficient accuracy for benchmarks\n");

 	CHLMemcpy(C, C_buf,  M * N *sizeof(float), C_loc, CHL_MEMLOCS -1);
	CHLFree(C_out, M * N*sizeof(float), CHL_MEMLOCS-1);
	CHLFree(C_out1, M * N*sizeof(float), CHL_MEMLOCS-1);
	CHLFree(C_buf, M * N*sizeof(float), CHL_MEMLOCS-1);
#endif

	cpu_timer = csecond();
	if (predef_control_values!= NULL) return_values = PARALiASgemmControled(TransA, TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C , ldC, predef_control_values);
	else return_values = PARALiASgemm(TransA, TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C , ldC);
	CHLSyncCheckErr();
	cpu_timer  = csecond() - cpu_timer;

#ifdef CHECKLOG
	CheckLogLvl3(filename, return_values, TransA, TransB, alpha, beta, M, N, K, A_loc, B_loc, C_loc, C_out_loc);
#endif
	// Store the time required for the first call (+ 1-time overheads etc)
	//StoreLogLvl3(filename, return_values, TransA, TransB, alpha, beta, M, N, K, A_loc, B_loc, C_loc, C_out_loc, cpu_timer);

	double first_over_t = cpu_timer;

	int warmup_bench_it = REP_TILE*2 + 1;
	for(int it = 0; it < warmup_bench_it; it++){
		if (predef_control_values!= NULL) return_values = PARALiASgemmControled(TransA, TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C , ldC, predef_control_values);
		else return_values = PARALiASgemm(TransA, TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C , ldC);
	}
	CHLSyncCheckErr();

	double min_t = first_over_t, max_t = 0, avg_t = 0;
	cpu_timer = csecond();
	int bench_it = 100;
	//TODO: bench if ( M >= 20000 && N >= 20000 && K >= 20000) bench_it = 20;
	bench_it = 10;
	for(int it = 0; it < bench_it; it++){
		cpu_timer = csecond();
		if (predef_control_values!= NULL) return_values = PARALiASgemmControled(TransA, TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C , ldC, predef_control_values);
		else return_values = PARALiASgemm(TransA, TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C , ldC);
		CHLSyncCheckErr();
		cpu_timer = csecond() - cpu_timer;
		StoreLogLvl3(filename, return_values, TransA, TransB, alpha, beta, M, N, K, A_loc, B_loc, C_loc, C_out_loc, cpu_timer, return_values->pred_t, return_values->pred_J);
		if ( cpu_timer < min_t ) min_t = cpu_timer;
		if ( cpu_timer > max_t ) max_t = cpu_timer;
		avg_t += cpu_timer;
	}
	avg_t/=bench_it;
	fprintf(stderr, "sgemm_runner (%s):\n\tfirst_it_t = %lf ms ( %lf Gflops/s )\n\tavg_t = %lf ms ( %lf Gflops/s )\n\tmin_t = %lf ms ( %lf Gflops/s )\n\tmax_t = %lf ms ( %lf Gflops/s )\n",
	return_values->print_csv(),
	first_over_t  * 1000, Gval_per_s(gemm_ops(M,N,K),first_over_t),
	avg_t  * 1000, Gval_per_s(gemm_ops(M,N,K),avg_t),
	min_t  * 1000, Gval_per_s(gemm_ops(M,N,K),min_t),
	max_t  * 1000, Gval_per_s(gemm_ops(M,N,K),max_t));

	for (int i = 0; i< CHL_MEMLOCS; i++) PARALiADevCacheFree((i));

	CHLSyncCheckErr();
	CHLFree(A, M * K*sizeof(float), A_loc);
	CHLFree(B, N * K*sizeof(float), B_loc);
	CHLFree(C, M * N*sizeof(float), C_loc);
	return 0;
}
