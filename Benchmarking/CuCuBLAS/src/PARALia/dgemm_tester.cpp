///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief The start of Zawarudo
///

#include "backend_wrappers.hpp"
#include "PARALiA.hpp"
#include "BackenedLibsWrapped.hpp"
#include "Testing.hpp"

#define CBLASXT_MAX_SAFE_TILE 10000

int main(const int argc, const char *argv[]) {

	int run_cpu_mem, run_gpu_mem, run_large;

	char TransA, TransB;
  	double alpha, beta;
	long int M, N, K, T;
	int A_loc, B_loc, C_loc, C_out_loc;
	double cache_limit = 0;
	long int ldA, ldB, ldC;

	if(argc != 4) error("Incorrect input arguments. Usage: ./correct_run run_cpu_mem run_gpu_mem run_large\n");
	// Control Parameters
	run_cpu_mem = atoi(argv[1]);
	run_gpu_mem = atoi(argv[2]);
	run_large = atoi(argv[3]);

	/// Local Timers
	double cpu_timer = csecond();
	fprintf(stderr, "dgemm_tester: Initallizing tests for PARALiADgemmTile with\
		run_cpu_mem = %d, run_gpu_mem = %d, run_large = %d\n", run_cpu_mem, run_gpu_mem, run_large);

	int dev_ids[CHL_WORKERS];
	for(int i = 0; i< CHL_WORKERS; i++) dev_ids[i] = i;
	double *A, *B, *C, *C_comp, *C_buf;

	ATC_p ret_autotune_val;
	ldA = ldB = ldC = M = N = K = 8192;

	fprintf(stdout, "\n-----------------------------------------Testing GEMM (dtype=double)-----------------------------------------\n");
	fprintf(stdout, "\n0) PARALiA build option string : ");
	char* dummy = CoCoImplementationPrint();
	free(dummy);

	if(run_cpu_mem){
		fprintf(stdout, "\n-----------------------------------------Testing for matrices on CPU-----------------------------------------\n");
		fprintf(stderr, "dgemm_tester: Allocating CPU buffers...->100 MB...\n");
		A_loc = B_loc = C_loc = CHL_MEMLOCS - 1;

		A = (double*) CHLMalloc(M * K*sizeof(double), A_loc, 1);
		B = (double*) CHLMalloc(N * K*sizeof(double), B_loc, 1);
		C = (double*) CHLMalloc(M * N*sizeof(double), C_loc, 1);

		CHLSyncCheckErr();
		cpu_timer  = csecond() - cpu_timer;
		fprintf(stderr, "done.\nAlloc time:\t%lf ms\n\n",  cpu_timer  * 1000);
		cpu_timer = csecond();
		fprintf(stderr, "Initializing to random values...");
		CHLVecInit(A, K * M, 42, A_loc);
		CHLVecInit(B, K * N, 43, B_loc);
		CHLVecInit(C, M * N, 44, C_loc);
		CHLSyncCheckErr();
		cpu_timer  = csecond() - cpu_timer ;
		fprintf(stderr, "done.\nInit time:\t%lf ms\n\n",  cpu_timer  * 1000);

		C_comp = (double*) CHLMalloc(M * N*sizeof(double), CHL_MEMLOCS - 1, 1);
		C_buf = (double*) CHLMalloc(M * N*sizeof(double), CHL_MEMLOCS - 1, 1);
		CHLMemcpy(C_comp, C,  M * N *sizeof(double), CHL_MEMLOCS - 1, C_loc);
		CHLMemcpy(C_buf, C,  M * N *sizeof(double), CHL_MEMLOCS - 1, C_loc);

		fprintf(stdout, "\n1) Testing Square Problems < 100 MB:\n\n");
		TransA = TransB = 'N';
		alpha = 1.23;
		beta = 0.9876;
		for (int dim = 256; dim <= M; dim*=2){
			fprintf(stdout, "[M, N, K] = [%d, %d, %d]\t-> ", dim, dim, dim);
			cpu_timer = csecond();
			ret_autotune_val = PARALiADgemm(TransA, TransB, dim, dim, dim, alpha, A, ldA, B, ldB, beta, C , ldC);
			cpu_timer  = csecond() - cpu_timer;
			double comp_flops = Gval_per_s(gemm_ops(dim,dim,dim),cpu_timer);
			fprintf(stderr, "M=N=K=%d: Gflops/s -> ", dim);
			fprintf(stderr, "PARALiA: %.1lf, ", comp_flops);
			cpu_timer = csecond();
			T = fmin(dim,fmin(dim,dim))/2;
			cuBLASXtDgemmWrap(TransA, TransB, dim, dim, dim, alpha, A, ldA, B, ldB, beta, C_comp, ldC,  T, cache_limit, CHL_WORKERS, dev_ids);
			CHLSyncCheckErr();
			cpu_timer  = csecond() - cpu_timer;
			fprintf(stderr, "cuBLASXT: %.1lf\n", Gval_per_s(gemm_ops(dim,dim,dim),cpu_timer));
			fprintf(stderr, "%s\n", ret_autotune_val->print_csv());
			if (comp_flops < Gval_per_s(gemm_ops(dim,dim,dim),cpu_timer)) warning("Inferior Perf to cublasXt\n");
			int succcess = Dtest_equality(C_comp, C, dim * dim); 
			if (succcess) fprintf(stdout, " | First run: OK");
			else fprintf(stdout, " | First run: RAN, WRONG C");
			CHLMemcpy(C, C_buf, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
			ret_autotune_val = PARALiADgemm(TransA, TransB, dim, dim, dim, alpha, A, ldA, B, ldB, beta, C , ldC);
			int succcess2 = Dtest_equality(C_comp, C, dim * dim);
			if (succcess2) fprintf(stdout, " | Metadata-reuse run: OK");
			else fprintf(stdout, "| Metadata-reuse run: RAN, WRONG C");
			if(succcess && succcess2) fprintf(stdout, " | -> Test PASSED\n");
			else fprintf(stdout, " | -> Test FAILED\n");
			CHLMemcpy(C, C_buf, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
			CHLMemcpy(C_comp, C_buf, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
			CHLSyncCheckErr();
		}

		CHLVecInit(C, M * N, 44, C_loc);
		CHLMemcpy(C_comp, C,  M * N *sizeof(double), C_loc, C_loc);

		fprintf(stdout, "\n2) Testing Non-Square Problems < 100 MB:\n\n");
		alpha = 1.23;
		beta = 0.9876;
		for (int dim1 = 256; dim1 <= M; dim1*=4) for (int dim2 = 256; dim2 <= N; dim2*=4) for (int dim3 = 256; dim3 <= K; dim3*=4) if ( dim1 != dim2 || dim2 != dim3 || dim1!= dim3){
			fprintf(stdout, "[M, N, K] = [%d, %d, %d]\t-> ", dim1, dim2, dim3);
			cpu_timer = csecond();
			ret_autotune_val = PARALiADgemm(TransA, TransB, dim1, dim2, dim3, alpha, A, ldA, B, ldB, beta, C , ldC);
			cpu_timer  = csecond() - cpu_timer;
			double comp_flops =  Gval_per_s(gemm_ops(dim1,dim2,dim3),cpu_timer);
			fprintf(stderr, "M=%d,N=%d,K=%d: Gflops/s -> ", dim1, dim2, dim3);
			fprintf(stderr, "PARALiA: %.1lf, ", comp_flops);
			cpu_timer = csecond();
			T = fmin(dim1,fmin(dim2,dim3))/2;
			cuBLASXtDgemmWrap(TransA, TransB, dim1, dim2, dim3, alpha, A, ldA, B, ldB, beta, C_comp, ldC,  T, cache_limit, CHL_WORKERS, dev_ids);
			CHLSyncCheckErr();
			cpu_timer  = csecond() - cpu_timer;
			fprintf(stderr, "cuBLASXT: %.1lf\n", Gval_per_s(gemm_ops(dim1,dim2,dim3),cpu_timer));
			fprintf(stderr, "%s\n", ret_autotune_val->print_csv());
			if (comp_flops < Gval_per_s(gemm_ops(dim1,dim2,dim3),cpu_timer)) warning("Inferior Perf to cublasXt\n");
			int succcess = Dtest_equality(C_comp, C, dim1 * dim2); 
			if (succcess) fprintf(stdout, " | First run: OK");
			else fprintf(stdout, " | First run: RAN, WRONG C");
			CHLMemcpy(C, C_buf, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
			ret_autotune_val = PARALiADgemm(TransA, TransB, dim1, dim2, dim3, alpha, A, ldA, B, ldB, beta, C , ldC);
			int succcess2 = Dtest_equality(C_comp, C, dim1 * dim2);
			if (succcess2) fprintf(stdout, " | Metadata-reuse run: OK");
			else fprintf(stdout, "| Metadata-reuse run: RAN, WRONG C");
			if(succcess && succcess2) fprintf(stdout, " | -> Test PASSED\n");
			else fprintf(stdout, " | -> Test FAILED\n");
			CHLMemcpy(C, C_buf, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
			CHLMemcpy(C_comp, C_buf, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
			CHLSyncCheckErr();
		}

		CHLVecInit(C, M * N, 44, C_loc);
		CHLMemcpy(C_comp, C,  M * N *sizeof(double), C_loc, C_loc);

		fprintf(stdout, "\n3) Testing Weird Problem dimensions < 100 MB:\n\n");
		alpha = 1.23;
		beta = 0.9876;
		for (int dim1 = 289; dim1 <= M; dim1*=4) for (int dim2 = 353; dim2 <= N; dim2*=4) for (int dim3 = 307; dim3 <= K; dim3*=4) if ( dim1 != dim2 || dim2 != dim3 || dim1!= dim3){
			fprintf(stdout, "[M, N, K] = [%d, %d, %d]\t-> ", dim1, dim2, dim3);
			cpu_timer = csecond();
			ret_autotune_val = PARALiADgemm(TransA, TransB, dim1, dim2, dim3, alpha, A, ldA, B, ldB, beta, C , ldC);
			cpu_timer  = csecond() - cpu_timer;
			double comp_flops =  Gval_per_s(gemm_ops(dim1,dim2,dim3),cpu_timer);
			fprintf(stderr, "M=%d,N=%d,K=%d: Gflops/s -> ", dim1, dim2, dim3);
			fprintf(stderr, "PARALiA: %.1lf, ", comp_flops);
			cpu_timer = csecond();
			T = fmin(dim1,fmin(dim2,dim3))/2;
			cuBLASXtDgemmWrap(TransA, TransB, dim1, dim2, dim3, alpha, A, ldA, B, ldB, beta, C_comp, ldC,  T, cache_limit, CHL_WORKERS, dev_ids);
			CHLSyncCheckErr();
			cpu_timer  = csecond() - cpu_timer;
			fprintf(stderr, "cuBLASXT: %.1lf\n", Gval_per_s(gemm_ops(dim1,dim2,dim3),cpu_timer));
			fprintf(stderr, "%s\n", ret_autotune_val->print_csv());
			if (comp_flops < Gval_per_s(gemm_ops(dim1,dim2,dim3),cpu_timer)) warning("Inferior Perf to cublasXt\n");
			int succcess = Dtest_equality(C_comp, C, dim1 * dim2); 
			if (succcess) fprintf(stdout, " | First run: OK");
			else fprintf(stdout, " | First run: RAN, WRONG C");
			CHLMemcpy(C, C_buf, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
			ret_autotune_val = PARALiADgemm(TransA, TransB, dim1, dim2, dim3, alpha, A, ldA, B, ldB, beta, C , ldC);
			int succcess2 = Dtest_equality(C_comp, C, dim1 * dim2);
			if (succcess2) fprintf(stdout, " | Metadata-reuse run: OK");
			else fprintf(stdout, "| Metadata-reuse run: RAN, WRONG C");
			if(succcess && succcess2) fprintf(stdout, " | -> Test PASSED\n");
			else fprintf(stdout, " | -> Test FAILED\n");
			CHLMemcpy(C, C_buf, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
			CHLMemcpy(C_comp, C_buf, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
			CHLSyncCheckErr();
		}

		CHLVecInit(C, M * N, 44, C_loc);
		CHLMemcpy(C_comp, C,  M * N *sizeof(double), C_loc, C_loc);

		fprintf(stdout, "\n4) Testing (Weird) Transpose Problems < 100 MB:\n\n");
		TransA = TransB = 'T';
		alpha = 1.23;
		beta = 0.9876;
		for  (int dim1 = 289; dim1 <= M; dim1*=4) for (int dim2 = 353; dim2 <= N; dim2*=4) for (int dim3 = 307; dim3 <= K; dim3*=4){
			fprintf(stdout, "[M, N, K] = [%d, %d, %d]\t-> ", dim1, dim2, dim3);
			cpu_timer = csecond();
			ret_autotune_val = PARALiADgemm(TransA, TransB, dim1, dim2, dim3, alpha, A, ldA, B, ldB, beta, C , ldC);
			cpu_timer  = csecond() - cpu_timer;
			double comp_flops =  Gval_per_s(gemm_ops(dim1,dim2,dim3),cpu_timer);
			fprintf(stderr, "M=%d,N=%d,K=%d: Gflops/s -> ", dim1, dim2, dim3);
			fprintf(stderr, "PARALiA: %.1lf, ", comp_flops);
			cpu_timer = csecond();
			T = fmin(dim1,fmin(dim2,dim3))/2;
			cuBLASXtDgemmWrap(TransA, TransB, dim1, dim2, dim3, alpha, A, ldA, B, ldB, beta, C_comp, ldC,  T, cache_limit, CHL_WORKERS, dev_ids);
			CHLSyncCheckErr();
			cpu_timer  = csecond() - cpu_timer;
			fprintf(stderr, "cuBLASXT: %.1lf\n", Gval_per_s(gemm_ops(dim1,dim2,dim3),cpu_timer));
			fprintf(stderr, "%s\n", ret_autotune_val->print_csv());
			if (comp_flops < Gval_per_s(gemm_ops(dim1,dim2,dim3),cpu_timer)) warning("Inferior Perf to cublasXt\n");
			int succcess = Dtest_equality(C_comp, C, dim1 * dim2); 
			if (succcess) fprintf(stdout, " | First run: OK");
			else fprintf(stdout, " | First run: RAN, WRONG C");
			CHLMemcpy(C, C_buf, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
			ret_autotune_val = PARALiADgemm(TransA, TransB, dim1, dim2, dim3, alpha, A, ldA, B, ldB, beta, C , ldC);
			int succcess2 = Dtest_equality(C_comp, C, dim1 * dim2);
			if (succcess2) fprintf(stdout, " | Metadata-reuse run: OK");
			else fprintf(stdout, "| Metadata-reuse run: RAN, WRONG C");
			if(succcess && succcess2) fprintf(stdout, " | -> Test PASSED\n");
			else fprintf(stdout, " | -> Test FAILED\n");
			CHLMemcpy(C, C_buf, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
			CHLMemcpy(C_comp, C_buf, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
			CHLSyncCheckErr();
		}
		CHLFree(A, M* K* sizeof(double), A_loc);
		CHLFree(B, N* K* sizeof(double), B_loc);
		CHLFree(C, N* M* sizeof(double), C_loc);
		CHLFree(C_comp, N* M* sizeof(double), CHL_MEMLOCS -1);
		CHLFree(C_buf, N* M* sizeof(double), CHL_MEMLOCS -1);
		CHLSyncCheckErr();
	}
	if(run_gpu_mem) {
		int ctr = 0; 
		fprintf(stdout, "\n-----------------------------------------Testing for matrices on GPU-----------------------------------------\n");
		for (int i = 0; i< CHL_WORKERS; i++){
			int dev_id = dev_ids[i];
			fprintf(stderr, "dgemm_tester: Allocating GPU buffers...->100 MB...");
			cpu_timer = csecond();
			A_loc = B_loc = C_loc = dev_id;
			A = (double*) CHLMalloc(M * K*sizeof(double), A_loc, 1);
			B = (double*) CHLMalloc(N * K*sizeof(double), B_loc, 1);
			C = (double*) CHLMalloc(M * N*sizeof(double), C_loc, 1);
			C_comp = (double*) CHLMalloc(M * N*sizeof(double), C_loc, 1);
			C_buf = (double*) CHLMalloc(M * N*sizeof(double), C_loc, 1);

			double* C_host_buf, * C_host_comp_buf;
			C_host_buf =  (double*) CHLMalloc(M * N*sizeof(double), CHL_MEMLOCS - 1, 1);
			C_host_comp_buf =  (double*) CHLMalloc(M * N*sizeof(double), CHL_MEMLOCS - 1, 1);
			CHLSyncCheckErr();
			cpu_timer  = csecond() - cpu_timer;
			fprintf(stderr, "done.\nAlloc time:\t%lf ms\n\n",  cpu_timer  * 1000);

			cpu_timer = csecond();
			fprintf(stderr, "Initializing to random values...");
			CHLVecInit(A, K * M, 42, A_loc);
			CHLVecInit(B, K * N, 43, B_loc);
			CHLVecInit(C, M * N, 44, C_loc);
			CHLMemcpy(C_host_comp_buf, C,  M * N *sizeof(double), CHL_MEMLOCS - 1, C_loc);
			CHLMemcpy(C_comp, C_host_comp_buf,  M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
			CHLMemcpy(C_buf, C_host_comp_buf,  M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);

			CHLSyncCheckErr();
			cpu_timer  = csecond() - cpu_timer ;
			fprintf(stderr, "done.\nInit time:\t%lf ms\n\n",  cpu_timer  * 1000);

			fprintf(stdout, "\n%d) Testing (weird) matrices in GPU(%d) mem < 100 MB:\n\n", ++ctr, dev_id);
			TransA = TransB = 'N';
			alpha = 1.23;
			beta = 0.9876;
			for (int dim1 = 289; dim1 <= M; dim1*=4) for (int dim2 = 353; dim2 <= N; dim2*=4) for (int dim3 = 307; dim3 <= K; dim3*=4){
				fprintf(stdout, "[M, N, K] = [%d, %d, %d]\t-> ", dim1, dim2, dim3);
				cpu_timer = csecond();
				ret_autotune_val = PARALiADgemm(TransA, TransB, dim1, dim2, dim3, alpha, A, ldA, B, ldB, beta, C , ldC);
				cpu_timer  = csecond() - cpu_timer;
				double comp_flops =  Gval_per_s(gemm_ops(dim1,dim2,dim3),cpu_timer);
				fprintf(stderr, "M=%d,N=%d,K=%d: Gflops/s -> ", dim1, dim2, dim3);
				fprintf(stderr, "PARALiA: %.1lf, ", comp_flops);
				cpu_timer = csecond();
				T = fmin(dim1,fmin(dim2,dim3))/2;
				cuBLASXtDgemmWrap(TransA, TransB, dim1, dim2, dim3, alpha, A, ldA, B, ldB, beta, C_comp, ldC,  T, cache_limit, CHL_WORKERS, dev_ids);
				CHLSyncCheckErr();
				cpu_timer  = csecond() - cpu_timer;
				fprintf(stderr, "cuBLASXT: %.1lf\n", Gval_per_s(gemm_ops(dim1,dim2,dim3),cpu_timer));
				fprintf(stderr, "%s\n", ret_autotune_val->print_csv());
				if (comp_flops < Gval_per_s(gemm_ops(dim1,dim2,dim3),cpu_timer)) warning("Inferior Perf to cublasXt\n");
				CHLMemcpy(C_host_buf, C,  dim1 * dim2 *sizeof(double), CHL_MEMLOCS - 1, C_loc);
				CHLMemcpy(C_host_comp_buf, C_comp,  dim1 * dim2 *sizeof(double), CHL_MEMLOCS - 1, C_loc);
				int succcess = Dtest_equality(C_host_comp_buf, C_host_buf, dim1 * dim2); 
				if (succcess) fprintf(stdout, " | First run: OK");
				else fprintf(stdout, " | First run: RAN, WRONG C");
				CHLMemcpy(C, C_buf, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
				ret_autotune_val = PARALiADgemm(TransA, TransB, dim1, dim2, dim3, alpha, A, ldA, B, ldB, beta, C , ldC);
				CHLMemcpy(C_host_buf, C,  dim1 * dim2 *sizeof(double), CHL_MEMLOCS - 1, C_loc);
				int succcess2 = Dtest_equality(C_host_comp_buf, C_host_buf, dim1 * dim2);
				if (succcess2) fprintf(stdout, " | Metadata-reuse run: OK");
				else fprintf(stdout, "| Metadata-reuse run: RAN, WRONG C");
				if(succcess && succcess2) fprintf(stdout, " | -> Test PASSED\n");
				else fprintf(stdout, " | -> Test FAILED\n");
				CHLMemcpy(C, C_buf, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
				CHLMemcpy(C_comp, C_buf, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
				CHLSyncCheckErr();
			}

			fprintf(stdout, "\n%d) Testing (weird) matrices in GPU(%d) mem + Transpose < 100 MB:\n\n", ++ctr, dev_id);
			TransA = TransB = 'T';
			alpha = 1.23;
			beta = 0.9876;
			for (int dim1 = 289; dim1 <= M; dim1*=4) for (int dim2 = 353; dim2 <= N; dim2*=4) for (int dim3 = 307; dim3 <= K; dim3*=4){
				fprintf(stdout, "[M, N, K] = [%d, %d, %d]\t-> ", dim1, dim2, dim3);
				cpu_timer = csecond();
				ret_autotune_val = PARALiADgemm(TransA, TransB, dim1, dim2, dim3, alpha, A, ldA, B, ldB, beta, C , ldC);
				cpu_timer  = csecond() - cpu_timer;
				double comp_flops =  Gval_per_s(gemm_ops(dim1,dim2,dim3),cpu_timer);
				fprintf(stderr, "M=%d,N=%d,K=%d: Gflops/s -> ", dim1, dim2, dim3);
				fprintf(stderr, "PARALiA: %.1lf, ", comp_flops);
				cpu_timer = csecond();
				T = fmin(dim1,fmin(dim2,dim3))/2;
				cuBLASXtDgemmWrap(TransA, TransB, dim1, dim2, dim3, alpha, A, ldA, B, ldB, beta, C_comp, ldC,  T, cache_limit, CHL_WORKERS, dev_ids);
				CHLSyncCheckErr();
				cpu_timer  = csecond() - cpu_timer;
				fprintf(stderr, "cuBLASXT: %.1lf\n", Gval_per_s(gemm_ops(dim1,dim2,dim3),cpu_timer));
				fprintf(stderr, "%s\n", ret_autotune_val->print_csv());
				if (comp_flops < Gval_per_s(gemm_ops(dim1,dim2,dim3),cpu_timer)) warning("Inferior Perf to cublasXt\n");
				CHLMemcpy(C_host_buf, C,  dim1 * dim2 *sizeof(double), CHL_MEMLOCS - 1, C_loc);
				CHLMemcpy(C_host_comp_buf, C_comp,  dim1 * dim2 *sizeof(double), CHL_MEMLOCS - 1, C_loc);
				int succcess = Dtest_equality(C_host_comp_buf, C_host_buf, dim1 * dim2); 
				if (succcess) fprintf(stdout, " | First run: OK");
				else fprintf(stdout, " | First run: RAN, WRONG C");
				CHLMemcpy(C, C_buf, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
				ret_autotune_val = PARALiADgemm(TransA, TransB, dim1, dim2, dim3, alpha, A, ldA, B, ldB, beta, C , ldC);
				CHLMemcpy(C_host_buf, C,  dim1 * dim2 *sizeof(double), CHL_MEMLOCS - 1, C_loc);
				int succcess2 = Dtest_equality(C_host_comp_buf, C_host_buf, dim1 * dim2);
				if (succcess2) fprintf(stdout, " | Metadata-reuse run: OK");
				else fprintf(stdout, "| Metadata-reuse run: RAN, WRONG C");
				if(succcess && succcess2) fprintf(stdout, " | -> Test PASSED\n");
				else fprintf(stdout, " | -> Test FAILED\n");
				CHLMemcpy(C, C_buf, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
				CHLMemcpy(C_comp, C_buf, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
				CHLSyncCheckErr();
			}

			CHLFree(A, M* K* sizeof(double), A_loc);
			CHLFree(B, N* K* sizeof(double), B_loc);
			CHLFree(C, N* M* sizeof(double), C_loc);
			CHLFree(C_comp, N* M* sizeof(double), C_loc);
			CHLFree(C_buf, N* M* sizeof(double), C_loc);

			CHLFree(C_host_buf, N* M* sizeof(double), CHL_MEMLOCS - 1);
			CHLFree(C_host_comp_buf, N* M* sizeof(double), CHL_MEMLOCS - 1);
		}
		A_loc = 0;
		if (CHL_WORKERS == 1) B_loc = C_loc = 0;
		else if (CHL_WORKERS == 2){
			B_loc = 1;
			C_loc = 1;
		}
		else if (CHL_WORKERS > 2){
			B_loc = 1;
			C_loc = 2;
		}

		fprintf(stderr, "dgemm_tester: Allocating Mixed GPU buffers...-> A(dev=%d) : %.3lf GB, B(dev=%d) : %.3lf GB, C(dev=%d) : %.3lf GB(x2 for check):", A_loc, M*K*sizeof(double)/1e9, B_loc, K*N*sizeof(double)/1e9, C_loc, M*N*sizeof(double)/1e9);
		cpu_timer = csecond();

		A = (double*) CHLMalloc(M * K*sizeof(double), A_loc, 1);
		B = (double*) CHLMalloc(N * K*sizeof(double), B_loc, 1);
		C = (double*) CHLMalloc(M * N*sizeof(double), C_loc, 1);
		C_comp = (double*) CHLMalloc(M * N*sizeof(double), C_loc, 1);
		C_buf = (double*) CHLMalloc(M * N*sizeof(double), C_loc, 1);

		double* C_host_buf, * C_host_comp_buf;
		C_host_buf =  (double*) CHLMalloc(M * N*sizeof(double), CHL_MEMLOCS - 1, 1);
		C_host_comp_buf =  (double*) CHLMalloc(M * N*sizeof(double), CHL_MEMLOCS - 1, 1);
		CHLSyncCheckErr();
		cpu_timer  = csecond() - cpu_timer;
		fprintf(stderr, "done.\nAlloc time:\t%lf ms\n\n",  cpu_timer  * 1000);

		cpu_timer = csecond();
		fprintf(stderr, "Initializing to random values...");
		CHLVecInit(A, K * M, 42, A_loc);
		CHLVecInit(B, K * N, 43, B_loc);
		CHLVecInit(C, M * N, 44, C_loc);
		CHLMemcpy(C_host_comp_buf, C,  M * N *sizeof(double), CHL_MEMLOCS - 1, C_loc);
		CHLMemcpy(C_comp, C_host_comp_buf,  M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
		CHLMemcpy(C_buf, C_host_comp_buf,  M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
		CHLSyncCheckErr();
		cpu_timer  = csecond() - cpu_timer ;
		fprintf(stderr, "done.\nInit time:\t%lf ms\n\n",  cpu_timer  * 1000);

		fprintf(stdout, "\n%d) Testing mixed GPU [A_loc, B_loc, C_loc] = [%d, %d, %d] matrices < 100 MB:\n\n", ++ctr, A_loc, B_loc, C_loc);
		TransA = TransB = 'N';
		alpha = 1.23;
		beta = 0.9876;
		for (int dim1 = 289; dim1 <= M; dim1*=4) for (int dim2 = 353; dim2 <= N; dim2*=4) for (int dim3 = 307; dim3 <= K; dim3*=4){
			fprintf(stdout, "[M, N, K] = [%d, %d, %d]\t-> ", dim1, dim2, dim3);
			cpu_timer = csecond();
			ret_autotune_val = PARALiADgemm(TransA, TransB, dim1, dim2, dim3, alpha, A, ldA, B, ldB, beta, C , ldC);
			cpu_timer  = csecond() - cpu_timer;
			double comp_flops =  Gval_per_s(gemm_ops(dim1,dim2,dim3),cpu_timer);
			fprintf(stderr, "M=%d,N=%d,K=%d: Gflops/s -> ", dim1, dim2, dim3);
			fprintf(stderr, "PARALiA: %.1lf, ", comp_flops);
			cpu_timer = csecond();
			T = fmin(dim1,fmin(dim2,dim3))/2;
			cuBLASXtDgemmWrap(TransA, TransB, dim1, dim2, dim3, alpha, A, ldA, B, ldB, beta, C_comp, ldC,  T, cache_limit, CHL_WORKERS, dev_ids);
			CHLSyncCheckErr();
			cpu_timer  = csecond() - cpu_timer;
			fprintf(stderr, "cuBLASXT: %.1lf\n", Gval_per_s(gemm_ops(dim1,dim2,dim3),cpu_timer));
			fprintf(stderr, "%s\n", ret_autotune_val->print_csv());
			if (comp_flops < Gval_per_s(gemm_ops(dim1,dim2,dim3),cpu_timer)) warning("Inferior Perf to cublasXt\n");
			CHLMemcpy(C_host_buf, C,  dim1 * dim2 *sizeof(double), CHL_MEMLOCS - 1, C_loc);
			CHLMemcpy(C_host_comp_buf, C_comp,  dim1 * dim2 *sizeof(double), CHL_MEMLOCS - 1, C_loc);
			int succcess = Dtest_equality(C_host_comp_buf, C_host_buf, dim1 * dim2); 
			if (succcess) fprintf(stdout, " | First run: OK");
			else fprintf(stdout, " | First run: RAN, WRONG C");
			CHLMemcpy(C, C_buf, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
			ret_autotune_val = PARALiADgemm(TransA, TransB, dim1, dim2, dim3, alpha, A, ldA, B, ldB, beta, C , ldC);
			CHLMemcpy(C_host_buf, C,  dim1 * dim2 *sizeof(double), CHL_MEMLOCS - 1, C_loc);
			int succcess2 = Dtest_equality(C_host_comp_buf, C_host_buf, dim1 * dim2);
			if (succcess2) fprintf(stdout, " | Metadata-reuse run: OK");
			else fprintf(stdout, "| Metadata-reuse run: RAN, WRONG C");
			if(succcess && succcess2) fprintf(stdout, " | -> Test PASSED\n");
			else fprintf(stdout, " | -> Test FAILED\n");
			CHLMemcpy(C, C_buf, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
			CHLMemcpy(C_comp, C_buf, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
			CHLSyncCheckErr();
		}
		CHLFree(A, M* K* sizeof(double), A_loc);
		CHLFree(B, N* K* sizeof(double), B_loc);
		CHLFree(C, N* M* sizeof(double), C_loc);
		CHLFree(C_comp, N* M* sizeof(double), C_loc);
		CHLFree(C_buf, N* M* sizeof(double), C_loc);

		CHLFree(C_host_buf, N* M* sizeof(double), CHL_MEMLOCS - 1);
		CHLFree(C_host_comp_buf, N* M* sizeof(double), CHL_MEMLOCS - 1);
	}
	if (run_cpu_mem && run_gpu_mem){
		A_loc = C_loc = CHL_MEMLOCS - 1;
		B_loc = 0;

		fprintf(stderr, "dgemm_tester: Allocating Mixed CPU/GPU buffers...-> A(dev=%d) : %.3lf GB, B(dev=%d) : %.3lf GB, C(dev=%d) : %.3lf GB(x2 for check):", A_loc, M*K*sizeof(double)/1e9, B_loc, K*N*sizeof(double)/1e9, C_loc, M*N*sizeof(double)/1e9);
		cpu_timer = csecond();

		A = (double*) CHLMalloc(M * K*sizeof(double), A_loc, 1);
		B = (double*) CHLMalloc(N * K*sizeof(double), B_loc, 1);
		C = (double*) CHLMalloc(M * N*sizeof(double), C_loc, 1);
		C_comp = (double*) CHLMalloc(M * N*sizeof(double), C_loc, 1);
		C_buf = (double*) CHLMalloc(M * N*sizeof(double), C_loc, 1);

		double* C_host_buf, * C_host_comp_buf;
		C_host_buf =  (double*) CHLMalloc(M * N*sizeof(double), CHL_MEMLOCS - 1, 1);
		C_host_comp_buf =  (double*) CHLMalloc(M * N*sizeof(double), CHL_MEMLOCS - 1, 1);
		CHLSyncCheckErr();
		cpu_timer  = csecond() - cpu_timer;
		fprintf(stderr, "done.\nAlloc time:\t%lf ms\n\n",  cpu_timer  * 1000);

		cpu_timer = csecond();
		fprintf(stderr, "Initializing to random values...");
		CHLVecInit(A, K * M, 42, A_loc);
		CHLVecInit(B, K * N, 43, B_loc);
		CHLVecInit(C, M * N, 44, C_loc);
		CHLMemcpy(C_host_comp_buf, C,  M * N *sizeof(double), CHL_MEMLOCS - 1, C_loc);
		CHLMemcpy(C_comp, C_host_comp_buf,  M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
		CHLMemcpy(C_buf, C_host_comp_buf,  M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
		CHLSyncCheckErr();
		cpu_timer  = csecond() - cpu_timer ;
		fprintf(stderr, "done.\nInit time:\t%lf ms\n\n",  cpu_timer  * 1000);

		fprintf(stdout, "\n1) Testing mixed CPU/GPU [A_loc, B_loc, C_loc] = [%d, %d, %d] matrices < 100 MB:\n\n", A_loc, B_loc, C_loc);
		TransA = TransB = 'N';
		alpha = 1.23;
		beta = 0.9876;
		for (int dim1 = 289; dim1 <= M; dim1*=4) for (int dim2 = 353; dim2 <= N; dim2*=4) for (int dim3 = 307; dim3 <= K; dim3*=4){
			fprintf(stdout, "[M, N, K] = [%d, %d, %d]\t-> ", dim1, dim2, dim3);
			cpu_timer = csecond();
			ret_autotune_val = PARALiADgemm(TransA, TransB, dim1, dim2, dim3, alpha, A, ldA, B, ldB, beta, C , ldC);
			cpu_timer  = csecond() - cpu_timer;
			double comp_flops =  Gval_per_s(gemm_ops(dim1,dim2,dim3),cpu_timer);
			fprintf(stderr, "M=%d,N=%d,K=%d: Gflops/s -> ", dim1, dim2, dim3);
			fprintf(stderr, "PARALiA: %.1lf, ", comp_flops);
			cpu_timer = csecond();
			T = fmin(dim1,fmin(dim2,dim3))/2;
			cuBLASXtDgemmWrap(TransA, TransB, dim1, dim2, dim3, alpha, A, ldA, B, ldB, beta, C_comp, ldC,  T, cache_limit, CHL_WORKERS, dev_ids);
			CHLSyncCheckErr();
			cpu_timer  = csecond() - cpu_timer;
			fprintf(stderr, "cuBLASXT: %.1lf\n", Gval_per_s(gemm_ops(dim1,dim2,dim3),cpu_timer));
			fprintf(stderr, "%s\n", ret_autotune_val->print_csv());
			if (comp_flops < Gval_per_s(gemm_ops(dim1,dim2,dim3),cpu_timer)) warning("Inferior Perf to cublasXt\n");
			CHLMemcpy(C_host_buf, C,  dim1 * dim2 *sizeof(double), CHL_MEMLOCS - 1, C_loc);
			CHLMemcpy(C_host_comp_buf, C_comp,  dim1 * dim2 *sizeof(double), CHL_MEMLOCS - 1, C_loc);
			int succcess = Dtest_equality(C_host_comp_buf, C_host_buf, dim1 * dim2); 
			if (succcess) fprintf(stdout, " | First run: OK");
			else fprintf(stdout, " | First run: RAN, WRONG C");
			CHLMemcpy(C, C_buf, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
			ret_autotune_val = PARALiADgemm(TransA, TransB, dim1, dim2, dim3, alpha, A, ldA, B, ldB, beta, C , ldC);
			CHLMemcpy(C_host_buf, C,  dim1 * dim2 *sizeof(double), CHL_MEMLOCS - 1, C_loc);
			int succcess2 = Dtest_equality(C_host_comp_buf, C_host_buf, dim1 * dim2);
			if (succcess2) fprintf(stdout, " | Metadata-reuse run: OK");
			else fprintf(stdout, "| Metadata-reuse run: RAN, WRONG C");
			if(succcess && succcess2) fprintf(stdout, " | -> Test PASSED\n");
			else fprintf(stdout, " | -> Test FAILED\n");
			CHLMemcpy(C, C_buf, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
			CHLMemcpy(C_comp, C_buf, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
			CHLSyncCheckErr();
		}
		CHLFree(A, M* K* sizeof(double), A_loc);
		CHLFree(B, N* K* sizeof(double), B_loc);
		CHLFree(C, N* M* sizeof(double), C_loc);
		CHLFree(C_comp, N* M* sizeof(double), C_loc);
		CHLFree(C_buf, N* M* sizeof(double), C_loc);

		CHLFree(C_host_buf, N* M* sizeof(double), CHL_MEMLOCS - 1);
		CHLFree(C_host_comp_buf, N* M* sizeof(double), CHL_MEMLOCS - 1);
	}

	if (run_large) error("dgemm_tester: PARALiA 3.0 (and/or dgemm_tester) not updated for large mem\n");
	if (run_cpu_mem && run_large){
		ldA = ldB = ldC = M = N = K = (long int) 1.5*CHLGetMaxDimSqAsset2D(3, sizeof(double), 256, 0);
		fprintf(stdout, "\n----------------------------------------------------------------------------------\n");
		fprintf(stdout, "dgemm_tester: Allocating CPU buffers...-> %.3lf GB:", (gemm_mem_ops(M,N,K) + M * N)* sizeof(double)/1e9);
		cpu_timer = csecond();

		A_loc = B_loc = C_loc = CHL_MEMLOCS - 1;

		double *A, *B, *C;
		A = (double*) CHLMalloc(M * K*sizeof(double), A_loc, 1);
		B = (double*) CHLMalloc(N * K*sizeof(double), B_loc, 1);
		C = (double*) CHLMalloc(M * N*sizeof(double), C_loc, 1);

		CHLSyncCheckErr();
		cpu_timer  = csecond() - cpu_timer;
		fprintf(stderr, "done.\nAlloc time:\t%lf ms\n\n",  cpu_timer  * 1000);
		cpu_timer = csecond();
		fprintf(stderr, "Initializing to random values...");
		CHLVecInit(A, K * M, 42, A_loc);
		CHLVecInit(B, K * N, 43, B_loc);
		CHLVecInit(C, M * N, 44, C_loc);
		CHLSyncCheckErr();
		cpu_timer  = csecond() - cpu_timer ;
		fprintf(stderr, "done.\nInit time:\t%lf ms\n\n",  cpu_timer  * 1000);

		double *C_comp = (double*) CHLMalloc(M * N*sizeof(double), CHL_MEMLOCS - 1, 1);;
		CHLMemcpy(C_comp, C,  M * N *sizeof(double), CHL_MEMLOCS - 1, C_loc);

		fprintf(stdout, "\n----------------------------------------------------------------------------------\n");
		fprintf(stdout, "dgemm_tester: Testing Square Problem: %.3lf GB:", gemm_mem_ops(M,N,K) * sizeof(double)/1e9);
		TransA = TransB = 'N';
		alpha = 1.23;
		beta = 0.9876;
		cpu_timer = csecond();
		ret_autotune_val = PARALiADgemm(TransA, TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C , ldC);
		cpu_timer  = csecond() - cpu_timer;
		for (int i = 0; i< CHL_MEMLOCS; i++) PARALiADevCacheFree(i);
		CHLSyncCheckErr();
		double comp_flops = Gval_per_s(gemm_ops(M,N,K),cpu_timer);
		fprintf(stderr, "M=%zu, N=%zu, K=%zu: Gflops/s -> ", M, N, K);
		fprintf(stderr, "PARALiA: %.1lf, ", comp_flops);
		cpu_timer = csecond();
		T = fmin(M,fmin(N,K))/4;
		cuBLASXtDgemmWrap(TransA, TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C_comp, ldC,  T, cache_limit, CHL_WORKERS, dev_ids);
		CHLSyncCheckErr();
		cpu_timer  = csecond() - cpu_timer;
		fprintf(stderr, "cuBLASXT: %.1lf\n", Gval_per_s(gemm_ops(M,N,K),cpu_timer));
		fprintf(stderr, "%s\n", ret_autotune_val->print_csv());
		if (comp_flops < Gval_per_s(gemm_ops(M,N,K),cpu_timer)) warning("Inferior Perf to cublasXt\n");
		Dtest_equality(C_comp, C, M * N);
		CHLMemcpy(C, C_comp, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
		CHLSyncCheckErr();

		CHLVecInit(C, M * N, 44, C_loc);
		CHLMemcpy(C_comp, C,  M * N *sizeof(double), CHL_MEMLOCS - 1, C_loc);

		M = (long int) M/1.24223;
		N = (long int) N/1.34645;
		K = (long int) K/2.18321;
		fprintf(stdout, "\n----------------------------------------------------------------------------------\n");
		fprintf(stdout, "dgemm_tester: Testing Weird Non-Square Problem: %.3lf GB:", gemm_mem_ops(M,N,K) * sizeof(double)/1e9);
		alpha = 1.23;
		beta = 0.9876;
		cpu_timer = csecond();
		ret_autotune_val = PARALiADgemm(TransA, TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C , ldC);
		cpu_timer  = csecond() - cpu_timer;
		for (int i = 0; i< CHL_MEMLOCS; i++) PARALiADevCacheFree(i);
		CHLSyncCheckErr();
		comp_flops = Gval_per_s(gemm_ops(M,N,K),cpu_timer);
		fprintf(stderr, "M=%zu, N=%zu, K=%zu: Gflops/s -> ", M, N, K);
		fprintf(stderr, "PARALiA: %.1lf, ", comp_flops);
		cpu_timer = csecond();
		T = fmin(M,fmin(N,K))/4;
		cuBLASXtDgemmWrap(TransA, TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C_comp, ldC,  T, cache_limit, CHL_WORKERS, dev_ids);
		CHLSyncCheckErr();
		cpu_timer  = csecond() - cpu_timer;
		fprintf(stderr, "cuBLASXT: %.1lf\n", Gval_per_s(gemm_ops(M,N,K),cpu_timer));
		fprintf(stderr, "%s\n", ret_autotune_val->print_csv());
		if (comp_flops < Gval_per_s(gemm_ops(M,N,K),cpu_timer)) warning("Inferior Perf to cublasXt\n");
		Dtest_equality(C_comp, C, M * N);
		CHLMemcpy(C, C_comp, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
		CHLSyncCheckErr();

		CHLVecInit(C, M * N, 44, C_loc);
		CHLMemcpy(C_comp, C,  M * N *sizeof(double), CHL_MEMLOCS - 1, C_loc);

		fprintf(stdout, "\n----------------------------------------------------------------------------------\n");
		fprintf(stdout, "dgemm_tester: Testing Large Transpose\n\n");
		TransA = TransB = 'T';
		alpha = 1.23;
		beta = 0.9876;
		cpu_timer = csecond();
		ret_autotune_val = PARALiADgemm(TransA, TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C , ldC);
		cpu_timer  = csecond() - cpu_timer;
		for (int i = 0; i< CHL_MEMLOCS; i++) PARALiADevCacheFree(i);
		CHLSyncCheckErr();
		comp_flops = Gval_per_s(gemm_ops(M,N,K),cpu_timer);
		fprintf(stderr, "M=%zu, N=%zu, K=%zu: Gflops/s -> ", M, N, K);
		fprintf(stderr, "PARALiA: %.1lf, ", comp_flops);
		cpu_timer = csecond();
		T = fmin(M,fmin(N,K))/4;
		cuBLASXtDgemmWrap(TransA, TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C_comp, ldC,  T, cache_limit, CHL_WORKERS, dev_ids);
		CHLSyncCheckErr();
		cpu_timer  = csecond() - cpu_timer;
		fprintf(stderr, "cuBLASXT: %.1lf\n", Gval_per_s(gemm_ops(M,N,K),cpu_timer));
		fprintf(stderr, "%s\n", ret_autotune_val->print_csv());
		if (comp_flops < Gval_per_s(gemm_ops(M,N,K),cpu_timer)) warning("Inferior Perf to cublasXt\n");
		Dtest_equality(C_comp, C, M * N);
		CHLMemcpy(C, C_comp, M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
		CHLSyncCheckErr();

		CHLFree(A, M* K* sizeof(double), A_loc);
		CHLFree(B, N* K* sizeof(double), A_loc);
		CHLFree(C, N* M* sizeof(double), A_loc);
		CHLFree(C_comp, N* M* sizeof(double), CHL_MEMLOCS - 1);
		CHLSyncCheckErr();
	}
	if (run_gpu_mem && run_large){
		A_loc = 0;
		if (CHL_WORKERS == 1){
			B_loc = C_loc = 0;
			ldA = ldB = ldC = M = N = K = (long int) CHLGetMaxDimSqAsset2D(4, sizeof(double), 256, 0);
		}
		else if (CHL_WORKERS == 2){
			B_loc = 0;
			C_loc = 1;
			ldA = ldB = ldC = M = N = K = (long int) CHLGetMaxDimSqAsset2D(2, sizeof(double), 256, 0);
		}
		else if (CHL_WORKERS > 2){
			B_loc = 1;
			C_loc = 2;
			ldA = ldB = ldC = M = N = K = (long int) CHLGetMaxDimSqAsset2D(2, sizeof(double), 256, 0);
		}

		fprintf(stdout, "\n----------------------------------------------------------------------------------\n");
		fprintf(stdout, "dgemm_tester: Allocating Mixed GPU buffers...-> A(dev=%d) : %.3lf GB, B(dev=%d) : %.3lf GB, C(dev=%d) : %.3lf GB(x2 for check):", A_loc, M*K*sizeof(double)/1e9, B_loc, K*N*sizeof(double)/1e9, C_loc, M*N*sizeof(double)/1e9);
		cpu_timer = csecond();

		A = (double*) CHLMalloc(M * K*sizeof(double), A_loc, 1);
		B = (double*) CHLMalloc(N * K*sizeof(double), B_loc, 1);
		C = (double*) CHLMalloc(M * N*sizeof(double), C_loc, 1);
		C_comp = (double*) CHLMalloc(M * N*sizeof(double), C_loc, 1);

		double* C_host_buf, * C_host_comp_buf;
		C_host_buf =  (double*) CHLMalloc(M * N*sizeof(double), CHL_MEMLOCS - 1, 1);
		C_host_comp_buf =  (double*) CHLMalloc(M * N*sizeof(double), CHL_MEMLOCS - 1, 1);
		CHLSyncCheckErr();
		cpu_timer  = csecond() - cpu_timer;
		fprintf(stderr, "done.\nAlloc time:\t%lf ms\n\n",  cpu_timer  * 1000);

		cpu_timer = csecond();
		fprintf(stderr, "Initializing to random values...");
		CHLVecInit(A, K * M, 42, A_loc);
		CHLVecInit(B, K * N, 43, B_loc);
		CHLVecInit(C, M * N, 44, C_loc);
		CHLMemcpy(C_host_comp_buf, C,  M * N *sizeof(double), CHL_MEMLOCS - 1, C_loc);
		CHLMemcpy(C_comp, C_host_comp_buf,  M * N *sizeof(double), C_loc, CHL_MEMLOCS - 1);
		CHLSyncCheckErr();
		cpu_timer  = csecond() - cpu_timer ;
		fprintf(stderr, "done.\nInit time:\t%lf ms\n\n",  cpu_timer  * 1000);

		fprintf(stdout, "\n----------------------------------------------------------------------------------\n");
		fprintf(stdout, "dgemm_tester: Testing Large Matrices In GPU\n\n");
		TransA = TransB = 'N';
		alpha = 1.23;
		beta = 0.9876;
		cpu_timer = csecond();
		ret_autotune_val = PARALiADgemm(TransA, TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C , ldC);
		cpu_timer  = csecond() - cpu_timer;
		for (int i = 0; i< CHL_MEMLOCS; i++) PARALiADevCacheFree(i);
		CHLSyncCheckErr();
		double comp_flops =  Gval_per_s(gemm_ops(M,N,K),cpu_timer);
		fprintf(stderr, "M=%zu,N=%zu,K=%zu: Gflops/s -> ", M, N, K);
		fprintf(stderr, "PARALiA: %.1lf, ", comp_flops);
		cpu_timer = csecond();
		T = fmin(M,fmin(N,K))/4;
		cuBLASXtDgemmWrap(TransA, TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C_comp, ldC,  T, cache_limit, CHL_WORKERS, dev_ids);
		CHLSyncCheckErr();
		cpu_timer  = csecond() - cpu_timer;
		fprintf(stderr, "cuBLASXT: %.1lf\n", Gval_per_s(gemm_ops(M, N, K),cpu_timer));
		fprintf(stderr, "%s\n", ret_autotune_val->print_csv());
		if (comp_flops < Gval_per_s(gemm_ops(M, N, K),cpu_timer)) warning("Inferior Perf to cublasXt\n");
		CHLMemcpy(C_host_buf, C,  M * N *sizeof(double), CHL_MEMLOCS - 1, C_loc);
		CHLMemcpy(C_host_comp_buf, C_comp,  M * N  *sizeof(double), CHL_MEMLOCS - 1, C_loc);
		Dtest_equality(C_host_comp_buf, C_host_buf, M * N);
		CHLMemcpy(C, C_comp, M * N *sizeof(double), C_loc, C_loc);
		CHLSyncCheckErr();

		fprintf(stdout, "\n----------------------------------------------------------------------------------\n");
		fprintf(stdout, "dgemm_tester: Testing Large Matrices In GPUmem + Transpose\n\n");
		TransA = TransB = 'T';
		alpha = 1.23;
		beta = 0.9876;
		cpu_timer = csecond();
		ret_autotune_val = PARALiADgemm(TransA, TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C , ldC);
		cpu_timer  = csecond() - cpu_timer;
		for (int i = 0; i< CHL_MEMLOCS; i++) PARALiADevCacheFree(i);
		CHLSyncCheckErr();
		comp_flops =  Gval_per_s(gemm_ops(M,N,K),cpu_timer);
		fprintf(stderr, "M=%zu,N=%zu,K=%zu: Gflops/s -> ", M, N, K);
		fprintf(stderr, "PARALiA: %.1lf, ", comp_flops);
		cpu_timer = csecond();
		T = fmin(M,fmin(N,K))/4;
		cuBLASXtDgemmWrap(TransA, TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C_comp, ldC,  T, cache_limit, CHL_WORKERS, dev_ids);
		CHLSyncCheckErr();
		cpu_timer  = csecond() - cpu_timer;
		fprintf(stderr, "cuBLASXT: %.1lf\n", Gval_per_s(gemm_ops(M, N, K),cpu_timer));
		fprintf(stderr, "%s\n", ret_autotune_val->print_csv());
		if (comp_flops < Gval_per_s(gemm_ops(M, N, K),cpu_timer)) warning("Inferior Perf to cublasXt\n");
		CHLMemcpy(C_host_buf, C,  M * N *sizeof(double), CHL_MEMLOCS - 1, C_loc);
		CHLMemcpy(C_host_comp_buf, C_comp,  M * N  *sizeof(double), CHL_MEMLOCS - 1, C_loc);
		Dtest_equality(C_host_comp_buf, C_host_buf, M * N);
		CHLMemcpy(C, C_comp, M * N *sizeof(double), C_loc, C_loc);
		CHLSyncCheckErr();

		CHLFree(A, M* K* sizeof(double), A_loc);
		CHLFree(B, N* K* sizeof(double), B_loc);
		CHLFree(C, N* M* sizeof(double), C_loc);
		CHLFree(C_comp, N* M* sizeof(double), C_loc);
		CHLFree(C_host_buf, N* M* sizeof(double), CHL_MEMLOCS - 1);
		CHLFree(C_host_comp_buf, N* M* sizeof(double), CHL_MEMLOCS - 1);
	}
	fprintf(stdout, "\n-----------------------------------------GEMM (dtype=double): Tests Finished-----------------------------------------\n");
	return 0;
}
