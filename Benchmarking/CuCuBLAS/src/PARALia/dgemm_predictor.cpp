///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief The start of Zawarudo
///

#include "PARALiA.hpp"
#include "BackenedLibsWrapped.hpp"
#include "Testing.hpp"

#include "backend_wrappers.hpp"

#define CBLASXT_MAX_SAFE_TILE 10000

int main(const int argc, const char *argv[]) {
	char TransA, TransB;
  	double alpha, beta;
	long int M, N, K;
	int A_loc, B_loc, C_loc, C_out_loc;
	ATC_p predef_control_values = NULL, return_values = NULL;
	ParseInputLvl3(argc, argv, &predef_control_values, &TransA, &TransB, &alpha, &beta, &M, &N, &K, &A_loc, &B_loc, &C_loc, &C_out_loc);

	char *filename = (char *) malloc(1024* sizeof(char));
	if (predef_control_values!= NULL){
		if(predef_control_values->T > 0) {
			if (predef_control_values->T > M || predef_control_values->T > N || predef_control_values->T > K)
				error("Given Tin=%ld bigger than problem dim\n", predef_control_values->T);
			else if (predef_control_values->T > M/1.5 && predef_control_values->T > N/1.5 && predef_control_values->T > K/1.5)
				warning("Given Tin=%ld bigger than all problem dims/1.5\n", predef_control_values->T);
		}
		sprintf(filename, "%s/dgemm_predictor_predefined_vals_%s_%s.log",
			TESTLIBDIR, CoCoImplementationPrint(), VERSION);
	}
	else sprintf(filename, "%s/dgemm_predictor_%s_%s.log",
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

	double *A_ptr, *B_ptr, *C_ptr;
	// allocate in device if loc = 0, otherwise allocate in pinned memory for benchmarks
	A_ptr = (double*) CHLMalloc(M*K*sizeof(double), A_loc, 0);
	B_ptr = (double*) CHLMalloc(N*K*sizeof(double), B_loc, 0);
	C_ptr = (double*) CHLMalloc(M*N*sizeof(double), C_loc, 1);

	gemm_backend_in* initial_gemm = (gemm_backend_in*) malloc(sizeof(struct gemm_backend_in));

	initial_gemm->TransA = TransA;
	initial_gemm->TransB = TransB;
	initial_gemm->M = M;
	initial_gemm->N = N;
	initial_gemm->K = K;
	initial_gemm->A = (void**) &A_ptr;
	initial_gemm->B = (void**) &B_ptr;
	initial_gemm->C = (void**) &C_ptr;
	initial_gemm->alpha = &alpha;
	initial_gemm->beta = &beta;
	initial_gemm->ldA = ldA;
	initial_gemm->ldB = ldB;
	initial_gemm->ldC = ldC;
	initial_gemm->dev_id = -1;


	ATC_p autotune_controller_gemm = NULL;
	if (autotune_controller_gemm == NULL) autotune_controller_gemm = new ATC();
	if (predef_control_values && autotune_controller_gemm->diff_intialized_params_ATC(predef_control_values))
		 autotune_controller_gemm->mimic_ATC(predef_control_values);

	double autotune_timer = autotune_controller_gemm->autotune_problem("MM_FP64", A_loc, B_loc, C_loc, C_out_loc, M, N, K, sizeof(double));

	lprintf(0, "dgemm_predictor: t_autotune = %lf ms\n", autotune_timer*1000);
	autotune_controller_gemm->print();
	StoreLogLvl3(filename, autotune_controller_gemm, TransA, TransB, alpha, beta, M, N, K, A_loc, B_loc, C_loc, C_out_loc, cpu_timer, autotune_controller_gemm->pred_t, autotune_controller_gemm->pred_J);

	CHLSyncCheckErr();
	CHLFree(A_ptr, A_loc, M*K*sizeof(double));
	CHLFree(B_ptr, B_loc, N*K*sizeof(double));
	CHLFree(C_ptr, C_loc, M*N*sizeof(double));
	return 0;
}
