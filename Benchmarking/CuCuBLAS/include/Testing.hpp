///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief The external wrapper for CoCoPelia + wrapped cuBLASXt
///
#ifndef TESTING_H
#define TESTING_H

/// Return a string with the active Cmake implemetation flag used
char* CoCoImplementationPrint();

/// Return a string with the active Cmake Subkernel Distribution flag
char* CoCoDistributionPrint();

void ParseInputLvl3(const int argc, const char *argv[], ATC_p* predef_control_values, char* TransA, char* TransB,
	double* alpha, double* beta, long int* D1, long int* D2, long int* D3, short* loc1, short* loc2, short* loc3, short* outloc);
void CheckLogLvl3(char* filename, ATC_p predef_control_values, char TransA, char TransB,
	double alpha, double beta, long int D1, long int D2, long int D3, short loc1, short loc2, short loc3, short outloc);
void StoreLogLvl3(char* filename, ATC_p predef_control_values, char TransA, char TransB, double alpha, double beta,
	long int D1, long int D2, long int D3, short loc1, short loc2, short loc3, short outloc, double timer, double pred_t, double pred_J);

#endif
