///
/// \author Anastasiadis Petros (panastas@cslab.ece.ntua.gr)
///
/// \brief Some convinient C/C++ utilities for CHL.
///

#include "chl_smart_wrappers.hpp"
#include "chl_grid_amalgamation.hpp"
#include "microbenchmarks.hpp"

#include <numa.h>

double latbw_linear_regression(double* y, double* x, int samples, double* ans_b, double* ans_a, int active_unit_num){
 	/*	
	double sumX=0, sumX2=0, sumY=0, sumXY=0;
 	for(int i=0; i<samples; i++){
		sumX = sumX + x[i];
		sumX2 = sumX2 + x[i]*x[i];
		sumY = sumY + y[i];
		sumXY = sumXY + x[i]*y[i];
	}
 	//Calculating a and b 
	//*ans_a = (samples*sumXY-sumX*sumY)/(samples*sumX2-sumX*sumX); /// TODO: These result in negative latencies...
	//*ans_b = (sumY*sumX2-sumX*sumXY)/(samples*sumX2-sumX*sumX); //(sumY - (*ans_b)*sumX)/samples; /// TODO: And NNLS seems like an overkill...
	*/
	*ans_a = y[samples-1]/x[samples-1]; // Assuming the last sample has near-0 latency
	double sumb = 0; 
	int active_samples = samples; 
	for(int i=0; i<samples; i++){
		double rem_lat = y[i] - (*ans_a)*x[i];
		if(rem_lat >= 0.0) sumb += rem_lat;
		else active_samples--;
	}
	*ans_a = (y[samples-1]- sumb/active_samples)/x[samples-1];
	*ans_b = sumb/active_samples;
	//*ans_b = abs(sumY - (*ans_a)*sumX)/samples;
	double t_lat = *ans_b, t_b = *ans_a;
	fprintf(stderr,"Latency estimated by linear regresion :\t t_lat = %e\n", t_lat);
	fprintf(stderr,"BW-inv estimated by linear regresion:\t t_b = %e ( %lf Gb/s)\n", t_b, (1/t_b)/1e9);
	long byte_breakpoint = (int) (t_lat/t_b); 
	fprintf(stderr,"Latency-BW estimated breakpoint:\t b_eq = %ld bytes\n", byte_breakpoint/active_unit_num);
	fprintf(stderr,"75%% BW estimated breakpoint:\t\t b_75 = %ld bytes\n", byte_breakpoint/active_unit_num*2);
	fprintf(stderr,"95%% BW estimated breakpoint:\t\t b_95 = %ld bytes\n\n", byte_breakpoint/active_unit_num*20);
	fprintf(stderr,"Latency-BW estimated breakpoint:\t b_eq = 8-bytes X (%d X %d)\n", (int) sqrt(byte_breakpoint/active_unit_num/8), (int) sqrt(byte_breakpoint/active_unit_num/8));
	fprintf(stderr,"75%% BW estimated breakpoint:\t\t b_75 = 8-bytes X (%d X %d)\n", (int) sqrt(byte_breakpoint/active_unit_num/4), (int) sqrt(byte_breakpoint/active_unit_num/4));
	fprintf(stderr,"95%% BW estimated breakpoint:\t\t b_95 = 8-bytes X (%d X %d)\n\n", (int) sqrt(byte_breakpoint/active_unit_num*20/8), (int) sqrt(byte_breakpoint/active_unit_num*20/8));
	double APE = 0; 
	long long byte_lbw_bp = 0, byte_75pbw_bp = 0, byte_95pbw_bp = 0;
	for(int i=0; i<samples; i++){
		double pred_y =  t_lat + x[i]*t_b; 
#ifdef PDEBUG
		fprintf(stderr, "Point (x = %d):\t y=%lf ms (%lf Gb/s), pred_y = %lf ms (%lf Gb/s) PE = %.1lf\n", 
		int(x[i]), y[i] * 1000, Gval_per_s(int(x[i]), y[i]), pred_y * 1000, Gval_per_s(int(x[i]), pred_y), (pred_y - y[i])/y[i]*100);
#endif
		APE += abs((pred_y - y[i])/y[i]*100);
		if (!byte_lbw_bp && y[i]/x[i] <= (*ans_a)/0.5) byte_lbw_bp = x[i];
		if (!byte_75pbw_bp && y[i]/x[i] <= (*ans_a)/0.75) byte_75pbw_bp = x[i];
		if (!byte_95pbw_bp && y[i]/x[i] <= (*ans_a)/0.95) byte_95pbw_bp = x[i];

	}
	
	fprintf(stderr,"Latency-BW empirical breakpoint:\t b_eq = 8-bytes X (%d X %d)\n", (int) sqrt(byte_lbw_bp/active_unit_num/8), (int) sqrt(byte_lbw_bp/active_unit_num/8));
	fprintf(stderr,"75%% BW empirical breakpoint:\t\t b_75 = 8-bytes X (%d X %d)\n", (int) sqrt(byte_75pbw_bp/active_unit_num/8), (int) sqrt(byte_75pbw_bp/active_unit_num/8));
	fprintf(stderr,"95%% BW empirical breakpoint:\t\t b_95 = 8-bytes X (%d X %d)\n\n", (int) sqrt(byte_95pbw_bp/active_unit_num/8), (int) sqrt(byte_95pbw_bp/active_unit_num/8));

	return APE/samples; 
}
int confidence_interval_5_percent(long int sample_sz, double cpu_timer, double* transfer_t_vals, double* transfer_t_sum_ptr, double* transfer_t_mean_ptr, double* error_margin_ptr){
	transfer_t_vals[sample_sz-1] = cpu_timer;
	(*transfer_t_sum_ptr) += transfer_t_vals[sample_sz-1];
	(*transfer_t_mean_ptr) = (*transfer_t_sum_ptr)/sample_sz;
	if (sample_sz < 2) return 0;
	double std_dev = 0; 
	for (int i = 0; i < sample_sz; i++) std_dev += pow(transfer_t_vals[i] - (*transfer_t_mean_ptr), 2);
	std_dev /= sample_sz;
		std_dev = sqrt(std_dev);
	boost::math::students_t dist(sample_sz - 1);
	double Td = boost::math::quantile(boost::math::complement(dist, alphaCI / 2));
	(*error_margin_ptr) = Td*std_dev/sqrt(sample_sz);
#ifdef DPDEBUG
	fprintf(stderr, "\tItter %ld:\t mean=%lf, std_dev = %lf, Error margin =%lf\n", sample_sz, (*transfer_t_mean_ptr) , std_dev, (*error_margin_ptr));
#endif
	if ((*error_margin_ptr)/(*transfer_t_mean_ptr)  * 100 <= 0.5) return 1;
	else return 0;
}