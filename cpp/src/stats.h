#pragma once
#include <vector>
#include "pop_data.h"


void reg_with_intercept(const std::vector<double>& x,
                        const std::vector<double>& y,
                        double& a,
                        double& b);

double mean(const std::vector<double>& x);
double sd(const std::vector<double>& x);
double t_pvalue(double t, int df);
double sd_pop(const std::vector<double>& x);  // ddof=0

double se_pop(const std::vector<double>& x);                 // sd_pop(x)/sqrt(n)
void ci95_mean_pop(const std::vector<double>& x, double& lo, double& hi);  // mean ± tcrit*SE
                                                            
//double corr_prs_pheno(const std::vector<PopPoint>& v);

double corr_sums(double n, double sumx, double sumy, double sumxx, double sumyy, double sumxy);




bool qc_middle_bins(const std::vector<BinStats>& bins, double& min_p, int& min_bin, double& bonf_p); 
