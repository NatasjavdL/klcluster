// header guard

#ifndef ALGO_DIST_H
#define ALGO_DIST_H

#include "util_algo.h"

// include rcpp

#include <Rcpp.h>

using namespace Rcpp;

// function prototypes

NumericVector get_free_interval_c(NumericVector avec, NumericVector uvec, NumericVector vvec, double eps);

FreeSpaceDiag comp_free_space(NumericMatrix p1, NumericMatrix p2, double eps);

bool use_interval_c(NumericMatrix p1, NumericMatrix p2, double eps);

bool check_eps_c(NumericMatrix p1_mat, NumericMatrix p2_mat, double eps);

double get_frechet_c(NumericMatrix p1_mat, NumericMatrix p2_mat, double precision);

double calc_ca(NumericMatrix ca, int i, int j, NumericMatrix p1, NumericMatrix p2);

NumericMatrix dtw(NumericMatrix p1, NumericMatrix p2, bool sqr_dist);



// end of header file
#endif
