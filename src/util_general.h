// header guard

#ifndef UTIL_GEN_H
#define UTIL_GEN_H

// include rcpp

#include <Rcpp.h>

using namespace Rcpp;

// function prototypes

int which_min(NumericVector x);

NumericVector subset_vec(NumericVector x, IntegerVector index);

double get_max(NumericVector x);

double get_min(NumericVector x);

double get_max2(double x, double y);

double get_min2(double x, double y);

NumericVector df1_to_nv(DataFrame df);

NumericMatrix df_to_mat(DataFrame x);

DataFrame mat_to_df(NumericMatrix mat, bool path);

NumericVector concat_vec(NumericVector vec1, NumericVector vec2);

NumericMatrix subset_mat(NumericMatrix mat, NumericVector index);


// end of header file
#endif