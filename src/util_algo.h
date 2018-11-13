// header guard

#ifndef UTIL_ALGO_H
#define UTIL_ALGO_H

// include rcpp

#include <Rcpp.h>

using namespace Rcpp;

struct FreeSpaceDiag {

  NumericMatrix L1;
  NumericMatrix L2;
  NumericMatrix B1;
  NumericMatrix B2;

};

// function prototypes

NumericVector get_interval_c();

bool is_empty_c(NumericMatrix xmin, NumericMatrix xmax, int row, int col);

NumericMatrix dist_matrix(NumericMatrix p1, NumericMatrix p2);

NumericMatrix backtrack(NumericMatrix center, NumericMatrix path, NumericMatrix dist_mat);

NumericVector get_coord(NumericVector point1, NumericVector point2, double dist);

IntegerVector index_shuffle(int start, int n);


// end of header file
#endif
