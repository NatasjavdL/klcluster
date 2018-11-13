// header guard

#ifndef UTIL_DIST_H
#define UTIL_DIST_H

// include rcpp

#include <Rcpp.h>

using namespace Rcpp;

// function prototypes

double dist_c(double ax, double ay, double bx, double by);

double squared_dist(double ax, double ay, double bx, double by);

double calc_dot(double ax, double ay, double bx, double by);

NumericVector coord_wise_mean(NumericMatrix coords, double nr_coords);

NumericVector calc_dist_c(NumericMatrix a, NumericMatrix b);

NumericVector subtract_xy_c(NumericVector u, NumericVector v);

NumericVector dot_product_c(NumericMatrix a, NumericMatrix b);

double dist_point_to_line_c(NumericVector a, NumericVector u, NumericVector v);


// end of header file
#endif