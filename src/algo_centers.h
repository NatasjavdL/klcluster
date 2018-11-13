// header guard

#ifndef ALGO_CENTERS_H
#define ALGO_CENTERS_H

// include rcpp

#include <Rcpp.h>

using namespace Rcpp;


// function prototypes

NumericMatrix segment_wise_coords(NumericMatrix center, NumericMatrix path, double eps);

NumericVector center_disk(NumericVector p1, NumericVector p2, NumericVector p3);

NumericVector sed(NumericMatrix points);

List count_sed(NumericMatrix points, int times);

List count_sed_random(NumericMatrix points, int times);


// end of header file
#endif
