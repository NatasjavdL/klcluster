// General functions

// include headers
#include "util_general.h"

// include rcpp

#include <Rcpp.h>

using namespace Rcpp;

// functions

// [[Rcpp::export]]
int which_min(NumericVector x) {

  int size_x = x.size();

  int index = 0;

  double value = x[0];

  for(int i = 0; i < size_x; i++) {

    if(x[i] < value) {

      value = x[i];
      index = i;

    }

  }

  return index;

}

// [[Rcpp::export]]
NumericVector subset_vec(NumericVector x, IntegerVector index) {

  int size_index = index.size();

  NumericVector new_x(size_index);

  index = index - 1;

  for (int i = 0; i < size_index; i++) {
    new_x[i] = x[index[i]];
  }

  return new_x;
}

// [[Rcpp::export]]
double get_max(NumericVector x) {

  int n = x.size();

  double max = 0;

  for(int i = 0; i < n; i++) {

    if(max < x[i]) {
      max = x[i];
    }
  }

  return(max);

}

// [[Rcpp::export]]
double get_min(NumericVector x) {

  int n = x.size();

  double min = x[0];

  for(int i = 0; i < n; i++) {

    if(min > x[i]) {
      min = x[i];
    }
  }

  return(min);

}

// [[Rcpp::export]]
double get_max2(double x, double y) {

  if(x < y) {
    return(y);
  }

  if(x > y) {
    return(x);
  }

  return(x);
}

// [[Rcpp::export]]
double get_min2(double x, double y) {

  if(x < y) {
    return(x);
  }

  if(x > y) {
    return(y);
  }

  return(x);

}

// [[Rcpp::export]]
NumericVector df1_to_nv(DataFrame df) {
  NumericVector nv;

  nv[0] = df["x"];
  nv[1] = df["y"];

  return(nv);

}

// transform dataframe to matrix
// [[Rcpp::export]]
NumericMatrix df_to_mat(DataFrame x) {
  int nRows=x.nrows();
  NumericMatrix y(nRows,x.size());

  for (int i=0; i<x.size();i++) {

    y(_,i)=NumericVector(x[i]);

  }

  return(y);
}

// transform matrix to dataframe
//
// [[Rcpp::export]]
DataFrame mat_to_df(NumericMatrix mat, bool path) {

  if(path) {

    return DataFrame::create(_["long"]= mat(_,0), _["lat"]= mat(_,1));

  } else {

    int size_m = mat.ncol();

    List tmp(size_m);

    for(int i = 0; i < size_m; i++) {

      tmp[i] = mat(_, i);

    }

    Rcpp::DataFrame df(tmp);
    return df;

  }

}

// [[Rcpp::export]]
NumericVector concat_vec(NumericVector vec1, NumericVector vec2) {

  int size_v1 = vec1.size();
  int size_v2 = vec2.size();

  NumericVector new_vec(size_v1 + size_v2);

  for(int i = 0; i < size_v1; i++) {

    new_vec[i] = vec1[i];

  }

  for(int j = 0; j < size_v2; j++) {

    new_vec[j + size_v1] = vec2[j];

  }

  return new_vec;

}

// [[Rcpp::export]]
NumericMatrix subset_mat(NumericMatrix mat, NumericVector index) {

  int size_ind = index.size();

  NumericMatrix new_mat(size_ind, mat.ncol());

  for(int i = 0; i < size_ind; i++) {

    new_mat(i,_) = mat(index[i], _);

  }

  return new_mat;

}
