// Supporting functions algorithms

// include headers
#include "util_general.h"
#include "util_dist.h"
#include "util_algo.h"

// include rcpp

#include <Rcpp.h>

using namespace Rcpp;


// functions

// [[Rcpp::export]]
NumericVector get_interval_c() {
  NumericVector vec(2);

  vec[0] = 1;
  vec[1] = 0;

  return(vec);
}

// [[Rcpp::export]]
bool is_empty_c(NumericMatrix xmin, NumericMatrix xmax, int row, int col) {

  return(xmin(row, col) > xmax(row, col));

}


// [[Rcpp::export]]
NumericMatrix dist_matrix(NumericMatrix p1, NumericMatrix p2) {

  // p1 corresponds to the rows, p2 to the columns
  int size_p1 = p1.nrow();
  int size_p2 = p2.nrow();

  NumericMatrix mat(size_p1, size_p2);

  for(int i = 0; i < size_p2; i++) {

    NumericVector tmp(size_p1);

    for(int j = 0; j < size_p1; j++) {

      tmp[j] = dist_c(p1(j,0), p1(j,1), p2(i,0), p2(i,1));

    }

    mat(_, i) = tmp;

  }

  return mat;

}

// [[Rcpp::export]]
NumericVector get_coord(NumericVector point1, NumericVector point2, double dist) {

  double ax = point1[0];
  double ay = point1[1];
  double bx = point2[0];
  double by = point2[1];

  double d = sqrt(pow(bx-ax, 2) + pow(by-ax, 2));

  NumericVector new_vec(2);

  new_vec[0] = ax + (dist/d) * (bx-ax);
  new_vec[1] = ay + (dist/d) * (by-ay);

  return new_vec;

}

// [[Rcpp::export]]
NumericMatrix backtrack(NumericMatrix center, NumericMatrix path,
  NumericMatrix dist_mat) {

  int size_c = center.nrow();
  int size_p = path.nrow();

  NumericVector vec_c(size_p + size_c);
  NumericVector vec_px(size_p + size_c);
  NumericVector vec_py(size_p + size_c);

  int counter = 0;

  vec_c[counter] = size_c - 1;
  vec_px[counter] = path(size_p - 1, 0);
  vec_py[counter] = path(size_p - 1, 1);


  int i = size_c - 1;
  int j = size_p - 1;

  while(i > 0 || j > 0) {

    counter += 1;

    if(i == 0) {

      j = j - 1;

    } else if(j == 0) {

      i = i - 1;

    } else {

      NumericVector tmp(3);
      tmp[0] = dist_mat(i-1, j-1);
      tmp[1] = dist_mat(i, j-1);
      tmp[2] = dist_mat(i-1, j);


      int index = which_min(tmp);

      if(index == 0) {

        i = i - 1;
        j = j - 1;

      } else if(index == 1) {

        j = j - 1;

      } else {

        i = i - 1;

      }

    }

    vec_c[counter] = i;
    vec_px[counter] = path(j, 0);
    vec_py[counter] = path(j, 1);

  }

  NumericMatrix points(counter + 1, 3);

  for(int i = 0; i < counter + 1; i++) {

    points(i,0) = vec_c[i] + 1;
    points(i,1) = vec_px[i];
    points(i,2) = vec_py[i];

  }


  return points;
}


inline int randWrapper(const int n) {

  GetRNGstate();
  int tmp = floor(unif_rand()*n);
  PutRNGstate();

  return tmp;

}


// Create index permutation
//
// @param integer x - number of elements in vectors
//
// [[Rcpp::export]]
IntegerVector index_shuffle(int start, int n) {

  IntegerVector shuffle_vec(n);


  for (int i = start; i < n; i++) shuffle_vec[i] = i;

  std::random_shuffle(shuffle_vec.begin(), shuffle_vec.end(), randWrapper);

  return shuffle_vec;
}

