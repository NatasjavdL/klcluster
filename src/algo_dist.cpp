// Distance algorithms - frechet - dtw

// include headers
#include "util_general.h"
#include "util_dist.h"
#include "util_algo.h"
#include "algo_dist.h"

// include rcpp

#include <Rcpp.h>

using namespace Rcpp;

// forward declarations
struct FreeSpaceDiag;

// functions


// subroutine computing free space diagram cont. frechet
//
// @param avec numericvector corresponding to one point on the plane
// avec[0] is the x-coordinate
// avec[1] is the y-coordinate
//
// @param uvec - see avec
// @param vvec - see avec
//
//
// [[Rcpp::export]]
NumericVector get_free_interval_c(NumericVector avec, NumericVector uvec,
  NumericVector vvec, double eps) {

  NumericVector intv = get_interval_c();

  double dist_a_to_uv = dist_point_to_line_c(avec, uvec, vvec);

  if(dist_a_to_uv < eps) {

    double uv_dist = dist_c(uvec[0], uvec[1], vvec[0], vvec[1]);

    double root = 0;

    if(!((vvec[0] == uvec[0]) && (vvec[1] == uvec[1]))) {

      NumericVector s1 = subtract_xy_c(avec, uvec);
      NumericVector s2 = subtract_xy_c(vvec, uvec);

      root = calc_dot(s1[0], s1[1], s2[0], s2[1]);
      root = root / uv_dist;

    }

    double offset = sqrt(get_max2((pow(eps, 2) - pow(dist_a_to_uv, 2)), 0));

    double minim = root - offset;
    if (minim < 0) {
      minim = 0;
    }

    double maxim = root + offset;
    if (maxim > uv_dist) {
      maxim = uv_dist;
    }

    intv[0] = minim;
    intv[1] = maxim;


  }

  return(intv);
}

// [[Rcpp::export]]
bool use_interval_c(NumericMatrix p1, NumericMatrix p2, double eps) {

  int size_p1 = p1.nrow();
  int size_p2 = p2.nrow();

  NumericMatrix LR_min(size_p1, size_p2);
  NumericMatrix LR_max(size_p1, size_p2);
  NumericMatrix BR_min(size_p1, size_p2);
  NumericMatrix BR_max(size_p1, size_p2);

  LR_min = LR_min + 1;
  BR_min = BR_min + 1;

  for(int i = 0; i < (size_p1 - 1); i++) {

    for(int j = 0; j < (size_p2 - 1); j++) {

      NumericVector L;
      NumericVector B;

      if(i == 0) {

        if(j == 0) {

          L = get_free_interval_c(p1(i,_), p2(j,_), p2(j+1,_), eps);

        } else {

          L = get_interval_c();

        }

      } else if(!is_empty_c(BR_min, BR_max, (i-1), j)) {

        L = get_free_interval_c(p1(i,_), p2(j,_), p2(j+1,_), eps);

      } else if(!is_empty_c(LR_min, LR_max, (i-1), j)) {

        L = get_free_interval_c(p1(i,_), p2(j,_), p2(j+1,_), eps);

        L[0] = get_max2(L[0], LR_min(i-1, j));

      } else {

        L = get_interval_c();

      }

      // add L to LR
      LR_min(i, j) = L[0];
      LR_max(i, j) = L[1];

      if(j == 0) {

        if(i == 0) {

          B = get_free_interval_c(p2(j,_), p1(i,_), p1(i+1,_), eps);

        } else {

          B = get_interval_c();

        }

      } else if(!is_empty_c(LR_min, LR_max, i, (j-1))) {

        B = get_free_interval_c(p2(j,_), p1(i,_), p1(i+1,_), eps);

      } else if(!is_empty_c(BR_min, BR_max, i, (j-1))) {

        B = get_free_interval_c(p2(j,_), p1(i,_), p1(i+1,_), eps);

        B[0] = get_max2(B[0], BR_min(i, j-1));

      } else {

        B = get_interval_c();

      }
      // add B to BR
      BR_min(i, j) = B[0];
      BR_max(i, j) = B[1];

    }
  }

  return(!is_empty_c(LR_min, LR_max, (size_p1 - 2), (size_p2 - 2)) ||
    !is_empty_c(BR_min, BR_max, (size_p1 - 2), (size_p2 - 2)));
}

// [[Rcpp::export]]
bool check_eps_c(NumericMatrix p1_mat, NumericMatrix p2_mat, double eps) {

  int size_p1 = p1_mat.nrow();
  int size_p2 = p2_mat.nrow();

  if(dist_c(p1_mat(0,0), p1_mat(0,1), p2_mat(0,0), p2_mat(0,1)) >= eps) {

    return(FALSE);

  }

  if(dist_c(p1_mat((size_p1-1),0), p1_mat((size_p1-1),1),
    p2_mat((size_p2-1),0), p2_mat((size_p2-1),1)) >= eps) {

    return(FALSE);

  }

  if(size_p1 == 1) {

    NumericVector dists(size_p2-1);

    for(int i = 0; i < size_p2; i++) {

      dists[i] = dist_c(p1_mat(0,0), p1_mat(0,1), p2_mat(i,0), p2_mat(i,1));

    }

    if(get_max(dists) >= eps) {
      return(FALSE);
    }

  }

  if(size_p2 == 1) {

    NumericVector dists(size_p1-1);

    for(int i = 0; i < size_p1; i++) {

      dists[i] = dist_c(p2_mat(0,0), p2_mat(0,1), p1_mat(i,0), p1_mat(i,1));

    }

    if(get_max(dists) >= eps) {
      return(FALSE);
    }

  }

  return(use_interval_c(p1_mat, p2_mat, eps));

}

// [[Rcpp::export]]
double get_frechet_c(NumericMatrix p1_mat, NumericMatrix p2_mat, double precision) {

  int size_p1 = p1_mat.nrow();
  int size_p2 = p2_mat.nrow();

  if(size_p1 == 1) {

    NumericVector dists = calc_dist_c(p1_mat, p2_mat);

    return get_max(dists);

  }

  if(size_p2 == 1) {

    NumericVector dists = calc_dist_c(p1_mat, p2_mat);

    return get_max(dists);

  }

  // initialize lowerbound and upperbound

  double lowerbound = get_max2(dist_c(p1_mat(0,0), p1_mat(0,1),
    p2_mat(0,0), p2_mat(0,1)), dist_c(p1_mat(size_p1,0), p1_mat(size_p1,1),
    p2_mat(size_p2,0), p2_mat(size_p2,1)));

  if(lowerbound < 0) lowerbound = 0;

  double upperbound = 0;

  for(int i = 0; i < get_max2(size_p1, size_p2); i++) {

    double it_p1 = get_min2(i, size_p1);
    double it_p2 = get_min2(i, size_p2);

    upperbound = get_max2(upperbound, dist_c(p1_mat(it_p1,0), p1_mat(it_p1,1),
      p2_mat(it_p2,0), p2_mat(it_p2,1)));

  }


  // binary search on epsilon

  while(lowerbound < upperbound - precision) {

    double mid = (lowerbound + upperbound) / 2;

    if(check_eps_c(p1_mat, p2_mat, mid)) {

      upperbound = mid;

    } else {

      lowerbound = mid;

    }

  }

  if(upperbound == INFINITY) {

    std::cout<<("warning: distance equals infinity")<<"\n";

  }
  return upperbound;

}

// [[Rcpp::export]]
NumericMatrix dtw(NumericMatrix p1, NumericMatrix p2, bool sqr_dist) {

  // p1 corresponds to the rows, p2 to the columns
  int size_p1 = p1.nrow();
  int size_p2 = p2.nrow();

  NumericMatrix mat(size_p1, size_p2);

  double inf = 10000000000000;

  // initialize values to infinity
  for(int i = 0; i < size_p1; i++) {

    for(int j = 0; j < size_p2; j++) {

      mat(i,j) = inf;

    }

  }

  for(int i = 0; i < size_p1; i++) {

    for(int j = 0; j < size_p2; j++) {

      double cost = squared_dist(p1(i,0), p1(i,1), p2(j,0), p2(j,1));

      if(!sqr_dist) {

        cost = sqrt(cost);

      }

      if(i == 0) {

        if(j == 0) {

          mat(i,j) = cost;

        } else {

          mat(i,j) = cost + mat(i, j-1);

        }

      } else if(j == 0) {

        mat(i,j) = cost + mat(i-1, j);

      } else {

        NumericVector values(3);

        values[0] = mat(i-1, j);
        values[1] = mat(i, j-1);
        values[2] = mat(i-1, j-1);

        mat(i,j) = cost + get_min(values);

      }

    }

  }

  return mat;

}

// [[Rcpp::export]]
double calc_ca(NumericMatrix ca, int i, int j, NumericMatrix p1, NumericMatrix p2) {

  if(ca(i,j) > -1) {

    return ca(i,j);

  } else if(i == 0 && j == 0) {

    ca(i,j) = dist_c(p1(i,0), p1(i,1), p2(j,0), p2(j,1));

  } else if(i > 0 && j == 0) {

    ca(i,j) = get_max2(calc_ca(ca, (i-1), j, p1, p2),
      dist_c(p1(i,0), p1(i,1), p2(j,0), p2(j,1)));

  } else if(i == 0 &&  j > 0) {

    get_max2(calc_ca(ca, i, (j-1), p1, p2),
      dist_c(p1(i,0), p1(i,1), p2(j,0), p2(j,1)));

  } else if(i > 0 && j > 0) {

    NumericVector values(3);
    values[0] = calc_ca(ca, (i-1), j, p1, p2);
    values[1] = calc_ca(ca, (i-1), (j-1), p1, p2);
    values[2] = calc_ca(ca, i, (j-1), p1, p2);

    ca(i,j) = get_max2(get_min(values), dist_c(p1(i,0), p1(i,1), p2(j,0), p2(j,1)));

  } else {

    ca(i,j) = 10000000000;

  }

  return ca(i,j);

}

FreeSpaceDiag comp_free_space(NumericMatrix p1, NumericMatrix p2, double eps) {

  int size_p1 = p1.nrow();
  int size_p2 = p2.nrow();

  NumericMatrix LR_min(size_p1, size_p2);
  NumericMatrix LR_max(size_p1, size_p2);
  NumericMatrix BR_min(size_p1, size_p2);
  NumericMatrix BR_max(size_p1, size_p2);

  LR_min = LR_min + 1;
  BR_min = BR_min + 1;

  for(int i = 0; i < (size_p1 - 1); i++) {

    for(int j = 0; j < (size_p2 - 1); j++) {

      NumericVector L;
      NumericVector B;

      if(i == 0) {

        if(j == 0) {

          L = get_free_interval_c(p1(i,_), p2(j,_), p2(j+1,_), eps);

        } else {

          L = get_interval_c();

        }

      } else if(!is_empty_c(BR_min, BR_max, (i-1), j)) {

        L = get_free_interval_c(p1(i,_), p2(j,_), p2(j+1,_), eps);

      } else if(!is_empty_c(LR_min, LR_max, (i-1), j)) {

        L = get_free_interval_c(p1(i,_), p2(j,_), p2(j+1,_), eps);

        L[0] = get_max2(L[0], LR_min(i-1, j));

      } else {

        L = get_interval_c();

      }

      // add L to LR
      LR_min(i, j) = L[0];
      LR_max(i, j) = L[1];

      if(j == 0) {

        if(i == 0) {

          B = get_free_interval_c(p2(j,_), p1(i,_), p1(i+1,_), eps);

        } else {

          B = get_interval_c();

        }

      } else if(!is_empty_c(LR_min, LR_max, i, (j-1))) {

        B = get_free_interval_c(p2(j,_), p1(i,_), p1(i+1,_), eps);

      } else if(!is_empty_c(BR_min, BR_max, i, (j-1))) {

        B = get_free_interval_c(p2(j,_), p1(i,_), p1(i+1,_), eps);

        B[0] = get_max2(B[0], BR_min(i, j-1));

      } else {

        B = get_interval_c();

      }
      // add B to BR
      BR_min(i, j) = B[0];
      BR_max(i, j) = B[1];

    }
  }

  FreeSpaceDiag fsp;
  fsp.L1 = LR_min;
  fsp.L2 = LR_max;
  fsp.B1 = BR_min;
  fsp.B2 = BR_max;

  return fsp;

}

