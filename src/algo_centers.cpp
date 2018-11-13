// Algorithms for optimizing centers

// include rcpp

#include <Rcpp.h>

using namespace Rcpp;

// include headers
#include "util_general.h"
#include "util_dist.h"
#include "util_algo.h"
#include "algo_dist.h"
#include "algo_centers.h"

// functions

// [[Rcpp::export]]
NumericMatrix segment_wise_coords(NumericMatrix center, NumericMatrix path, double eps) {

  int size_c = center.nrow();
  int size_p = path.nrow();
  NumericMatrix new_coords(size_c*2, 3);

  FreeSpaceDiag fsp_cp = comp_free_space(center, path, eps);

  for(int i = 0; i < (size_c - 1); i++) {

    NumericVector keep_empty(size_p);
    NumericVector x(2);
    NumericVector y(2);

    for(int j = 0; j < (size_p - 1); j++) {

      if(is_empty_c(fsp_cp.L1, fsp_cp.L2, i, j)) {

        keep_empty[j] = 1;

        if(j > 0 && keep_empty[j-1] == 0) {

          break;

        }

      } else {

        keep_empty[j] = 0;

        if(j == 0 || keep_empty[j-1] == 0) {

          if(fsp_cp.L1(i,j) == 0) {

            x[0] = path(j,0);
            y[0] = path(j,1);

          } else {

            NumericVector tmp = get_coord(path(j,_), path(j+1,_), fsp_cp.L1(i,j));
            x[0] = tmp[0];
            y[0] = tmp[1];

          }

        }

        double uvdist = dist_c(path(j, 0), path(j, 1), path((j+1), 0), path((j+1), 1));

        if(uvdist == fsp_cp.L2(i,j)) {

          x[1] = path((j+1), 0);
          y[1] = path((j+1), 1);

        } else {

          NumericVector tmp = get_coord(path(j,_), path(j+1,_), fsp_cp.L2(i,j));

          x[1] = tmp[0];
          y[1] = tmp[1];

          break;

        }

      }

    }

    // fill with coordinates
    new_coords((i+1)*2-2, 0) = i + 1;
    new_coords((i+1)*2-1, 0) = i + 1;
    new_coords((i+1)*2-2, 1) = x[0];
    new_coords((i+1)*2-1, 1) = x[1];
    new_coords((i+1)*2-2, 2) = y[0];
    new_coords((i+1)*2-1, 2) = y[1];
  }

  new_coords(2*size_c - 2, 0) = size_c;
  new_coords(2*size_c - 1, 0) = size_c;
  new_coords(2*size_c - 2, 1) = path(size_p - 1, 0);
  new_coords(2*size_c - 1, 1) = path(size_p - 1, 0);
  new_coords(2*size_c - 2, 2) = path(size_p - 1, 1);
  new_coords(2*size_c - 1, 2) = path(size_p - 1, 1);

  return new_coords;

}

// calculate center of disk
//
// @param points a numeric matrix with an x and y column, filled with points
//
// [[Rcpp::export]]
NumericVector center_disk(NumericVector p1, NumericVector p2, NumericVector p3) {

  // init center
  NumericVector center(2);

  // p1 = p2 if by chance we get two times the same points
  if(p2[0] == p1[0] && p2[1] == p1[1]) {

    center[0] = p1[0];
    center[1] = p1[1];

    return center;
  }

  // p2 = p3 if 2 we have only boundary points
  if(p2[0] == p3[0] && p2[1] == p3[1]) {

    center[0] = (p2[0] + p1[0])/2.0;
    center[1] = (p2[1] + p1[1])/2.0;

  } else {

    double ma = (p2[1] - p1[1])/(p2[0] - p1[0]);
    double mb = (p3[1] - p2[1])/(p3[0] - p2[0]);

    double x = (ma*mb*(p1[1] - p3[1]) + mb*(p1[0] + p2[0]) - ma*(p2[0] + p3[0]))/(2*(mb-ma));
    double y = (-1/ma)*(x - (p1[0] + p2[0])/2.0) + (p1[1] + p2[1])/2.0;

    center[0] = x;
    center[1] = y;
  }

  return center;

}

//' Ramdomized smallest enclosing disk
//
//' @param points a numeric matrix with an x and y column, filled with points
//
//' @export
//
// [[Rcpp::export]]
NumericVector sed_randomized(NumericMatrix points) {

  int n = points.nrow();

  // check nr points
  if(n < 3) {

    if(n == 1) {

      NumericVector c(2);
      c[0] = points(0,0);
      c[1] = points(0,1);

      return c;

    } else {

      return center_disk(points(0,_), points(1,_), points(1,_));
    }

  }

  // Create permutation Pi
  IntegerVector Pi = index_shuffle(0, n);

  // initialize center and radius
  NumericVector center = center_disk(points(Pi[0],_), points(Pi[1],_), points(Pi[1],_));

  double radius = dist_c(center[0], center[1], points(Pi[0],0), points(Pi[0],1));


  // start algo
  for(int it_i = 2; it_i < n; it_i++) {

    int i = Pi[it_i];

    // if p_i not in circle (otherwise circle stays the same)
    if(dist_c(center[0], center[1], points(i,0), points(i,1)) > radius) {

      // set center SED p_i and p_0
      NumericVector center_p = center_disk(points(i,_), points(Pi[0],_), points(Pi[0],_));

      double radius_p = dist_c(center_p[0], center_p[1], points(Pi[0],0), points(Pi[0],1));

      // create permutation
      IntegerVector Pj = index_shuffle(1, it_i);

      // check points 1 to i-1 if in SED(p0, p_i)
      for(int it_j = 0; it_j < (it_i - 1); it_j++) {

        int j = Pi[Pj[it_j]];

        // if p_j not in SED(p0, p_i) (otherwise circle stays the same)
        if(dist_c(center_p[0], center_p[1], points(j,0), points(j,1)) > radius_p) {

          // set center SED p_i and p_j
          NumericVector center_q = center_disk(points(i,_), points(j,_), points(j,_));

          double radius_q = dist_c(center_q[0], center_q[1], points(j,0), points(j,1));

          // create permutation
          IntegerVector Pk = index_shuffle(0, it_j);

          for(int it_k = 0; it_k < it_j; it_k++) {

            int k = Pi[Pj[Pk[it_k]]];

            // if p_k not in SED(p_i, p_j) (otherwise circle stays the same)
            if(dist_c(center_q[0], center_q[1], points(k,0), points(k,1)) > radius_q) {

              // center SED(p_i, p_j, p_k)
              center_q = center_disk(points(i,_), points(j,_), points(k,_));
              radius_q = dist_c(center_q[0], center_q[1], points(j,0), points(j,1));

            }

          }

          center_p = center_q;
          radius_p = radius_q;
        }

      }

      center = center_p;
      radius = radius_p;

    }

  }

  return center;
}




//' Smallest enclosing disc
//'
//' @param points a numeric matrix with an x and y column, filled with points
//'
//' @export
//'
// [[Rcpp::export]]
NumericVector sed(NumericMatrix points) {

  int n = points.nrow();

  // check nr points
  if(n < 3) {

    if(n == 1) {

      NumericVector c(2);
      c[0] = points(0,0);
      c[1] = points(0,1);

      return c;

    } else {

      return center_disk(points(0,_), points(1,_), points(1,_));
    }

  }

  // initialize center and radius
  NumericVector center = center_disk(points(0,_), points(1,_), points(1,_));

  double radius = dist_c(center[0], center[1], points(0,0), points(0,1));


  // start algo
  for(int i = 2; i < n; i++) {

    // if p_i not in circle (otherwise circle stays the same)
    if(dist_c(center[0], center[1], points(i,0), points(i,1)) > radius) {

      // set center SED p_i and p_0
      NumericVector center_p = center_disk(points(i,_), points(0,_), points(0,_));

      double radius_p = dist_c(center_p[0], center_p[1], points(0,0), points(0,1));

      // check points 1 to i-1 if in SED(p0, p_i)
      for(int j = 1; j < i; j++) {

        // if p_j not in SED(p0, p_i) (otherwise circle stays the same)
        if(dist_c(center_p[0], center_p[1], points(j,0), points(j,1)) > radius_p) {

          // set center SED p_i and p_j
          NumericVector center_q = center_disk(points(i,_), points(j,_), points(j,_));

          double radius_q = dist_c(center_q[0], center_q[1], points(j,0), points(j,1));

          for(int k = 0; k < j; k++) {

            // if p_k not in SED(p_i, p_j) (otherwise circle stays the same)
            if(dist_c(center_q[0], center_q[1], points(k,0), points(k,1)) > radius_q) {

              // center SED(p_i, p_j, p_k)
              center_q = center_disk(points(i,_), points(j,_), points(k,_));
              radius_q = dist_c(center_q[0], center_q[1], points(j,0), points(j,1));

            }

          }

          center_p = center_q;
          radius_p = radius_q;
        }

      }

      center = center_p;
      radius = radius_p;

    }

  }

  return center;
}




//counter debug
//
// @param points a numeric matrix with an x and y column, filled with points
//
//
// [[Rcpp::export]]
List count_sed(NumericMatrix points, int times) {

  NumericVector sub1(times);
  NumericVector sub2(times);

  for(int z = 0; z < times; z++) {

    int counter1 = 0;
    int counter2 = 0;

    int n = points.nrow();

    // check nr points
    if(n < 3) {

      if(n == 1) {

        NumericVector c(2);
        c[0] = points(0,0);
        c[1] = points(0,1);

      }

    }

    // initialize center and radius
    NumericVector center = center_disk(points(0,_), points(1,_), points(1,_));

    double radius = dist_c(center[0], center[1], points(0,0), points(0,1));


    // start algo
    for(int i = 2; i < n; i++) {

      // if p_i not in circle (otherwise circle stays the same)
      if(dist_c(center[0], center[1], points(i,0), points(i,1)) > radius) {

        counter1 = counter1 + 1;
        // set center SED p_i and p_0
        NumericVector center_p = center_disk(points(i,_), points(0,_), points(0,_));

        double radius_p = dist_c(center_p[0], center_p[1], points(0,0), points(0,1));

        // check points 1 to i-1 if in SED(p0, p_i)
        for(int j = 1; j < i; j++) {

          // if p_j not in SED(p0, p_i) (otherwise circle stays the same)
          if(dist_c(center_p[0], center_p[1], points(j,0), points(j,1)) > radius_p) {

            counter2 = counter2 + 1;
            // set center SED p_i and p_j
            NumericVector center_q = center_disk(points(i,_), points(j,_), points(j,_));

            double radius_q = dist_c(center_q[0], center_q[1], points(j,0), points(j,1));

            for(int k = 0; k < j; k++) {

              // if p_k not in SED(p_i, p_j) (otherwise circle stays the same)
              if(dist_c(center_q[0], center_q[1], points(k,0), points(k,1)) > radius_q) {

                // center SED(p_i, p_j, p_k)
                center_q = center_disk(points(i,_), points(j,_), points(k,_));
                radius_q = dist_c(center_q[0], center_q[1], points(j,0), points(j,1));

              }

            }

            center_p = center_q;
            radius_p = radius_q;
          }

        }

        center = center_p;
        radius = radius_p;

      }

    }

    sub1[z] = counter1;
    sub2[z] = counter2;

  } // end iterator

  return List::create(sub1, sub2);

}



// counter debug random
//
// @param points a numeric matrix with an x and y column, filled with points
//
//
// [[Rcpp::export]]
List count_sed_random(NumericMatrix points, int times) {

  IntegerVector sub1(times);
  IntegerVector sub2(times);

  for(int z = 0; z < times; z++) {

    int counter1 = 0;
    int counter2 = 0;

    int n = points.nrow();

    // check nr points
    if(n < 3) {

      if(n == 1) {

        NumericVector c(2);
        c[0] = points(0,0);
        c[1] = points(0,1);

      }

    }

    // Create permutation Pi
    IntegerVector Pi = index_shuffle(0, n);

    // initialize center and radius
    NumericVector center = center_disk(points(Pi[0],_), points(Pi[1],_), points(Pi[1],_));

    double radius = dist_c(center[0], center[1], points(Pi[0],0), points(Pi[0],1));


    // start algo
    for(int it_i = 2; it_i < n; it_i++) {

      int i = Pi[it_i];

      // if p_i not in circle (otherwise circle stays the same)
      if(dist_c(center[0], center[1], points(i,0), points(i,1)) > radius) {

        counter1 = counter1 + 1;

        // set center SED p_i and p_0
        NumericVector center_p = center_disk(points(i,_), points(Pi[0],_), points(Pi[0],_));

        double radius_p = dist_c(center_p[0], center_p[1], points(Pi[0],0), points(Pi[0],1));

        // create permutation
        IntegerVector Pj = index_shuffle(1, it_i);

        // check points 1 to i-1 if in SED(p0, p_i)
        for(int it_j = 0; it_j < (it_i - 1); it_j++) {

          int j = Pi[Pj[it_j]];

          // if p_j not in SED(p0, p_i) (otherwise circle stays the same)
          if(dist_c(center_p[0], center_p[1], points(j,0), points(j,1)) > radius_p) {

            counter2 = counter2 + 1;

            // set center SED p_i and p_j
            NumericVector center_q = center_disk(points(i,_), points(j,_), points(j,_));

            double radius_q = dist_c(center_q[0], center_q[1], points(j,0), points(j,1));

            // create permutation
            IntegerVector Pk = index_shuffle(0, it_j);

            for(int it_k = 0; it_k < it_j; it_k++) {

              int k = Pi[Pj[Pk[it_k]]];

              // if p_k not in SED(p_i, p_j) (otherwise circle stays the same)
              if(dist_c(center_q[0], center_q[1], points(k,0), points(k,1)) > radius_q) {

                // center SED(p_i, p_j, p_k)
                center_q = center_disk(points(i,_), points(j,_), points(k,_));
                radius_q = dist_c(center_q[0], center_q[1], points(j,0), points(j,1));

              }

            }

            center_p = center_q;
            radius_p = radius_q;
          }

        }

        center = center_p;
        radius = radius_p;

      }

    }

    sub1[z] = counter1;
    sub2[z] = counter2;

  } // end iterator

  return List::create(sub1, sub2);

}
