// General distance functions

// include headers
#include "util_dist.h"
#include "util_general.h"

// include rcpp

#include <Rcpp.h>

using namespace Rcpp;

// functions

// calculate euclidean distance between two points (a,b) on a plane
//
// @param ax the x coordinate of point a
//
// @param ay the y coordinate of point a
//
// @param bx the x coordinate of point b
//
// @param by the y coordinate of point b
//
// [[Rcpp::export]]
double dist_c(double ax, double ay, double bx, double by) {

  double dx = ax - bx;
  double dy = ay - by;

  return(sqrt(pow(dx, 2) + pow(dy, 2)));
}

// calculate the squared distance between two points (a,b) on a plane
//
// @param ax the x coordinate of point a
//
// @param ay the y coordinate of point a
//
// @param bx the x coordinate of point b
//
// @param by the y coordinate of point b
//
// [[Rcpp::export]]
double squared_dist(double ax, double ay, double bx, double by) {

  double dx = ax - bx;
  double dy = ay - by;

  return(pow(dx, 2) + pow(dy, 2));
}

// calculate the dot product between two points (a,b) on a plane
//
// @param ax the x coordinate of point a
//
// @param ay the y coordinate of point a
//
// @param bx the x coordinate of point b
//
// @param by the y coordinate of point b
//
// [[Rcpp::export]]
double calc_dot(double ax, double ay, double bx, double by) {

  return(ax * bx + ay * by);

}

// calculate the coordinate wise mean between a set of points on the plane
//
// @param coords a matrix, consisting of 2 columns, one with x coordinates,
// the second with the y coordinates
//
// @param nr_coords, specify the number of coordinates to calculate the mean
// over. Usefull in case this does not equal the number of points in coords.
//
// [[Rcpp::export]]
NumericVector coord_wise_mean(NumericMatrix coords, double nr_coords) {

  NumericVector x = coords(_, 0);
  NumericVector y = coords(_, 1);

  double new_x = std::accumulate(x.begin(), x.end(), 0) / nr_coords;
  double new_y = std::accumulate(y.begin(), y.end(), 0) / nr_coords;

  NumericVector new_coords(2);
  new_coords[0] = new_x;
  new_coords[1] = new_y;

  return new_coords;

}

// [[Rcpp::export]]
NumericVector calc_dist_c(NumericMatrix a, NumericMatrix b) {
  NumericVector ax = a(_,0);
  NumericVector ay = a(_,1);
  NumericVector bx = b(_,0);
  NumericVector by = b(_,1);

  int size_a = ax.size();
  int size_b = bx.size();

  int n = get_max2(size_a, size_b);

  NumericVector new_vec(n);

  for(int i = 0; i < n; i++) {
    int it_a = get_min2(size_a, i);
    int it_b = get_min2(size_b, i);

    new_vec[i] = dist_c(ax[it_a], ay[it_a], bx[it_b], by[it_b]);

  }

  return(new_vec);

}

// [[Rcpp::export]]
NumericVector subtract_xy_c(NumericVector u, NumericVector v) {
  NumericVector new_vec;

  new_vec = u-v;

  return(new_vec);

}

// [[Rcpp::export]]
NumericVector dot_product_c(NumericMatrix a, NumericMatrix b) {
  NumericVector ax = a(_,0);
  NumericVector ay = a(_,1);
  NumericVector bx = b(_,0);
  NumericVector by = b(_,1);

  int n = ax.size();

  NumericVector new_vec(n);

  for(int i = 0; i < n; i++) {

    new_vec[i] = calc_dot(ax[i], ay[i], bx[i], by[i]);

  }

  return(new_vec);

}

// [[Rcpp::export]]
double dist_point_to_line_c(NumericVector a, NumericVector u, NumericVector v) {

  if((a[0] == u[0] && a[1] == u[1]) || (a[0] == v[0] && a[1] == v[1])) {

    return(0);

  }

  if(v[0] == u[0] && v[1] == u[1]) {

    return(dist_c(a[0], a[1], u[0], u[1]));

  }

  // calculate sides triangle
  double uv = dist_c(u[0], u[1], v[0], v[1]);
  double au = dist_c(a[0], a[1], u[0], u[1]);
  double av = dist_c(a[0], a[1], v[0], v[1]);

  double error = 0.00000000001;

  // calculate area
  double s = 0.5 * (au + av + uv);
  double area = sqrt(s * (s-au + error) * (s-av + error) * (s-uv + error));

  double height = 2*area/uv;

  return(height);

}
