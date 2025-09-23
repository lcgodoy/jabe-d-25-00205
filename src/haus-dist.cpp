#include <RcppArmadillo.h>
#include <cmath> // For sin, cos, atan2, sqrt, pow

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

const double EARTH_RADIUS = 6371.0;
const double PI = 3.14159265358979323846;

// Helper function to convert degrees to radians
inline double deg2rad(double degrees) {
  return degrees * PI / 180.0;
}

// [[Rcpp::export]]
double gc_dist(const arma::rowvec& p1, const arma::rowvec& p2) {
  // Assumes p1 and p2 are rows with [longitude, latitude] in degrees
  double lon1_rad = deg2rad(p1[0]);
  double lat1_rad = deg2rad(p1[1]);
  double lon2_rad = deg2rad(p2[0]);
  double lat2_rad = deg2rad(p2[1]);
  double d_lon = lon2_rad - lon1_rad;
  double d_lat = lat2_rad - lat1_rad;
  double a = std::pow(std::sin(d_lat / 2.0), 2.0) +
    std::cos(lat1_rad) * std::cos(lat2_rad) *
    std::pow(std::sin(d_lon / 2.0), 2.0);
  double c = 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));
  return EARTH_RADIUS * c;
}

// [[Rcpp::export]]
arma::mat crossdist(const arma::mat& x, const arma::mat& y) {
  int nrow1 = x.n_rows;
  int nrow2 = y.n_rows;
  arma::mat out(nrow1, nrow2);

  for (int i = 0; i < nrow1; i++) {
    for (int j = 0; j < nrow2; j++) {
      out(i, j) = gc_dist(x.row(i), y.row(j));
    }
  }
  return out;
}

// [[Rcpp::export]]
double haus_aux(const arma::mat& x, const arma::mat& y) {
  arma::mat aux_mat = crossdist(x, y);
  double h_X_Y = arma::max(arma::min(aux_mat, 1));
  double h_Y_X = arma::max(arma::min(aux_mat, 0));
  return std::max(h_X_Y, h_Y_X);
}

// [[Rcpp::export]]
arma::mat dist_haus(const std::vector<arma::mat>& poly_list) {
  int n = poly_list.size();
  arma::mat out(n, n, arma::fill::zeros);

  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      out(i, j) = haus_aux(poly_list[i], poly_list[j]);
      out(j, i) = out(i, j); // The matrix is symmetric
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::mat cdist_haus_mp(const std::vector<arma::mat>& poly_list,
                        const arma::mat& pred_mat) {
  int n1 = poly_list.size();
  int n2 = pred_mat.n_rows;
  arma::mat out(n1, n2);

  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      // Treat each point in pred_mat as a polygon with one vertex
      out(i, j) = haus_aux(poly_list[i], pred_mat.row(j));
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::mat cdist_haus_lists(const std::vector<arma::mat>& poly_list,
                           const std::vector<arma::mat>& poly_pred) {
  int n1 = poly_list.size();
  int n2 = poly_pred.size();
  arma::mat out(n1, n2);

  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      out(i, j) = haus_aux(poly_list[i], poly_pred[j]);
    }
  }
  return out;
}
