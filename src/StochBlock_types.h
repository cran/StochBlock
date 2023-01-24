#include <memory>
#include <float.h>
#include <algorithm>
#include <RcppArmadillo.h>

using DVector = Rcpp::NumericVector;
using IVector = Rcpp::IntegerVector;
using Array   = arma::Cube<double>; // [ Rows, Cols, Slices ]
using DMatrix = arma::Mat<double>;

enum Diagonale
{
    Same     = 0,
    Ignore   = 1,
    Seperate = 2
};

enum BorderType
{
    None   = 0,
    Inside = 1,
    Outside = 2,
};
