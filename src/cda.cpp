//' @useDynLib trlasso
//' @importFrom Rcpp sourceCpp evalCpp

#include <Rcpp.h>
#include <Rinternals.h>
#include <RcppCommon.h>
#include <Rcpp/Benchmark/Timer.h>
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins("cpp11")]]
// #include <boost/numeric/ublas/matrix.hpp>
using namespace Rcpp;

//' calculate covariance matrix
//'
//' @param X design matrix
//' @return covariance matrix
//' @keywords internal
//'
// [[Rcpp::export]]
NumericMatrix covC(NumericMatrix X) {
  int n = X.nrow();
  int p = X.ncol();
  NumericMatrix out(p, p);
  for (int i = 0; i < p; i++) {
    for (int j = i; j < p; j++) {
      double total = 0;
      int num = 0;
      for (int k = 0; k < n; k++) {
        if (!R_IsNA(X(k, i)) & !R_IsNA(X(k, j))) {
          total += X(k, i) * X(k, j);
          num += 1;
        }
      }
      out(i, j) = total / num;
    }
  }
  for (int i = 0; i < p; i++) {
    for (int j = 0; j < i; j++) {
      out(i, j) = out(j, i);
    }
  }
  return out;
}

//' soft thresholding function for transfer lasso
//'
//' @param z z
//' @param g1 gamma1
//' @param g2 gamma2
//' @param b b
//' @return value
//' @keywords internal
//'
// [[Rcpp::export]]
double softThresholdC2(double z, double g1, double g2, double b) {
  if (b >= 0) {
    if ((z >= -g1) && (z <= g2)) {
      return 0;
    } else if ((z >= (g2 + b)) && (z <= (g1 + b))) {
      return b;
    } else if ((z >= g2) && (z <= (g2 + b)) ) {
      return z - g2;
    } else {
      if (z >= 0) {
        return z - g1;
      } else {
        return z + g1;
      }
    }
  } else {
    if ((z >= -g2) && (z <= g1)) {
      return 0;
    } else if ((z >= (-g1 + b)) && (z <= (-g2 + b))) {
      return b;
    } else if ((z >= (-g2 + b)) && (z <= -g2)) {
      return z + g2;
    } else {
      if (z >= 0) {
        return z - g1;
      } else {
        return z + g1;
      }
    }
  }
}

//' Optimize a linear regression model by coordinate descent algorithm using a covariance matrix
//'
//' @param Gamma covariance matrix of explanatory variables
//' @param gamma covariance vector of explanatory and objective variables
//' @param lambda lambda sequence
//' @param alpha alpha
//' @param init_beta initial values of beta
//' @param beta_tilde initial estimates (source parameters) of beta
//' @param maxit max iteration
//' @param eps convergence threshold for optimization
//' @param strong whether use strong screening or not
//' @return standardized beta
//' @keywords internal
//' @export
// [[Rcpp::export]]
NumericMatrix covCdaC(NumericMatrix Gamma_, NumericVector gamma_, 
                      NumericVector lambda_, 
                      NumericMatrix init_beta, NumericVector beta_tilde_, 
                      NumericVector penalty_factor_lasso_,
                      NumericVector penalty_factor_trlasso_,
                      double alpha = 1,
                      double maxit = 1e+4, double eps = 1e-04,
                      bool strong = true) {
  // initialize
  int p = Gamma_.ncol();
  int nlambda = lambda_.size();
  double z, b, prod;
  double *beta = R_Calloc(p * nlambda, double);
  double *g1 = R_Calloc(p, double);
  double *g2 = R_Calloc(p, double);
  double *beta_prev = R_Calloc(p, double);
  bool *candid = R_Calloc(p, bool);
  bool *candid2 = R_Calloc(p, bool);
  double *Gamma_beta = R_Calloc(p, double);
  for (int j=0; j<p; j++) for (int k=0; k<nlambda; k++) beta[j + k * p] = 0;
  for (int j=0; j<p; j++) beta_prev[j] = 0;
  for (int j=0; j<p; j++) candid[j] = false;
  for (int j=0; j<p; j++) candid2[j] = false;
  for (int j=0; j<p; j++) Gamma_beta[j] = 0;
  double *Gamma = REAL(Gamma_);
  double *gamma = REAL(gamma_);
  double *beta_tilde = REAL(beta_tilde_);
  double *lambda = REAL(lambda_);
  double *penalty_factor_lasso = REAL(penalty_factor_lasso_);
  double *penalty_factor_trlasso = REAL(penalty_factor_trlasso_);
  
  
  // CDA for all lambda
  for (int k = 0; k < nlambda; k++) {
    int kk = p * k;
    int kkk = p * (k-1);
    for (int j=0; j < p; j++) g1[j] = lambda[k] * (alpha * (penalty_factor_lasso[j] - penalty_factor_trlasso[j]) + penalty_factor_trlasso[j]);
    for (int j=0; j < p; j++) g2[j] = lambda[k] * (alpha * (penalty_factor_lasso[j] + penalty_factor_trlasso[j]) - penalty_factor_trlasso[j]);

    // warm-start
    if (k != 0) {
      for (int j=0; j < p; j++) {
        beta[j + kk] = beta[j + kkk];
      }
    }
    
    // strong screening
    if (strong && (k != 0)) {
      for (int j = 0; j < p; j++) {
        int jj = j * p;
        if (!candid2[j]) {
          prod = 0;
          for (int l = 0; l < p; l++) {
            prod += Gamma[l + jj] * beta[j + kk];
          }
          candid2[j] = (fabs(gamma[j] - prod) >= 2 * lambda[k] - lambda[k-1]);
        }
      }
    } else if (strong) {
      double lmax = 0;
      for (int j = 0; j < p; j++) {
        if (lmax < fabs(gamma[j])) {
          lmax = fabs(gamma[j]);
        }
      }
      for (int j = 0; j < p; j++) {
        prod = 0;
        int jj = j * p;
        for (int l = 0; l < p; l++) {
          prod += Gamma[l + jj] * beta[j + kk];
        }
        candid2[j] = (fabs(gamma[j] - prod) >= 2 * lambda[k] - lmax);
      }
    } else {
      for (int j = 0; j < p; j++) {
        candid2[j] = true;
      }
    }
    
    // CDA
    int i = 0;
    int violations = 0;
    double resid = 0;
    while (i < maxit) {
      while (i < maxit) {
        while (i < maxit) {
          i++;
          
          // first loop
          resid = 0;
          for (int j = 0; j < p; j++) {
            if (candid[j]) {
              prod = 0;
              int jj = j * p;
              for (int l = 0; l < p; l++) {
                if ((beta[l + kk] != 0) && (l != j)) {
                  prod += Gamma[l + jj] * beta[l + kk];
                }
              }
              z = gamma[j] - prod;
              if (z == 0) {
                beta[j + kk] = 0;
              } else {
                b = beta_tilde[j];
                beta[j + kk] = softThresholdC2(z, g1[j], g2[j], b) /  Gamma[j + jj];
              }

              double shift = beta[j + kk] - beta_prev[j];
              if (shift != 0) {
                beta_prev[j] = beta[j + kk];
                if (fabs(shift) > resid) {
                  resid = fabs(shift);
                }
              }
            }
          }
          if (resid < eps) {
            break;
          }
        }
        
        // second loop
        violations = 0;
        for (int j = 0; j < p; j++) {
          if ((!candid[j]) && candid2[j]) {
            prod = 0;
            int jj = j * p;
            for (int l = 0; l < p; l++) {
              if ((beta[l + kk] != 0) && (l != j)) {
                prod += Gamma[l + jj] * beta[l + kk];
              }
            }
            z = gamma[j] - prod;
            if (z == 0) {
              beta[j + kk] = 0;
            } else {
              b = beta_tilde[j];
              beta[j + kk] = softThresholdC2(z, g1[j], g2[j], b) /  Gamma[j + jj];
            }

            if (beta[j + kk] != 0) {
              if (!candid[j]) {
                violations++;
                candid[j] = true;
              } else {
              }
              
              double shift = beta[j + kk] - beta_prev[j];
              if (shift != 0) {
                beta_prev[j] = beta[j + kk];
              }
            } else {
            }
          }
        }
        if (violations == 0) {
          break;
        }
      }
      
      // third loop
      violations = 0;
      for (int j = 0; j < p; j++) {
        if (!candid2[j]) {
          prod = 0;
          int jj = j * p;
          for (int l = 0; l < p; l++) {
            if ((beta[l + kk] != 0) && (l != j)) {
              prod += Gamma[l + jj] * beta[l + kk];
            }
          }
          z = gamma[j] - prod;
          if (z == 0) {
            beta[j + kk] = 0;
          } else {
            b = beta_tilde[j];
            beta[j + kk] = softThresholdC2(z, g1[j], g2[j], b) /  Gamma[j + jj];
          }

          if (beta[j + kk] != 0) {
            candid[j] = candid2[j] = true;
            violations++;
            
            double shift = beta[j + kk] - beta_prev[j];
            if (shift != 0) {
              beta_prev[j] = beta[j + kk];
            }
          } else {
          }
        }
      }
      if (violations == 0) {
        break;
      }
    }
    
  }
  
  R_Free(beta_prev); R_Free(candid); R_Free(candid2); R_Free(Gamma_beta);
  
  SEXP result = PROTECT(Rf_allocMatrix(REALSXP, p, nlambda));
  for (int k=0; k<nlambda; k++) {
    int kk = p * k;
    for (int j=0; j<p; j++){
      SET_REAL_ELT(result, j + kk, beta[j + kk]);
    }
  }
  UNPROTECT(1);
  
  return result;
}
