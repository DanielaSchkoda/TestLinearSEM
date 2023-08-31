#include <Rcpp.h>
#include <stdexcept>
#include <iostream>

using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

NumericVector concatenateVectors(NumericVector vec1, NumericVector vec2) {
  int size1 = vec1.size();
  int size2 = vec2.size();
  
  NumericVector result(size1 + size2);
  
  // Copy elements from vec1 to result
  for (int i = 0; i < size1; i++) {
    result[i] = vec1[i];
  }
  
  // Copy elements from vec2 to result
  for (int i = 0; i < size2; i++) {
    result[size1 + i] = vec2[i];
  }
  
  return result;
}


// Function to estimate a moment
// [[Rcpp::export]]
double estimateMoment(NumericVector moment, NumericMatrix X) {
  int numIndices = moment.size();
  if (numIndices == 0) {
    return 0.0;
  } else {
    int n = X.nrow();
    NumericVector product = Rcpp::NumericVector(n, 1.0);
    for (int k = 0; k < moment.size(); k++) {
      // -1 since C++ indices start from 0
      //Rcout << mean(X(_,moment[k] - 1));
      product = product * X(_,moment[k] - 1);
    };
    return(mean(product));
  }
}

// Function to estimate a cumulant. Works up to order four
// [[Rcpp::export]]
double estimateCumulant(NumericVector cumulant, NumericMatrix X) {
  int numIndices = cumulant.size();
  if (numIndices == 0) {
    return 0.0;
  // for orders 1, 2 and 3, cumulant coincides with moment
  } else if (numIndices <= 3) {
    return(estimateMoment(cumulant, X));
  } else if (numIndices == 4) {
    return(estimateMoment(cumulant, X) 
             - estimateMoment({cumulant[0], cumulant[1]}, X) * estimateMoment({cumulant[2], cumulant[3]}, X)
             - estimateMoment({cumulant[0], cumulant[2]}, X) * estimateMoment({cumulant[1], cumulant[3]}, X)
             - estimateMoment({cumulant[0], cumulant[3]}, X) * estimateMoment({cumulant[1], cumulant[2]}, X));
  } else {
    throw std::invalid_argument("Function estimateCumulant only implemented for cumulants with order <= 4");
  }
}

// Function to estimate covariance between two moments
// [[Rcpp::export]]
double estimateCovBtwMoments(List moment_1, List moment_2, NumericMatrix X) {
  NumericVector indices_1 = moment_1["mom"];
  NumericVector indices_2 = moment_2["mom"];
  double coeff_1 = moment_1["coeff"];
  double coeff_2 = moment_2["coeff"];
  int numIndices_1 = indices_1.size();
  int numIndices_2 = indices_2.size();
  
  if (numIndices_1 == 0 || numIndices_2 == 0) {
    return 0.0;
  } else {
    int n = X.nrow();
    NumericVector combined_indices = concatenateVectors(indices_1, indices_2);
    double moment_1_mean = estimateMoment(indices_1, X);
    double moment_2_mean = estimateMoment(indices_2, X);
    double moment_product_mean = estimateMoment(combined_indices, X);
    double result = coeff_1 * coeff_2 * n / (n - 1) * (moment_product_mean - moment_1_mean * moment_2_mean);
    return result;
  } 
}

// // Function to estimate the covariance of two cumulants
// // [[Rcpp::export]]
// double estimateCovBtwCumulants(List moment_1, List moment_2, NumericMatrix X) {
//   NumericVector indices_1 = moment_1["mom"];
//   NumericVector indices_2 = moment_2["mom"];
//   double coeff_1 = moment_1["coeff"];
//   double coeff_2 = moment_2["coeff"];
//   int numIndices_1 = indices_1.size();
//   int numIndices_2 = indices_2.size();
//   
//   if (numIndices_1 == 0 || numIndices_2 == 0) {
//     return 0.0;
//   } else if (numIndices_1 <= 3 && numIndices_2 <= 3)
//     return(estimateCovBtwMoments(moment_1, moment_2, X));
//   } else if (numIndices_1 <= 3 || numIndices_2 <= 3) {
//     return(estimateCovBtwMoments(moment_1, moment_2, X) - );
//   }
// }

// [[Rcpp::export]]
NumericMatrix estimate_Pi(Rcpp::List Pi, NumericMatrix X) {
  int numCols = Pi.size();
  Rcpp::List firstCol = as<List>(Pi[0]);
  int numRows = firstCol.size();

  NumericMatrix result(numRows, numCols);
      
  for (int j = 0; j < numCols; j++) {
    Rcpp::List col = as<List>(Pi[j]);
    for (int i = 0; i < numRows; i++) {
      Rcpp::List entry = as<List>(col[i]);
      double coefficient = Rcpp::as<double>(entry["coeff"]);
      Rcpp::NumericVector moment = Rcpp::as<Rcpp::NumericVector>(entry["mom"]);
      result(i, j) = coefficient * estimateMoment(moment, X);
    }
  }
  return(result);
}


// [[Rcpp::export]]
NumericMatrix estimate_W(List Pi_vectorized, NumericMatrix X) {
  int length = Pi_vectorized.size();

  NumericMatrix W(length, length);
  for (int i = 0; i < length; i++) {
    for (int j = 0; j < length; j++) {
      List moment_1 = as<List>(Pi_vectorized[i]);
      List moment_2 = as<List>(Pi_vectorized[j]);
      W(i, j) = estimateCovBtwMoments(moment_1, moment_2, X);
    }
  }
  
  return W;
}



















    
      



      
    
  


