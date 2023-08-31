// The code is adapted from Nils Sturma, https://github.com/NilsSturma/TestGGM
#include <Rcpp.h>
#include <set>
#include <array>
#include <vector>

using namespace std;
using namespace Rcpp;

// Helpers
// [[Rcpp::export]]
Rcpp::List random_combs_c(int n, int k, int nr){
  set<set<int>> sub_sets;
  
  Rcpp::IntegerVector v(n);
  
  for(int i=0; i < n; i++){
    v[i] = i+1;
  }
  
  set<int> s;
  Rcpp::IntegerVector w(k);
  
  while(sub_sets.size() < (unsigned)nr){
    w = sample(v,k);
    for (int i = 0; i < k; i++) {
      s.insert(w[i]);
    }
    sub_sets.insert(s);
    s.clear();
  }
  
  Rcpp::List res = Rcpp::List::create(sub_sets);
  
  return(res);
}

// [[Rcpp::export]]
IntegerMatrix permutations(int n){
  
  int factorial = 1;
  for (int i=1; i <= n; i++){
    factorial = factorial * i;
  }
  
  IntegerVector v(n);
  for (int i=1; i <=n; i++){
    v[i-1] = i;
  }
  
  IntegerMatrix res(factorial,n);
  for (int i=0; i<factorial; i++){
    res(i,_) = v; 
    std::next_permutation(v.begin(), v.end());
  }
  return res;
}

// [[Rcpp::export]]
double estimate_polynomial(List poly, NumericMatrix L){
  double poly_est = 0;
  for(int i=0; i < poly.size(); i++){ 
    List current_mon = as<List>(poly[i]);
    double monomial_est = current_mon["coeff"];
    List moms = current_mon["moms"];
    int m_deg = 0;
    for(int j=0; j < moms.size(); j++){
      NumericVector curr_mom = moms[j];
      //std::cout << curr_mom << ' ';
      // start with l=1 since first char is m
      for(int l=0; l < curr_mom.size(); l++) {
        // -1 since array indexing starts with 0 in C++
        monomial_est = monomial_est * L(m_deg, curr_mom[l]-1);
      }
      m_deg++;
    }
    poly_est += monomial_est;
  }
  return(poly_est);
}

// Calculate H and G 
NumericVector h_c(NumericMatrix L, List polynomials){
  NumericVector h(polynomials.size());
  IntegerMatrix perm = permutations(L.nrow());
  perm = perm - 1;
  
  for (int per = 0; per < perm.nrow(); per++){
    for (int poly = 0; poly < polynomials.size(); poly++){
      List curr_poly = polynomials[poly];
      NumericMatrix L_permuted(L.nrow(), L.ncol());
      for (int j = 0; j < L.nrow(); j++) {
        L_permuted(j, _) = L(perm(per, j), _);
      }
      h[poly] = h[poly] + estimate_polynomial(curr_poly, L_permuted);
    }
  }
  return(h/perm.nrow());
}

// [[Rcpp::export]]
Rcpp::NumericMatrix H_c(NumericMatrix X, IntegerMatrix indices, List polynomials){
  indices = indices - 1;
  NumericMatrix H(indices.nrow(), polynomials.length());
  for (int i = 0; i < indices.nrow(); i++){
    NumericMatrix L(indices.ncol(), X.ncol());
    for (int j=0; j<indices.ncol();j++) {
      L(j,_) = X(indices(i, j), _);
    }
    H(i,_) = h_c(L, polynomials);
  }
  
  return(H);
}

IntegerMatrix compute_S(int n, int i, int L){
  int K = floor((n-1)/L);
  IntegerVector v(n-1);
  for (int j = 1; j < (n+1); ++j){
    if (j < i){
      v[j-1] = j;
    } else if (j==i) {
      continue;
    } else {
      v[j-2] = j;
    }
  }
  IntegerVector sequence(L);
  IntegerMatrix res(K,L);
  for (int k = 0; k < K; ++k){
    sequence = seq((k*L), (((k+1)*L)-1));
    for (int l = 0; l < L; l++){
      res(k,l)=v[sequence[l]];
    }
  }
  return res;
}

NumericVector g_c(NumericMatrix X, int i, int order_kernel, List polynomials){
  int n = X.nrow();
  NumericVector g(polynomials.size());
  int K = floor((n-1)/(order_kernel-1));
  
  IntegerMatrix S = compute_S(n,i,order_kernel-1);
  S = S - 1;
  
  for (int k = 0; k < K; k++){
    NumericMatrix L(order_kernel, X.ncol());
    L(0, _) = X((i-1),_);
    for (int j = 0; j < order_kernel-1; j++){
      L(j+1, _) = X(S(k,j),_);
    }
    g = g + h_c(L, polynomials);
  }
  return(g/K); 
}

// [[Rcpp::export]]
NumericMatrix G_c(NumericMatrix X, int order_kernel, List polynomials){
  int n = X.nrow();
  NumericMatrix G(n, polynomials.length());
  for (int i = 1; i <= n; i++){
    G((i-1),_) = g_c(X, i, order_kernel, polynomials);
  }
  return(G);
}
