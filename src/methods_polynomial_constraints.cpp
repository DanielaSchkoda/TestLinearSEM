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

// Estimate Polynomials 
// The R function esimate_poly is translated to C++ as a member function of the C++ object
// "Polynomial". For defining the class Polynomial, the class "Monomial" will be used.
class Monomial {
  public:
    double coef;
    // moms has 3 cols, as many rows as moms occur
    IntegerMatrix moms;
    int degree;

    // constructor
    Monomial(double coef_in, IntegerMatrix moms_in){
      coef = coef_in;
      moms = moms_in;
      degree = moms_in.nrow();
    }; 
};

int get_degree(Monomial mon){
  return(mon.degree);
}

class Polynomial {
  public:
    int degree;
    vector<Monomial> monomials;

    // constructor
    Polynomial(vector<Monomial> monomials_in){
      monomials = monomials_in;
      vector<int> degrees_mons;
      transform(monomials.begin(), monomials.end(), back_inserter(degrees_mons), get_degree);
      degree = *max_element(degrees_mons.begin(), degrees_mons.end());
    }; 

    double estimate_poly(NumericMatrix L){
      double poly_est = 0;
      for(int i=0; i < monomials.size(); i++){ 
        Monomial current_mon = monomials[i];
        double monomial_est = current_mon.coef;
        // -1 since array indexing starts with 0 in C++
        IntegerMatrix moms = current_mon.moms - 1;
        int m_deg = 0;
        for(int j=0; j < moms.nrow(); j++){
          // If last entry is -1, the row represents a second moment
          if (moms(j, 2) == -1){
            monomial_est = monomial_est * L(m_deg, moms(j, 0)) * L(m_deg, moms(j, 1));
          } else{
            monomial_est = monomial_est * L(m_deg, moms(j, 0)) * L(m_deg, moms(j, 1)) * L(m_deg, moms(j, 2));
          }
          m_deg++;
        }
        poly_est += monomial_est;
    }
    return(poly_est);
  };
};

// Function to create C++ Polynomial objects from a list that represenets a polynomial. Function is needed
// as the C++ functions should be invoked from R code, where the C++ object "Polynomial" does
// not exist. The argument R_poly needs to be a list of monomials where each monomial is 
// a vector of the form c(coef=1, moms=matrix(c(1, 2, 1), 1, 3))
// moms needs to have k columns. Each row represents one
// second respectively third moment appearing in the monomial. For second moments, the *last* entry in the row needs to be 0.
Polynomial create_c_polynomial(List R_poly){
  vector<Monomial> C_monomials;

  for (int j = 0; j < R_poly.length(); j++){
    List R_mon = as<List>(R_poly[j]);
    Monomial C_mon();
    double coef = as<double>(R_mon["coef"]);

    IntegerMatrix moms = as<IntegerMatrix>(R_mon["moms"]);
    C_monomials.push_back(Monomial(coef, moms));
  }
  return(C_monomials);
}

// [[Rcpp::export]]
double estimate_polynomial_c(List R_polynomial, NumericMatrix L){
  Polynomial C_poly = create_c_polynomial(R_polynomial);
  return(C_poly.estimate_poly(L));
}

// Calculate H and G 
NumericVector h_c(NumericMatrix L, vector<Polynomial> C_polynomials){
  NumericVector h(C_polynomials.size());
  IntegerMatrix perm = permutations(L.nrow());
  perm = perm - 1;
  
  for (int per = 0; per < perm.nrow(); per++){
    for (int poly = 0; poly < C_polynomials.size(); poly++){
      Polynomial curr_poly = C_polynomials[poly];
      NumericMatrix L_permuted(L.nrow(), L.ncol());
      for (int j = 0; j < L.nrow(); j++) {
        L_permuted(j, _) = L(perm(per, j), _);
      }
      h[poly] = h[poly] + curr_poly.estimate_poly(L_permuted);
    }
  }
  return(h/perm.nrow());
}

// [[Rcpp::export]]
Rcpp::NumericMatrix H_c(NumericMatrix X, IntegerMatrix indices, List R_polynomials){
  indices = indices - 1;
  NumericMatrix H(indices.nrow(), R_polynomials.length());
  vector<Polynomial> C_polynomials;
  for (int i = 0; i < R_polynomials.length(); i++){
    List R_poly = as<List>(R_polynomials[i]);
    C_polynomials.push_back(create_c_polynomial(R_poly));
  }
  for (int i = 0; i < indices.nrow(); i++){
    NumericMatrix L(indices.ncol(), X.ncol());
    for (int j=0; j<indices.ncol();j++) {
      L(j,_) = X(indices(i, j), _);
    }
    H(i,_) = h_c(L, C_polynomials);
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

NumericVector g_c(NumericMatrix X, int i, int order_kernel, vector<Polynomial> C_polynomials){
  int n = X.nrow();
  NumericVector g(C_polynomials.size());
  int K = floor((n-1)/(order_kernel-1));
  
  IntegerMatrix S = compute_S(n,i,order_kernel-1);
  S = S - 1;
  
  for (int k = 0; k < K; k++){
    NumericMatrix L(order_kernel, X.ncol());
    L(0, _) = X((i-1),_);
    for (int j = 0; j < order_kernel-1; j++){
      L(j+1, _) = X(S(k,j),_);
    }
    g = g + h_c(L, C_polynomials);
  }
  return(g/K); 
}

// [[Rcpp::export]]
NumericMatrix G_c(NumericMatrix X, int order_kernel, List R_polynomials){
  int n = X.nrow();
  NumericMatrix G(n, R_polynomials.length());
  vector<Polynomial> C_polynomials;
  for (int i = 0; i < R_polynomials.length(); i++){
    List R_poly = as<List>(R_polynomials[i]);
    C_polynomials.push_back(create_c_polynomial(R_poly));
  }
  
  for (int i = 1; i <= n; i++){
    G((i-1),_) = g_c(X, i, order_kernel, C_polynomials);
  }
  return(G);
}
