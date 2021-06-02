
#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
typedef std::vector<double> stdvec;

//' Rcpp function to get the gene profile of ambient droplets 
//' @param m a sparse matrix of ambient drpolets
//' @param nrow the rows of m
//' @param ncol the cols of m
// [[Rcpp::export]]
NumericVector getAmbientProp(arma::sp_mat m,int nrow,int ncol) {
 
  mat ones(1, nrow, fill::ones);
  
  mat  bcsum = ones*m;
  mat  ref(1,ncol);
  
  ref.fill(100.0);
  mat maxUMI = min(bcsum,ref);
  arma::mat adds(1,ncol);
  adds.fill(1.0);
  mat adds_ = maxUMI+adds;
  mat logAdd = log(adds_);
  mat gene(nrow,1);
  int i=0;
  mat calsum(nrow,1,fill::zeros);
  for (i =0; i<ncol; ++i){
    double rem = i%10000;
    if (rem==0){cout<<i<<endl;}
    mat prop_i=mat(m.col(i))/bcsum(i)*logAdd(i);
    calsum = prop_i + calsum;
  }
  gene = calsum/ncol;
  vec out=gene.col(0);
  return (Rcpp::NumericVector(out.begin(),out.end()));
  
}



