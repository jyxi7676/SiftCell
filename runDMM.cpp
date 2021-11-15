#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <ctime>
#include <cstdlib>
#include <time.h>
#include <set>
#include <Rcpp.h>
//#include <bits/stdc++.h>
using namespace Rcpp;
using namespace std;



//' Rcpp function to get the coefficieint of DMM
//' @param geneProfile matrix of gene profile for all cell types and soup
//' @param mtx a DGE matrix
//' @param ub_p a vector of upper bound for p
//' @param lb_p a vector of lower bound for p
//' @return a matrix with convergence and coefficients
// [[Rcpp::export]]
NumericMatrix runDMM(NumericMatrix geneProfile, NumericMatrix mtx , NumericVector ub_p, NumericVector lb_p)
{
  int d = geneProfile.nrow();
  int n = mtx.ncol();
  NumericVector pu = runif(d,0,1);
  NumericVector p0(d);
  NumericMatrix dmm(d+1,n);

  for(int i=0;i<d;i++)
  {
    p0[i] = pu[i]/sum(pu);
    //std::cout<<p0[i]<<std::endl;
  }

  for ( int i=0;i<n;i++)
  {
    //std::cout << i <<std::endl;
    NumericVector umi = mtx(_,i);
    Function f("solveDMM");
    NumericVector dmmi = f(Named("p0")=p0, Named("umi") = umi, Named("geneProfile") = geneProfile, Named("lb_p") = lb_p,Named("ub_p")= ub_p);
    //print(dmmi);
    dmm(_,i) = dmmi;
  }
      //solnp(p0,Named("fun") = obj,Named("eqfun") = cons,Named("eqB")=1, Named("LB") = lb_p,Named("UB")= ub_p);
//NumericMatrix dmmi = runSolnp();
    //R::print(dmmi);
//(p0, fun = function(x) {objectiveFn(x,umi,geneProfile)}, eqfun = constraintFn,eqB = 1,LB = rep(lb_p,d),UB = rep(ub_p,d),control = list(trace = 0));return(c(dmm$pars,dmm$convergence))} )

  //  NumericMatrix dmm = f(p0, fun = function(x) {objectiveFn(x,umi,geneProfile)}, eqfun = constraintFn,eqB = 1,LB = rep(lb_p,d),UB = rep(ub_p,d),control = list(trace = 0));return(c(dmm$pars,dmm$convergence))} )



return (dmm);
}
