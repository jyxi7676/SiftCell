#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <ctime>
#include <cstdlib>
#include <time.h>
#include <set>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;


/////////////////////////////////////////////////////////////////////////
// dge-shuffle : Shuffle digital expression matrix
// Author: Hyunmin Kang&Jingyue Xi
////////////////////////////////////////////////////////////////////////
//' This function generate shuffled DGE
//' @param n_genes number of genes
//' @param n_bcds number of barcodes
//' @param num_nonzeros number of non-zero umis
//' @param numis vector of umis
//' @param geneInd vector of gene index
//' @param barcodeInd index of barcoders
//' @param totalumi sum of umis
//' @return a list
//' @export
// [[Rcpp::export]]
List dgeShuffle(int n_genes, int n_bcds, int num_nonzeros, NumericVector numis,NumericVector geneInd, NumericVector barcodeInd,int totalumi)
{
  int iumi=0;
  //int seed = 0;
  //int sum_umis = num_nonzeros;
  int N = numis.size();
  std::vector<int> igenes(totalumi);
  std::vector<int> ibcds(totalumi);
  for(int i=0;i<N;i++)
  {
    int numi = numis[i];

    for(int j=0; j < numi; ++j)
      {
        igenes[iumi]=geneInd[i];
        ibcds[iumi] = barcodeInd[i];
        ++iumi;
      }

  }


  // set the random seed
  //if ( seed == 0 ) srand(time(0));
  //else srand((uint32_t) seed);

  //notice("Randomizing UMIs with seed %d", seed);
  for(int i=0; i < totalumi-1; ++i) {
    int r = i + round(R::runif(0,totalumi-i-1));
    //int r = i + ( rand() % (totalumi-i) );
    if ( i != r ) { // swap
      int tmp = igenes[i];
      igenes[i] = igenes[r];
      igenes[r] = tmp;
    }
  }

  // List ret;
  // ret["igenes"] = igenes;
  // ret["ibcds"] = ibcds;
  // return(ret);

std::map< int, std::map<int, int> > dge;
int new_num_nonzeros = 0;
for(int i=0; i < totalumi; ++i)
 {
  if ( ++dge[ibcds[i]][igenes[i]] == 1 )
    ++new_num_nonzeros;
 }


std::vector<int> bcd2;
std::vector<int> gene2;
std::vector<int> umi2;
// For accessing outer map
map<int, map<int, int> >::iterator itr;

// For accessing inner map
map<int, int>::iterator ptr;

for (itr = dge.begin(); itr != dge.end(); itr++) {

  for (ptr = itr->second.begin(); ptr != itr->second.end(); ptr++) {

    bcd2.push_back (itr->first);
    gene2.push_back (ptr->first);
    umi2.push_back (ptr->second);
  }
}

List ret;
ret["igenes"] = gene2;
ret["ibcds"] = bcd2;
ret["umi"] = umi2;

return(ret);
}


