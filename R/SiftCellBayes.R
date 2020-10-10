
#'This function get the gene profile
#' @param mtx DGE matrix
#' @param celltype dataframe with celltype information, 1st coloum is barcode, 2nd column is celltype
#' @param alpha scale between 0 and 1 used to adjust 0 value in gene profile for droplets with celltype infomation
#' @return a matrix with gene profile for all the celltypes and soup
#' @export
getGeneProfile = function(mtx, celltype,alpha=0.01)
{
 colnames(celltype) = c("BARCODE","TYPE")
 if(any(is.na(match(celltype$BARCODE,colnames(mtx)))))
 {
   stop ("The Barcodes in celltype cannot be be found in orginal DGE!")
 }
 cell.type = unique(celltype$TYPE)
 cell.type = c(cell.type[order(cell.type)],"AMB")
 p.AMB =rowSums(mtx) / sum(mtx)
 d = length(cell.type)
 p = sapply((1:(d-1)),function(x) {m = mtx[,match(celltype[celltype$TYPE==cell.type[x],]$BARCODE,colnames(mtx))];return((1-alpha)*rowSums(m)/sum(m)+alpha*p.AMB) })
 p = cbind(p,p.AMB)
 colnames(p) = paste0("p.",cell.type)
 return(p)
}

# obj =function(p)
# {
#   return(-sum(umi*log(p%*%geneProfile)))
#
# }
#
# cons = function(p)
# {
#   return(sum(p))
# }


#' This function defines the objective function
#' @param p proportion coeefficiento of DMM
#' @param umi a vector of UMIs
#' @param geneProfile this is a transposed geneProfile(as a matrix)
#' @return negative log likelihood of DMM
objectiveFn = function(p,umi,geneProfile)
{
  return(-sum(umi*log(p%*%geneProfile)))
}
#' This defines the equality constraint for the vector of proportion coefficients
#' @param p proportion coeefficiento of DMM
#' @return summation
constraintFn =function(p)
{
  return(sum(p))
}

#' The function solves the p by solnp in Rsolnp package
#' @param p0 The startin points
#' @param umi vector of UMIs in a droplet
#' @param geneProfile matrix of gene profile for all cell types and soup
#' @param lb_p vector of lower bound
#' @param ub_p vecetor of upper bound
#' @return vector of estimated proportion and convergence value
#' @import Rsolnp
solveDMM = function(p0,umi,geneProfile,lb_p,ub_p)
{

  dmm=solnp(p0, fun = function(x) {objectiveFn(x,umi,geneProfile)}, eqfun = constraintFn,eqB = 1,LB = lb_p,UB = ub_p,control = list(trace = 0))
  return(c(dmm$pars,dmm$convergence))
}



#' This function get the proportion coefficient for Multinomial Mixture Model
#' @param geneProfile matrix of gene profile for all cell types and soup
#' @param mtx DGE matrix
#' @param threshold focus on droplets with UMI>thershold
#' @param seed seed for generating starting point
#' @param ub_p a scalar of upper bound of p
#' @param lb_p a scalar of lower bound of p
#' @return matrix of estimated proportion of DMM
getP = function(geneProfile,mtx,threshold=100,seed = 0, ub_p=1, lb_p=0)
{
 #gp.t = t(as.matrix(geneProfile))
 m = as.matrix(mtx[,colSums(mtx)>threshold])
 d = nrow(geneProfile)
 #p0 = (p0=runif(d,0,1))/sum(p0)
 #umi = m[,i]
 set.seed(seed)
 result = runDMM(geneProfile,m,rep(ub_p,d),rep(lb_p,d)) #speed up using Rcpp
 #result = sapply(1:ncol(m),function(i) {umi = m[,i];dmm=solnp(p0, fun = function(x) {objectiveFn(x,umi,geneProfile)}, eqfun = constraintFn,eqB = 1,LB = rep(lb_p,d),UB = rep(ub_p,d),control = list(trace = 0));return(c(dmm$pars,dmm$convergence))} )
 rownames(result) = c(rownames(geneProfile),"convergence")
 colnames(result) = colnames(m)
 #result = sapply(1:ncol(m),function(i) {umi = m[,i];dmm=solnp(p0, fun = function(x) {objectiveFn(x,umi,geneProfile)}, eqfun = constraintFn,eqB = 1,LB = rep(lb_p,d),UB = rep(ub_p,d),control = list(trace = 0));conv=c(conv,dmm$convergence);return(dmm$pars)} )
 if(sum(result["convergence",])==0)
 {
   print("Finished with convergence = 0")
 }
 else
 {
   stop ("No convergence!")
 }
 return(result[1:d,])
 #solnp(p0, fun = function(x) {objectiveFn(x,umi,gp.t)}, eqfun = constraintFn,eqB = 1,LB = rep(0,d),UB = rep(1,d))
 }

#' This function runs SiftCellBayes
#' @param workingdir directory with DGE files
#' @param celltype dataframe with two columns, 1st column is cell barcode, 2nd column is cell type
#' @param alpha scalar between 0 and 1
#' @param threshold focus on droplets with UMI>threshold
#' @param ub_p upper bound of proportion. Default is 1
#' @param lb_p lower bound of proportion. Default is 0
#' @return estimated proportion
#' @export
SiftCellBayes = function(workingdir,celltype,alpha=0.01,threshold=100,ub_p=1, lb_p=0)
{
  if (dir.exists(workingdir))
  {setwd(workingdir)
  } else
  {stop ("Working directory does not exist")
  }
  orgDGE = readDGE(paste0(workingdir,"/DGE/"))
  orgDGE = orgDGE [rowSums(orgDGE) > 0,]
  orgDGE = orgDGE[,colSums(orgDGE)>0]
  geneProfile = t(as.matrix(getGeneProfile(orgDGE, celltype,alpha)))
  prop = getP(geneProfile,orgDGE,threshold=100,seed = 0, ub_p, lb_p)
  return(t(prop))
}
