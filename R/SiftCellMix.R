
#'This function get the gene profile
#' @param mtx DGE matrix
#' @param mtx_threshold DGE matrix with colsums > threshold
#' @param celltype dataframe with celltype information, 1st coloum is barcode, 2nd column is celltype
#' @param alpha scale between 0 and 1 used to adjust 0 value in gene profile for droplets with celltype infomation
#' @param threshold umi cut off
#' @return a matrix with gene profile for all the celltypes and soup
#' @export
getGeneProfile = function(mtx,mtx_threshold,celltype,alpha=0.01,threshold)
{
  cell.type = unique(celltype$TYPE)
  cell.type = c(cell.type[order(cell.type)],"AMB")
  n1 = dim(mtx)[1]
  n2 = dim(mtx)[2]
  bcsum=colSums(mtx)
  #weight=pmin(100,bcsum)+1
  #logweight=log(weight)
  #denom=sum(log(weight))
  #comp1=sweep(mtx,2,bcsum,'/')
  #comp2=sweep(comp1,1,logweight,'*')
  #sweep(matrix(1:9,3,3), 2, c(1,2,3), "/")
  #nom = rowSums(mtx / matrix(bcsum*logweight,nrow=n1,ncol=n2,byrow = T))
  #nom = rowSums(mtx / (bcsum*logweight))
  #amb=nom/denom
  #amb=getAmbientProp(mtx,n1,n2)
  #p.AMB=amb/sum(amb)

  p.AMB = getAmbientProp(mtx,n1,n2)
  p.AMB = p.AMB/sum(p.AMB)
  d = length(cell.type)
  p = sapply((1:(d-1)),function(x) {m = mtx_threshold[,match(celltype[celltype$TYPE==cell.type[x],]$BARCODE,colnames(mtx_threshold))];return((1-alpha)*rowSums(m)/sum(m)+alpha*p.AMB) })

  remove(mtx)
  remove(mtx_threshold)
  p = cbind(p,p.AMB)
  colnames(p) = paste0("p.",cell.type)
  return(p)
}

# getGeneProfile = function(mtx,mtx_threshold,celltype,alpha=0.01,threshold)
# {
#  cell.type = unique(celltype$TYPE)
#  cell.type = c(cell.type[order(cell.type)],"AMB")
#  n1 = dim(mtx)[1]
#  n2 = dim(mtx)[2]
#  p.AMB = getAmbientProp(mtx,n1,n2)
#  p.AMB = p.AMB/sum(p.AMB)
#  d = length(cell.type)
#  p = sapply((1:(d-1)),function(x) {m = mtx_threshold[,match(celltype[celltype$TYPE==cell.type[x],]$BARCODE,colnames(mtx_threshold))];return((1-alpha)*rowSums(m)/sum(m)+alpha*p.AMB) })
#
#  remove(mtx)
#  remove(mtx_threshold)
#  p = cbind(p,p.AMB)
#  colnames(p) = paste0("p.",cell.type)
#  return(p)
# }




#' This function defines the objective function
#' @param p proportion coeefficient of Multinomial Mixture (MM) model
#' @param umi a vector of UMIs
#' @param geneProfile this is a transposed geneProfile(as a matrix)
#' @return negative log likelihood of MM
objectiveFn= function(p,umi,geneProfile)
{
  logGP=log(geneProfile)
  umi_mt=matrix(umi,nrow=length(p),ncol=length(umi),byrow = T)
  product=logGP*umi_mt
  avg=mean(rowSums(product))

  logp=log(p)
  sum_genes=rowSums(product)-avg+logp
  #sum_genes=rowSums(product)+logp
  logmm=-logsum(sum_genes)
  #multinomial=sapply(1:length(p), function(i) { geneProfile_k=geneProfile[i,]; m_k= dmultinom(umi, prob=geneProfile_k,log=TRUE)})
  #logmm=-logsum(p*multinomial)
  return (logmm)
}
#' This defines the equality constraint for the vector of proportion coefficients
#' @param p proportion coeefficiento of DMM
#' @return summation
constraintFn =function(p)
{
  return(sum(p))
}





#' This function get the proportion coefficient for Multinomial Mixture Model
#' @param geneProfile matrix of gene profile for all cell types and soup
#' @param m DGE matrix
#' @param seed seed for generating starting point
#' @param ub_p a scalar of upper bound of p
#' @param lb_p a scalar of lower bound of p
#' @return matrix of estimated proportion of MM
#' @export
getP = function(geneProfile,m,seed = 0, ub_p=1, lb_p=0)
{

 #m = as.matrix(m)
 d = nrow(geneProfile)
 lb_p=rep(lb_p,d)
 ub_p=rep(ub_p,d)
 umi=colSums(m)
 set.seed(seed)
 #p0=rep(1/d,d)


 print('Start Parralel Computing')
 #parralel computing the fraction
 num_cores = detectCores()
 # Register the parallel backend
 cl = makeCluster(num_cores)
 registerDoParallel(cl)



# clusterExport(cl, c("solnp", "logsum","objectiveFn", "geneProfile","constraintFn", "lb_p", "ub_p","d"),envir = .GlobalEnv)
# clusterExport(cl, c("solnp", "logsum","objectiveFn","constraintFn"),envir = .GlobalEnv)
 #clusterExport(cl, c("logsum","objectiveFn","constraintFn"),envir = .GlobalEnv)

 # Define the function to be applied in parallel
 parallel_function <- function(umi) {

    out = solnp(rep(1/d,d), fun = function(x) {objectiveFn(x, umi, geneProfile)},
                 eqfun = constraintFn, eqB = 1,
                 LB = lb_p, UB = ub_p,
                 control = list(trace = 0))
   return(cbind(t(out$pars),out$convergence))
 }

 result = parApply(cl,m, 2, parallel_function)

 # Stop the parallel backend
 stopCluster(cl)
 registerDoSEQ()

 result=t(result)
 colnames(result) = c(rownames(geneProfile),"convergence")
 result=as.data.frame(result)
 remove(m)

 if(sum(result$convergence)==0)
 {
   print("Finished with convergence = 0")
 }
 else
 {
   stop ("No convergence!")
 }
 return(result[,1:d])
 }

#' This function runs SiftCellMix
#' @param workingdir directory with DGE files
#' @param celltype dataframe with two columns, 1st column is cell barcode, 2nd column is cell type
#' @param alpha scalar between 0 and 1
#' @param threshold focus on droplets with UMI>threshold
#' @param ub_p upper bound of proportion. Default is 1
#' @param lb_p lower bound of proportion. Default is 0
#' @return estimated proportion
#' @export
SiftCellMix = function(workingdir,celltype,alpha=0.01,threshold=100,ub_p=1, lb_p=0)
{
  if (dir.exists(workingdir))
  {setwd(workingdir)
  } else
  {stop ("Working directory does not exist")
  }
  colnames(celltype) = c("BARCODE","TYPE")
  print("Read and preprocess DGE")
  #tic()
  orgDGE = readDGE(paste0(workingdir,"/DGE/"))
  if(any(is.na(match(celltype$BARCODE,colnames(orgDGE)))))
  {
    stop ("The Barcodes in celltype cannot be be found in orginal DGE!")
  }
  orgDGE = orgDGE [,colSums(orgDGE)>0]
  orgDGE = orgDGE [rowSums(orgDGE)>0,]
  orgDGE_threshold=orgDGE[,colSums(orgDGE)>threshold]
  #toc()

  print('Get gene profile')
  #tic()
  geneProfile = t(as.matrix(getGeneProfile(orgDGE,orgDGE_threshold, celltype,alpha,threshold)))
  remove(orgDGE)
  #toc()

  print('Estimatae fraction')
  tic()
  prop = getP(geneProfile,orgDGE_threshold, seed = 0, ub_p, lb_p)
  toc()
  remove(orgDGE_threshold)
  write.csv(prop,paste0(workingdir,'/fraction.csv'))
}
