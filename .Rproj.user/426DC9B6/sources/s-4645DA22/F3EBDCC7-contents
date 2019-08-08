#' Conversion of non-zero gene FPKM value into p-value using each cells mean and standard deviation
#'
#' @param x gene expression matrix
#'
#' @return P-value matrix
#' @export
#'
#' @examples
#' binorm()
binorm <- function(x)
{
  x = log2(100*x + 1) ;
xdim = dim(x) ;
xpval = matrix(1 , xdim[1] , xdim[2],dimnames=list(toupper(rownames(x)),colnames(x)))
for (i in 1: xdim[2])
{
pos = which(x[,i] > 0.00000000001) ;
dpos = dim(as.matrix(pos)) ;
pnz = dpos[1]/ xdim[1] ;
meannz = mean(x[pos, i]) ;
sdnz = sd(x[pos,i]) ;
xpval[pos, i] =  (1 - pnorm(x[pos, i], meannz, sdnz)) ;
}
return(xpval) ;
}
