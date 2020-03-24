#' Adjusting of combine p-values using null model
#'
#' @param combp combined p-value matrix obtained using gene expression matrix
#' @param combp_ref combined p-value matrix obtained using null model
#'
#' @return A list of 3 matrices, One is absolute p-value matrix, second is raw adjusted p-value matrix and third is log transformed adjusted p-value matrix
#' @export
#'
#' @examples
#' adjust()
adjust <- function(combp,combp_ref){

adjpva <- matrix(0,nrow(combp),ncol(combp),dimnames=list(rownames(combp),colnames(combp)))
for( i in 1:nrow(combp)){
  for(j in 1:ncol(combp)){
    pos = which(combp_ref[i,] < combp[i,j]) ;
    adjpva[i,j] = nrow(as.matrix(pos))/1000 ;
  }
}
adjpva1 = (1-adjpva)
adjpva2 = -log2(adjpva1 +.0001)
scores = list(adjpva,adjpva1,adjpva2)
names(scores) <- c("adjpva","adjpvaraw","adjpvalog")

return(scores)
}



