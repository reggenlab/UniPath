#' Getting cluster numbers for each of the nearest neighbor of a cell
#'
#' @param pathwayscores Raw adjusted p-value matrix
#' @param index indices of nearest neighbours obtained from index function
#' @param clusters Clusters obtained from hierarchal clustering
#'
#' @return Matrix having cluster or class number for each of the top nearest neighbor of individual cell
#' @export
#'
#' @examples
#' KNN()
KNN <- function(pathwayscores,index,clusters){
KNN = matrix(0, nrow = ncol(pathwayscores), ncol = ncol(index))
Clusters = as.vector(unlist(clusters))
for (j in 1:ncol(index)){
  Index.Col = index[, j]
  for (i in 1:nrow(index)){
    temp = Index.Col[i]
    KNN[i, j] = Clusters[temp]
  }
}
return(KNN)
}

