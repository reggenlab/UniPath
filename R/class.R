
#' Finding how many times nearest neighbor of cells in each class are belonging to different cluster or class.
#'
#' @param  Clusters  clusters obtained from hierarchal clustering
#' @param  KNN Matrix obtained from KNN function with cluster number for top nearest neighbors for each of the cell
#'
#' @return Matrix with number of times cells in same class have top k neighbors in other classes
#' @export
#'
#' @examples
#' class1()
class1 <- function(clusters,KNN){
clusters = unlist(clusters)
class1 = matrix(0,max(clusters),max(clusters))
for (i in 1:nrow(KNN)){
  idx1 = clusters[i] ;
  for (j in 1:ncol(KNN)) {
    idx2 = KNN[i,j] ;
    class1[idx1, idx2] = class1[idx1, idx2] + 1 ;
  }
}
return(class1)
}
