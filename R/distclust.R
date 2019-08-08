#' Hierarchal clustering of pathway score matrix
#'
#' @param pathwayscores Raw adjusted p-value matrix
#' @param n Number of clusters required for pseudo temporal ordering
#'
#' @return Distance matrix and Clusters
#' @export
#'
#' @examples
#' dist_clust()
dist_clust <- function(pathwayscores,n){
  dis1 = vegan::vegdist(t(pathwayscores),method="euclidean")
  hc1 = hclust(dis1,method="ward.D2")
  clust1 = cutree(hc1,n)
     distclust  = list(dis1,clust1)
     names(distclust) = c("distance","clusters")
  return(distclust)

}
