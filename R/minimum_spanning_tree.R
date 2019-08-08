#' Finding minimum spanning tree using shrinked distance matrix
#'
#' @param distance shrinked distance matrix
#'
#' @return Minimum spanning tree
#' @export
#'
#' @examples
#' minimum_spanning_tree()
minimum_spanning_tree <- function(distance){
  corr_g<- igraph::graph.adjacency(as.matrix(distance)  ,mode="undirected",weighted=TRUE,diag = FALSE,add.rownames = NA,add.colnames = NA)
  corr_mst <- igraph::mst(corr_g)

  return(corr_mst)
}
