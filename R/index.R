#' Finding K nearest neighbour
#' @param pathwayscores Raw adjusted p-value matrix
#' @param k Number of k nearest neighbour for each sample or cell
#'
#' @return matrix of indices of nearest neighbour
#' @export
#'
#' @examples
#'
#' index()
  index <- function(pathwayscores,k=5){
  k = FNN::get.knn(t(pathwayscores),k)
  index = k$nn.index

  return(index)

}
