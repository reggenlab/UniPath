#' Two level shrinkage of distance matrix based on nearest neighbour indices and belongingness of cells to same class
#'
#' @param dist Distance matrix used for hierarchical clustering
#' @param class Matrix with number of times cells in same class have top k neighbors in other classes
#' @param clusters clusters obtained from hierarchical clustering
#'
#' @return shrinked distance matrix
#' @export
#'
#' @examples
#' distance()
distance <- function(dist,class,clusters){

correlation = as.matrix(1/dist)
distance_matrix <- as.matrix(dist)
lambda = 3

for (i in 1:nrow(class))
{
  class[i,i] = 0 ;
}
class1 = class / (rowSums(class) + 0.01) ;

for (i in 1:nrow(distance_matrix))
{
  for (j in 1:ncol(distance_matrix))
  {
    if(clusters[i] == clusters[j])
    {
      correlation[i,j] <- (correlation[i,j]*lambda)
    }

    else
      correlation[i,j] <- (correlation[i,j]*(0.25+(class1[clusters[i],clusters[j]])*lambda*0.5)) ;

  }
}

dist1 = (1-correlation)
shrinkdist = as.dist(dist1, diag = FALSE, upper = FALSE)
return(shrinkdist)
}




