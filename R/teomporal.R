#' Differential Pathways
#'
#' @param data Adjusted raw p-values matrix
#' @param group group of cell types among which differential pathway analysis needs to be performed.
#'
#' @return A list containing p-value based on wilcoxon rank sum test, fold change based on mean and median
#' @export
#'
#' @examples
#' temporaldif()
temporaldif <- function( data , group)
{
clases = unique(group) ;
data[is.na(data)] = 1 ;
data = -log2(data + 1E-20) ;
wilk = matrix(1, nrow(data) , length(clases)) ;
FC = wilk ; FC1 = wilk ;
rownames(wilk) = rownames(data) ;
for( i in 1:length(clases) )
{
pos = which(group == clases[i]) ;
pos1 = which(group != clases[i]) ;

for (j in 1:length(data[,1]))
{
tri = wilcox.test( data[j, pos] , data[j, pos1] ) ;
wilk[j,i] = tri$p.value ;
FC[j,i] = mean(data[j,pos]) - mean(data[j, pos1]) ;
FC1[j,i] = median(data[j,pos]) - median(data[j, pos1]) ;
}


}

return( list(wilk=wilk, FC=FC, FC1 = FC1)) ;

}

