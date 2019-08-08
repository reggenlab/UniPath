#' Differential coenriched pathways
#'
#' @param data Raw p-value pathway score matrix
#' @param group group 1 for cells of interest and rest in group 2
#'
#' @return A list containing matrix of p-values and difference among the groups
#' @export
#'
#' @examples
#' difcoccur()
difcoccur <- function( data , group )
{
data = -log2(data + 1E-20) ;
cols = ncol(data) ;
pos = which(group == 1) ;
pos1 = which(group== 2) ;
co1 = cor( t(data[ , pos]) , method="spearman") ;
co2 = cor( t(data[, pos1]) , method="spearman") ;
dif = abs(co1 - co2) ;
dif1 = co1 - co2 ;
rcounts = 0 * co1 ;

for (j in 1:1000)
{
npos = sample( 1:cols) ;
co1 = cor(t(data[, npos[pos]]) , method="spearman") ;
co2 = cor(t(data[, npos[pos1]]) , method="spearman") ;
rdif = abs( co1 - co2) ;
difdif = (rdif - dif)  ; pos3 = which(difdif > 0) ;
rcounts[pos3] = rcounts[pos3] + 1 ;
}
rcounts = rcounts / 1000 ;


toret = list( pval = rcounts, dif = dif1 ) ;
return(toret) ;

}






