#' Pathway co occurrence
#'
#' @param nulldata Pathway score matrix corresponding to null model
#' @param PathwayScores  Adjusted raw p-values matrix
#'
#' @return A list containing matrices of p-values and correlation for pathway pairs in particular cell type
#' @export
#'
#' @examples
#' coocurrence()
cooccurrence = function(nulldata,PathwayScores) {
  PathwayScores = -log2(PathwayScores + 1E-20) ;

  corr1 = cor(t(PathwayScores) , method="spearman")
  counts = 0*corr1

  for (j in 1:1000)
  {
    npos = sample( 1:ncol(nulldata),500)
    corr2 = cor(t(nulldata[, npos]) , method="spearman")
    diff = abs(corr1) - abs(corr2)
    pos = which(diff < 0) ;
    counts[pos] = counts[pos] + 1
  }
  counts = counts / 1000
  output = list( pval = counts, correlation = corr1 )

  return (output)
}
