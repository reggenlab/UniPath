#' Creating gradient of colors for showing continuum of co-occurrence of pathways on minimum spanning tree
#'
#' @param score1 scores of first pathway term
#' @param score2 scores of second pathway term
#'
#' @return matrix of gradient of colors for two pathway terms to be visualized on minimum spanning tree
#' @export
#'
#' @examples
#' makecolor()
makecolor  = function (score1, score2)
{
  score1 = as.matrix(score1)
  ro = ncol(score1) ;

  sc1 = as.integer(150 *(score1 - min(score1))/(max(score1) - min(score1) ) + 100 )
  sc2 = as.integer(150 *(score2 - min(score2))/(max(score2) - min(score2) ) + 100 )
  hsc1 = as.character(as.hexmode(sc1)) ;
  hsc2 = as.character(as.hexmode(sc2)) ;
  mystring = matrix("#000000", ro, 1) ;
  for (i in 1:ro)
  {
  mystring[i] = paste("#" , hsc1[i]  , hsc2[i], "00",sep="") ;
}
  return(mystring)
}






