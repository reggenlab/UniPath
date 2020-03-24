#' Generation of foreground file
#'
#' @param arg1 nearestGenes.pl script
#' @param arg2 genomic coordinate file
#' @param arg3 human reference genome file
#' @param arg4 output file
#'
#' @return foreground file
#' @export
#'
#' @examples
#' nearest_gene()
nearest_gene = function(arg1,arg2,arg3,arg4){
  cmd <- paste("perl", arg1, arg2, arg3, arg4)
  return(cmd)
}
