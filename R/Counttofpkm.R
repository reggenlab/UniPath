#' Title Conversion of read count data to fpkm
#'
#' @param countMatrix Raw read count matrix
#' @param length length of genes. Order of length of genes and row names of count data should be same
#'
#' @return fpkm converted matrix
#' @export
#'
#' @examples
#' Counttofpkm()
Counttofpkm = function(countMatrix,length){
fpkm <- countMatrix/length
return(fpkm)
}

