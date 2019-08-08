#' Global accessibility score calculation
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @param testfile 3 column peak file of genomic coordinates for test data
#' @param referencefile 3 column peak file of genomic coordinates for reference data
#' @param globalaccess_scores pre-calculated global accessibility scores
#'
#' @return global accessibility scores for test data
#' @export
#'
#' @examples
#' global_access()
global_access = function(testfile,referencefile,globalaccess_scores){
g1 = read.table(testfile,sep="\t",header=F,stringsAsFactors = F)
g2 = read.table(referencefile,sep="\t",header=F,stringsAsFactors = F)
globalaccess = read.table(globalaccess_scores,sep=",")

colnames(g1) = c('chr','start','end')
colnames(g2) = c('chr','start','end')

ref <- makeGRangesFromDataFrame(g2)
test <- makeGRangesFromDataFrame(g1)
overlaps <- findOverlaps(test, ref)
overlaps = as.matrix(overlaps)


access = matrix( 0 , length(g1[,1]) , 1)

access[(overlaps[,1])] = globalaccess[(overlaps[,2]),]
access = access[,1]
return(access)
}









