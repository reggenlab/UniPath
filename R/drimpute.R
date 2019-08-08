#' Imputation of scATAC-seq profiles
#'
#' @param countFile scATAC-seq count matrix
#'
#' @return imputed count matrix
#' @export
#'
#' @examples
#' drimpute()
drimpute = function(countFile){
X = read.table(countFile,sep=",",header=T,stringsAsFactors = F,row.names = 1)
X = as.matrix(X)
X.log <- log(X + 1)
X.log = X.log[,apply(X.log, 2, var, na.rm=TRUE) != 0]
X.imp <- DrImpute::DrImpute(X.log)
return(X.imp)
}

