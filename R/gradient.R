#' Showing gradient or continuum of pathways on minimum spanning tree
#'
#' @param pathwayfile adjusted p-value matrix
#' @param term Pathway term for which gradient needs to be plotted
#'
#' @return gradient of colors for specific pathway term
#' @export
#'
#' @examples
#' gradient()
gradient = function(pathwayfile,term){


#pathways = read.table(pathwayfile,sep=",",header=T,stringsAsFactors = F,row.names = 1)
tri = as.matrix(pathwayfile[term,])

colfunc<-colorRampPalette(c( 'red',"yellow"))
ntri = as.integer(9*(tri - min(tri))/ (max(tri) - min(tri)) +1)
col10 = as.matrix(colfunc(10))
coltri = col10[ntri]
return(coltri)
}


