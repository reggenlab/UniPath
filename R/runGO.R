#' Calculating pathway enrichment scores for scATAC-seq profiles
#'
#' @param gmtFile pathway annotation/gene-set file
#' @param BGfile background file
#' @param countFile scATAC-seq count matrix
#' @param method If method is chosen as 1, data normalization is performed using global accessibility scores. If selected method is 2 then local accessibility score-based normalization is performed
#' @param globalaccessibility_scores  global accessibility scores for input count data matrix
#' @param FGfile foreground file
#' @param promoters whether promoters to be used or not for conversion of scATAC-seq profiles to pathway scores. Default is false
#' @param dist  distance to be used for considering nearest gene to a peak
#' @param threshold Peaks above the given threshold value are chosen to be in a foreground set
#' @return A list containing two matrices, One matrix contains p-values based on hypergeometric test and other one is p-values based on binomial test
#' @export
#'
#' @examples
#' runGO()
runGO  <- function ( gmtFile,BGfile,countFile,method,globalaccessibility_scores,FGfile,promoters = FALSE,dist=1000000,threshold=1.25)
{
  ganot = gmtFile
  #ganot = read.table( gmtFile,fill=TRUE,sep="\t",flush=TRUE,stringsAsFactors=F,row.names = 1)

  bgfile = read.table( BGfile , header=FALSE, stringsAsFactors = F) ;
  bgc1 = as.matrix(bgfile[,1]) ;
  bgd1 = as.matrix(bgfile[,2]) ;
  bgc2 = as.matrix(bgfile[,3]) ;
  bgd2 = as.matrix(bgfile[,4]) ;

  for (i in 1:length(bgc1))
  {
    if(bgd1[i] < bgd2[i]) { bgc1[i] = bgc2[i] ; bgd1[i] = bgd2[i]  ; }
  }

  if(promoters == TRUE ) { pos = which(bgd1 < dist) ; }
  if(promoters == FALSE) { pos = which((bgd1 < dist)& (bgd1 > 1000) ) ;}

  bfc = bgc1[pos]
  bfc = as.matrix(bfc)

  p_pi <- matrix(0,nrow(ganot),ncol=1, dimnames=list(rownames(ganot)))

  freq = table(bfc)
  for (i in 1:nrow(ganot)){
    pathway <- ganot[i,]
    genes = pathway[2: ncol(pathway)] ;

    sum_path = 0
    genes = as.matrix(genes)

    for (gene in 1:length(genes)){
      val = genes[1,gene]
      freq_gene = as.integer(freq[val])
      if(is.na(freq_gene) == FALSE){
        sum_path = sum_path+freq_gene
      }
    }
    prob = sum_path/(nrow(bfc))
    p_pi[i] = prob
  }



  #########################  for Foreground analysis

  ######################reading data #################
  if (method==1){
  tagcount = read.table(countFile ,sep=",",header=T,stringsAsFactors=F,row.names = 1)
  #globalaccess = read.table(globalaccess,sep=",")[,1]


  tagcount = (tagcount/(globalaccess+.01))

  fgfile = read.table( FGfile , header=FALSE) ;
  fgc1 = as.matrix(fgfile[,1])
  fgd1 = as.matrix(fgfile[,2])
  fgc2 = as.matrix(fgfile[,3])
  fgd2 = as.matrix(fgfile[,4])

  for (i in 1:length(fgc1))
  {
    if(fgd1[i] < fgd2[i]) { fgc1[i] = fgc2[i] ; fgd1[i] = fgd2[i]  ; }
  }

  if(promoters == TRUE) { pos = which(fgd1 < dist) ;}
  if(promoters == FALSE) { pos = which((fgd1 < dist)& (fgd1 > 1000) ) ; }

  fgd1 = fgd1[pos] ; fgc1 = fgc1[pos] ;
  tagcount = tagcount[pos,] ;


  apvals = matrix( 1, nrow=nrow(ganot), ncol=ncol(tagcount),dimnames=list(rownames(ganot),colnames(tagcount)))

  apvals1 = matrix( 1, nrow=nrow(ganot), ncol=ncol(tagcount),dimnames=list(rownames(ganot),colnames(tagcount)))

  ##################### pathway enrichment ###################
  lengthGO = length(ganot[,1]) ;

  alln = matrix(0, ncol(tagcount), 1) ;
  thr = matrix(0, ncol(tagcount), 1 ) ;

  for (fgs in 1:ncol(tagcount))
  {
    #thr[fgs] = quantile(tagcount[, fgs] , probs=0.80) ;
    #pos1 = which(tagcount[,fgs] > thr[fgs]) ;
    pos1 = which(tagcount[,fgs] > threshold) ;
    alln[fgs] = length(pos1) ;
  }
  fgNum = nrow(tagcount) ;

  for (i in 1:lengthGO )
  {
    mat1 = as.matrix( ganot[i,] ) ;
    mat1 = mat1[ c(-1 , -2)] ;
    pos = which(mat1 != "") ;

    mat1 = as.matrix(mat1[pos]) ;
    match = which( fgc1 %in% mat1 ) ;
    p = p_pi[i] ;
    lfgd1 = fgd1[match] ;
    ltagcount = tagcount[match,] ;
    matchNum = length(match) ;


    for( fgs in 1:ncol(tagcount)){
      pos = which(ltagcount[,fgs] > threshold) ;
      #pos = which(ltagcount[,fgs] > thr[fgs]) ;

      k = length(pos) ;
      n = alln[fgs] ;
      apvals[i,fgs] <- (1-pbinom(k , n, p));
      apvals1[i,fgs] <- (1-phyper(k , matchNum, (fgNum - matchNum) , n ));

    }

  }
  }

if (method==2){
  data = read.table(countFile ,sep=",",header=T,stringsAsFactors=F,row.names = 1)
  data = as.matrix(data)
  tagcount = normalize.quantiles(data)
  colnames(tagcount) = colnames(data)
  rownames(tagcount) = rownames(data)



  mean = apply(tagcount,1,mean)
  tagcount = (tagcount/(mean+.01))

  fgfile = read.table( FGfile , header=FALSE) ;
  fgc1 = as.matrix(fgfile[,1])
  fgd1 = as.matrix(fgfile[,2])
  fgc2 = as.matrix(fgfile[,3])
  fgd2 = as.matrix(fgfile[,4])

  for (i in 1:length(fgc1))
  {
    if(fgd1[i] < fgd2[i]) { fgc1[i] = fgc2[i] ; fgd1[i] = fgd2[i]  ; }
  }

  if(promoters == TRUE) { pos = which(fgd1 < dist) ;}
  if(promoters == FALSE) { pos = which((fgd1 < dist)& (fgd1 > 1000) ) ; }

  fgd1 = fgd1[pos] ; fgc1 = fgc1[pos] ;
  tagcount = tagcount[pos,] ;


  apvals = matrix( 1, nrow=nrow(ganot), ncol=ncol(tagcount),dimnames=list(rownames(ganot),colnames(tagcount)))

  apvals1 = matrix( 1, nrow=nrow(ganot), ncol=ncol(tagcount),dimnames=list(rownames(ganot),colnames(tagcount)))

  ##################### pathway enrichment ###################
  lengthGO = length(ganot[,1]) ;

  alln = matrix(0, ncol(tagcount), 1) ;
  thr = matrix(0, ncol(tagcount), 1 ) ;

  for (fgs in 1:ncol(tagcount))
  {
    #thr[fgs] = quantile(tagcount[, fgs] , probs=0.80) ;
    #pos1 = which(tagcount[,fgs] > thr[fgs]) ;
    pos1 = which(tagcount[,fgs] > threshold) ;
    alln[fgs] = length(pos1) ;
  }
  fgNum = nrow(tagcount) ;

  for (i in 1:lengthGO )
  {
    mat1 = as.matrix( ganot[i,] ) ;
    mat1 = mat1[ c(-1 , -2)] ;
    pos = which(mat1 != "") ;

    mat1 = as.matrix(mat1[pos]) ;
    match = which( fgc1 %in% mat1 ) ;
    p = p_pi[i] ;
    lfgd1 = fgd1[match] ;
    ltagcount = tagcount[match,] ;
    matchNum = length(match) ;


    for( fgs in 1:ncol(tagcount)){
      pos = which(ltagcount[,fgs] > threshold) ;
      #pos = which(ltagcount[,fgs] > thr[fgs]) ;

      k = length(pos) ;
      n = alln[fgs] ;
      apvals[i,fgs] <- (1-pbinom(k , n, p));
      apvals1[i,fgs] <- (1-phyper(k , matchNum, (fgNum - matchNum) , n ));

    }
  }

}





  scores=list(apvals,apvals1)
  names(scores) <- c("binomial","hypergeometeric")
  return(scores)

}

