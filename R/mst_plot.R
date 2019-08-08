#' Plotting minimum spanning tree
#' @import netbiov
#' @import igraph
#' @param x
#' @param layout.function
#' @param colors
#' @param mst.edge.col
#' @param vertex.color
#' @param tkplot
#' @param expression
#' @param v.size
#' @param e.size
#' @param mst.e.size
#' @param edge.col.wt
#' @param v.lab
#' @param e.lab
#' @param bg
#' @param v.lab.cex
#' @param e.lab.cex
#' @param v.lab.col
#' @param lab.dist
#' @param e.lab.col
#' @param v.sf
#' @param sf
#' @param e.arrow
#' @param layout.overall
#'
#' @return mst plot
#' @export
#'
#' @examples
mst.plot.mod <- function(x, layout.function=NULL,colors=NULL,
                         mst.edge.col="white", vertex.color = "skyblue",tkplot=FALSE,
                         expression=NULL, v.size=FALSE, e.size=FALSE, mst.e.size=1,
                         edge.col.wt=NULL, v.lab=FALSE, e.lab=NULL, bg="black",v.lab.cex=0.5,
                         e.lab.cex=0.5,v.lab.col="blue",lab.dist=0, e.lab.col="blue",v.sf=c(3,12),
                         sf=0, e.arrow=.2, layout.overall=NULL){
  if(is.null(x)){
    stop("please input a valid igraph object")
  }
  if(!is.igraph(x)){
    stop("please input a valid igraph object")
  }
  gtmp1 <- .set.mst.graph.attributes_mod(x, v.lab)

  g <- gtmp1[[1]]
  v.label <- gtmp1[[2]]
  temp1 <- .set.mst.size_mod(g, e.size, v.size, v.sf)

  e.s <- temp1[[1]]
  v.s <- temp1[[2]]
  gtemp <- .get.mst.info(g)
  mst <- gtemp[[2]]
  gt <- gtemp[[1]]
  V(mst)$name <- paste("g", c(1:vcount(mst)), sep="")
  #########print(ecount(mst))
  gtmp1 <- .set.module.info_mst(mst)
  split.vector <- gtmp1[[1]]
  #####print(split.vector)
  v.id <-  gtmp1[[2]]
  mod.list <- gtmp1[[3]]
  layout.function<-.set.layout.function_mst(layout.function,split.vector)
  ###print(layout.function)
  ##########print(mod.list)
  gtmp2<-.set.graph.attributes_mst(mst,layout.function,split.vector,v.id)
  layouts <- gtmp2[[1]]
  #######print(layouts)
  lgsplitx <- gtmp2[[2]]
  abstract.graph <- .absg_mst(mst, mod.list)

  crd <- .getcrd.mod_mst(layouts,layout.overall,sf,
                         abstract.graph=abstract.graph)

  crd <- crd[V(mst)$name,]
  ##########print(nrow(crd))
  crd <- cbind(crd, c(1:(vcount(mst))))
  ##########print(crd)

  gtemp <- g
  if(is.null(edge.col.wt) && (length(edge.col.wt)<ecount(g))){
    temp2 <- .set.edge.color_mstmod(crd, mst, gt, colors)
    colid <- temp2[[1]]
    gt <- temp2[[2]]
    if(length(colors)==0){
      cole <- colors()[colid]
    }
    if(length(colors)>0){cole <- colid}
    e.s[which(cole=="grey19")] <- mst.e.size
    if(!is.null(mst.edge.col)){
      cole[which(cole=="grey19")] <- mst.edge.col
    }
  }
  else{
    cole <- .set.edge.col.wt(g,gt,mst,edge.col.wt,mst.e.size,e.size)
    gt <- cole[[3]]
    e.s <- cole[[2]]
    #########print(e.s)
    cole <- cole[[1]]
  }
  v.colors <- .set.mst.node.col_mod(vertex.color, gt, expression)
  if(is.null(e.lab)|| (e.lab==FALSE)){
    e.lab <- NA
  }
  else{
    if(e.lab==TRUE){
      e.lab <- NULL
    }
  }
  gparm <- list(g= gt, layout = crd[,1:2], vertex.color=v.colors,
                edge.color=cole,vertex.size=v.s,edge.arrow.size=e.arrow,
                vertex.label.cex= v.lab.cex, vertex.label=v.label, lab.color=v.lab.col,
                lab.dist=lab.dist, vertex.frame.color=v.colors, edge.width=e.s, bg=bg,
                edge.label=e.lab)


  class(gparm) <- "netbiov"
  if(tkplot){
    tkplot.netbiov(gparm)
  }
  else{
    plot.netbiov(gparm)
  }
  gparm
}

.dist.mst <- function(x1,y1,x2,y2){
  return(sqrt( (x1 - x2 )^2 + (y1 - y2)^2 ))
}

#n = 500
#g <- watts.strogatz.game(1, n, 5, 0.05)
#g <-as.directed(g, mode = "arbitrary")
#g <- simplify(g)


.rgbToHex.mod <- function(n){
  #R <-  G
  #B  <- 255
  #G <- 255
  #R=(n-1)
  if(n<255){
    R = 255
    G = 0
    B = n
  }
  else{
    R=255-(n-255)
    G=0
    B=255
  }

  #G=(n-20)
  #R=255-(n-10);
  #B=13


  k1 <- .toHex.mod(R)
  k2 <- .toHex.mod(G)
  k3 <- .toHex.mod(B)
  k <- paste("#",k1,k2,k3, sep="")
  return(k)
}

.toHex.mod <- function(n){
  hexstr <- "0123456789ABCDEF"
  if(is.na(n)) {return("00")}
  n = max(0,min(n,255));
  k1 <- strsplit(hexstr,"")[[1]][((n-n%%16)/16)+1]
  k2 <- strsplit(hexstr,"")[[1]][(n%%16)+1]
  return(paste(k1,k2,sep=""))
}
############################################################

.set.mst.graph.attributes_mod <- function(g,v.lab ){
  if(is.directed(g)){
    mode = "directed"
  }
  else{
    mode="undirected"
  }
  adj <- get.adjacency(g, names=FALSE)
  gname <- V(g)$name
  if(!is.null(V(g)$label)){
    v.label <- V(g)$label
  }
  if(vcount(g)==length(v.lab)){
    v.label <- v.lab
  }
  if(class(v.lab)=="logical"){
    if(v.lab==TRUE){
      if(!is.null(V(g)$label)){
        v.label <- V(g)$label
      }
      else{
        v.label <- V(g)$name
        if(length(gname)<vcount(g) )
        {
          v.label <- c(1:(vcount(g)))
        }
      }

    }
    else{
      v.label <- NA
    }
  }
  ##########print(v.label)
  #g <- graph.adjacency(adj,mode=mode)
  #g <- simplify(g)
  if(!is.null(v.label)&&(length(v.label)==vcount(g))){
    V(g)$label <- v.label
  }
  E(g)$label <-E(g)
  #########print(g)
  list(g, v.label)
}

.set.mst.size_mod <- function(g, e.size, v.size, sf){
  v.s = 2.5
  #print(v.size)
  if(class(v.size)=="logical"){
    v.s=2.5
  }
  if((class(v.size)=="numeric")||(class(v.size)=="integer")){
    if(length(v.size)< vcount(g)){
      v.s = v.size[1]
    }
    else{
      if(length(sf)==1){
        v.s = (v.size/max(v.size))*sf
        v.s[which(v.s<=.5)] <- .5
      }
      else{
        v.s <- v.size
        #v.s=min(sf)+(((v.s-min(v.s))/(max(v.s)-min(v.s)))*(max(sf)-min(sf)))
        v.s = min(sf)+((v.s-min(v.s))/(max(v.s)-min(v.s)))*(max(sf)-min(sf))

        ##########print(v.s)
      }
    }
  }
  if(class(e.size)=="logical"){
    e.s=rep(.5,ecount(g))
  }
  if(class(e.size)=="numeric"){
    if(length(e.size)< ecount(g)){
      e.s = rep(e.size[1],ecount(g))
    }
    else{
      e.s = (e.size[1:ecount(g)]/max(e.size[1:ecount(g)]))*12
      e.s[which(e.s<=.5)] <- .5
    }
  }

  list(e.s, v.s)

}
.get.mst.info <- function(g){
  mst <- minimum.spanning.tree(g)
  ##########print(ecount(mst))
  gt <- graph.difference(g,mst)
  ##########print(ecount(gt))
  lab <- sort(setdiff(E(g)$label,E(mst)$label))
  ##########print(ecount(gt))
  list(gt, mst)

}

.set.edge.color_mstmod <- function(crd=NULL, mst=NULL, gt=NULL, colors=NULL){
  ##########print(vcount(g))
  if(length(colors)==0){
    colcode<- grep("orange",colors())
    #colcode <- (heat.colors(20))
    ##########print(colcode)
    colid1 <-rep(280, ecount(mst))
  }
  if(length(colors)>0){
    colcode <- c(1:length(colors))
    names(colcode) <- colors
    colid1 <-rep("grey19", ecount(mst))
  }
  colid <- c()
  if(ecount(gt)>0){
    edgelist <- get.edgelist(gt, names=FALSE)
    d <-c()
    colid<-c()
    for(i in 1:length(edgelist[,1])){
      e1 <- which(crd[,3]==edgelist[i,1])
      e2 <- which(crd[,3]==edgelist[i,2])
      d <-append(d,.dist.mst(crd[e1,1],crd[e1,2],crd[e2,1],crd[e2,2]))
    }
    ##########print(d)
    col <- rep(round(((max(d)) - (min(d)))/
                       (length(colcode)),digits=5),(length(colcode)))
    ##########print(col)
    col <- cumsum(col)
    col <- append(min(d),col+min(d))
    ##########print(col)
    if(col[length(col)]< max(d))
    {
      #col[length(col)] <- col[length(col)]+2
      col[length(col)] <- max(d)+1
    }
    for(i in 1:length(d))
    {
      for(j in 1:(length(col)-1)){
        if((d[i]>=col[j])&&(d[i]<col[j+1])){
          if(length(colors)==0){
            colid <- append(colid, colcode[j])
          }
          if(length(colors)>0){
            colid <-append(colid,names(colcode)[j])
          }
        }
      }
    }
  }
  tmpelab <- c(E(gt)$label, E(mst)$label)
  gt <- add.edges(gt, t(get.edgelist(mst, names=FALSE)))
  E(gt)$label <- tmpelab
  list(c(colid, colid1), gt)
}

.set.mst.node.col_mod <- function(vertex.color, gt, expression){
  if(length(vertex.color)<vcount(gt)){
    v.colors <-  rep(vertex.color, vcount(gt))[1:vcount(gt)]
  }
  else{
    v.colors <- vertex.color
  }
  if(!is.null(expression)&&(length(expression)==vcount(gt))){
    sft <- c(1,510)
    rank.wt = min(sft)+((expression-min(expression))/
                          (max(expression)-min(expression)))*(max(sft)-min(sft))
    #rank.wt<-((rank.wt-min(rank.wt))/(max(rank.wt) - min(rank.wt)))*509 + 1
    rank.wt <- round(rank.wt, digits=0)
    #########print(as.vector(rank.wt))
    exp.col <- sapply(rank.wt,.rgbToHex.mod)
    vertex.color <- exp.col
    if(is.null(vertex.color)||(length(vertex.color)<2))
    {vertex.color <- heat.colors(max(rank.wt))}
    v.colors <- vertex.color
  }

  v.colors
}
.set.edge.col.wt <- function(g, gt, mst, edge.col.wt, mst.e.size, e.size){
  gt <- add.edges(gt, t(get.edgelist(mst, names=FALSE)))
  #########print(get.edgelist(g, names=FALSE))
  l1 <-  apply(get.edgelist(g, names=FALSE), 1,
               function(x)paste(as.vector((x)), collapse="_"))
  l2 <-  apply(get.edgelist(gt, names=FALSE), 1,
               function(x)paste(as.vector((x)), collapse="_"))
  l3 <-  apply(get.edgelist(mst, names=FALSE), 1,
               function(x)paste(as.vector((x)), collapse="_"))
  names(edge.col.wt) <- l1
  edge.col.wt <- edge.col.wt[l2]
  rank.wt <- (edge.col.wt)
  rank.wt <- ((rank.wt-min(rank.wt))/(max(rank.wt)-min(rank.wt)))*509+1
  rank.wt <- round(rank.wt, digits=0)
  #########print(as.vector(rank.wt))
  e.wt.col <- sapply(rank.wt,.rgbToHex.mod)
  e.col <- e.wt.col
  if(class(e.size)=="logical"){
    e.s=rep(.5,ecount(g))
  }
  if(class(e.size)=="numeric"){
    if(length(e.size)< vcount(g)){
      e.s = rep(e.size[1],ecount(g))
    }
    else{
      e.s = (e.s/max(e.s))*12
      e.s[which(e.s<=.5)] <- .5
    }
  }
  names(e.s) <- l1
  e.s <- e.s[l2]
  e.s[l3] <- mst.e.size
  list(e.col, e.s, gt)
}


.set.graph.attributes_mst <- function(g,layout.function, split.vector, v.id,
                                      exp.by.module=NULL, expression=NULL, color.in=NULL, col.grad=NULL,
                                      arrange.bydegree=FALSE){
  #print(split.vector)
  lgsplit <-list()
  lgsplitx <- list();
  for(i in 1:(length(split.vector)-1))
  {
    #print(split.vector[i]+1)
    #print(split.vector[i+1])
    if((split.vector[i]+1)<split.vector[i+1]){
      g.temp<-induced.subgraph(g,v.id[(split.vector[i]+1):(split.vector[i+1])])
    }
    else{
      g.temp<-induced.subgraph(g,v.id[(split.vector[i]+1):(split.vector[i]+1)])
    }
    lgsplitx[[i]] <- g.temp
    dg = degree(g, V(g.temp)$name)
    exp.mod = FALSE
    #########################################################
    if(class(exp.by.module)=="logical"){
      if(exp.by.module && !(is.null(expression))){
        dg <- expression[V(g.temp)$name]
        if(length(dg)< vcount(g.temp)){
          dg <- degree(g, V(g.temp)$name)
        }
        #dg <- (abs(dg)*510)/max(abs(dg))
        dg<-(((510 - 1)*(abs(dg)-min(abs(dg))))/
               (max(abs(dg))-min(abs(dg))))+1
        ##########print(head(dg))
        exp.mod = TRUE
      }
    }
    if(class(exp.by.module)=="numeric"){
      if((length(exp.by.module)>0) && !(is.null(expression))){
        if(i %in% exp.by.module){
          dg <- expression[V(g.temp)$name]
          exp.mod=TRUE
        }
        if(length(dg)< vcount(g.temp)){
          dg <- degree(g, V(g.temp)$name)
        }
        #dg <- (abs(dg)*510)/max(abs(dg))
        dg <- (((510 - 1)*(abs(dg)-min(abs(dg))))/
                 (max(abs(dg))-min(abs(dg)))) + 1

      }
    }
    ##########print(exp.mod)
    ###########################################################
    lgsplit[[i]]<-list(g =g.temp,dg = dg,color.input=color.in[i],
                       col.grad = col.grad[[i]], layout.function = layout.function[[i]],
                       arrange.bydegree=arrange.bydegree, exp.mod = exp.mod)
  }
  ###print(lgsplit[[1]]$layout.function)
  layouts <- lapply(lgsplit, .get.coord.mod_mst)

  list(layouts, lgsplitx)
}



.set.module.info_mst <- function(g){
  memb <- .mod.function(g)
  un <- unique(memb)
  #print(memb)
  mod.list <- sapply(sort(un), function(x)which(x==memb))
  #mod.list <- sapply(mod.list, function(x)x-1)
  ############print(mod.list)
  vid <- as.numeric(unlist(mod.list))
  v.rest <- setdiff(as.vector(V(g)), vid)
  v.id <- c(vid, v.rest)
  if(length(v.rest)>0){
    split.vector<-c(0,unlist(lapply(mod.list,length)),length(v.rest))
  }
  else{
    split.vector <- c(0,unlist(lapply(mod.list, length)))
  }
  split.vector <- cumsum(split.vector)
  list(split.vector, v.id, mod.list)

}



.get.coord.mod_mst <- function(g){
  layout.function <- g$layout.function
  g <- g[[1]]
  ###print(is.null(layout.function))
  ###print(class(layout.function))
  if(is.null(layout.function)||(class(layout.function)!="function")){
    crd1 <- layout_as_tree(g,
                           root=as.vector(V(g)[which.max(degree(g))]), circular=TRUE)
    crd1[(is.nan(crd1))] <- 0
  }
  else{
    tryCatch({

      crd1 <- layout.function(g)
    }, error = function(ex) {
      crd1 <- layout_as_tree(g,
                             root=as.vector(V(g)[which.max(degree(g))]),circular=TRUE)
      crd1[(is.nan(crd1))] <- 0
    })

  }
  rownames(crd1) <- V(g)$name

  crd1
}


.getcrd.mod_mst <- function(lst, layout.overall=NULL,sf=50, scale.module=NULL,
                            abstract.graph=NULL){
  ###############print(is.null(abstract.graph))
  #print(sf)
  cnt <- c()
  for(i in 1:length(lst)){
    cnt <- c(cnt, nrow(lst[[i]]))
  }
  cnt <- cbind(c(1:length(cnt)), cnt)

  if(is.null(abstract.graph)){
    cnt <- cnt[order(cnt[,2], decreasing = TRUE),]
    lstnew <- list()
    for(i in 1:length(lst)){
      lstnew[[i]] <- lst[[cnt[i,1]]]
    }
    lst <- lstnew

  }
  scale.module <- scale.module[cnt[,1]]
  ####################print(head(scale.module))

  dst <- function(x){return(dist(rbind(c(0,0), x)))}
  k  <- c()
  minn <- function(x){return(min(x,50))}
  mxx <- function(x){if(x<.99){return(1)}else{return(x)}}
  if(is.null(scale.module)){
    for(i in 1:length(lst)){
      k <- append(k, nrow(lst[[i]]))
    }
    k <- sqrt(k)
    #k <- sapply(k,minn)
    ####################print(k)
  }
  else{
    if(length(scale.module)!=length(lst)){
      k <- rep(scale.module, length(lst))
      k <- k[1:length(lst)]
    }
    else{
      k <- scale.module
    }
  }
  k <- sapply(k,minn)

  max.k <- max(k)
  #####print(is.null(abstract.graph))
  if(is.null(layout.overall)||(class(layout.overall)!="function")){
    if(!is.null( abstract.graph)){

      gx <- graph.adjacency( abstract.graph, mode="undirected")
      crdf <- layout_with_fr(gx)
      #####print("hello")
      #plot(gx, layout=crdf)

    }
    else{
      crdf<-layout_with_fr(graph.empty(length(lst)))
    }
  }
  else{
    if(!is.null( abstract.graph)){
      gx<-graph.adjacency(abstract.graph,mode="undirected")
      crdf <- layout.overall(gx)
      #########print(crdf)
    }
    else{crdf <- layout.overall(graph.empty(length(lst)))}
  }
  scale.x <- max(abs(crdf[,1])); scale.y <- min(abs(crdf[,1]))
  scale.x <- scale.y <- max(scale.x, scale.y)
  crdf <- layout.norm(crdf, xmax=scale.x+max.k+20+sf,
                      xmin=-1*(scale.x+max.k+20+sf), ymax=scale.y+max.k+20+sf,
                      ymin=-1*(scale.y+max.k+20+sf) )
  #crdf <- layout.norm(crdf, xmax=scale.x,
  #xmin=-1*(scale.x), ymax=scale.y, ymin=-1*(scale.y) )

  rownames(crdf) <- paste("g", c(1:length(lst)),sep="")
  d  <- apply(crdf,1, dst)
  names(d) <- rownames(crdf)
  if(is.null( abstract.graph)){
    d <- sort(d, decreasing=FALSE)
    crdf <- crdf[names(d),]

  }
  crd <- list()
  crdall <- c()
  for(i in 1:length(lst)){
    tmp <- lst[[i]]
    tmp1 <- matrix(as.numeric(tmp[,c(1:2)]), nrow=nrow(tmp))
    rownames(tmp1) <- rownames(tmp)
    tmp1 <- layout.norm(tmp1, (crdf[i,1]-k[i]),(crdf[i,1]+k[i]),
                        (crdf[i,2]-k[i]),(crdf[i,2]+k[i]))

    tmp1[(is.nan(tmp1[,1])),1] <- min(tmp1[!is.nan(tmp1[,1]),1])
    tmp1[(is.nan(tmp1[,2])),2] <- min(tmp1[!is.nan(tmp1[,2]),2])

    kk1 <- which(tmp1[,1]==Inf)
    kk2 <- which(tmp1[,2]==Inf)
    if(length(kk1)>0){
      tmp1[kk1,1] <- crdf[i,1]+rnorm(1)
    }
    if(length(kk2)>0){
      tmp1[kk2,2] <- crdf[i,2]+rnorm(1)
    }
    ####################print(tmp1)
    crd[[i]] <- tmp1
    crdall <- rbind(crdall, tmp1)
  }
  #crd <- cbind(crdall, c(0:(length(crdall[,1])-1)))
  crdall
}






.absg_mst <- function(g, mod.list){
  adj <- get.adjacency(g)
  ###########print(mod.list)
  m <- matrix(0, nrow=length(mod.list), ncol=length(mod.list))
  for(i in 1:(length(mod.list)-1)){
    for(j in (i+1):length(mod.list)){
      if(!is.directed(g)){
        xx <- adj[(mod.list[[i]]), (mod.list[[j]])]
        if(sum(xx) > 0){
          m[i,j]<- m[j,i] <- 1
        }
      }
      else{
        xx1 <-  adj[(mod.list[[i]]), (mod.list[[j]])]
        xx2 <-  adj[(mod.list[[j]]), (mod.list[[i]])]
        if(sum(xx1) > 0){
          m[i,j] <- 1
        }
        if(sum(xx2) > 0){
          m[j, i] <- 1
        }
      }
    }
  }
  m
}

.set.layout.function_mst <- function(layout.function, split.vector){
  ##########print(v.id)
  if(length(layout.function)<(length(split.vector)-1)){
    if(class(layout.function)=="function"){
      layout.function <- c(layout.function)
    }
    ly.tmp <- as.list(layout.function)
    ##########print(length(ly.tmp))
    layout.function <- rep(ly.tmp, (length(split.vector)-1))
    #  #########print(length(layout.function))
    layout.function <- layout.function[1:(length(split.vector)-1)]

  }
  layout.function

}


.mod.function <- function(g){
  if(!is.directed(g)){
    fc <- multilevel.community(g)
    memb <- fc$membership
  }
  else{
    memb <- walktrap.community(g)$memb
  }

  memb
}
