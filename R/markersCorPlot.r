markersCorPlot <- function(cross, population, map=c("genetic","physical"), cmBetween=25, verbose=TRUE){
  map <- defaultCheck.internal(map,"map",2,"genetic")
	if(map=="genetic"){
    cur_map <- population$maps$genetic
  }else{
    cur_map <- population$maps$physical
  }
  if(is.null(cur_map)) stop("no ",map," map provided!")  
  offsets1 <- getPopulationOffsets.internal(population,cur_map,cmBetween)
  offsets2 <- getChrOffsets.internal(cross,cmBetween)

  global_offset <- NULL
  for(x in 1:length(offsets1)){
    global_offset <- c(global_offset,max(offsets1[x],offsets2[x]))
  }
  global_offset<-c(0,global_offset)
  
  sum_gl_off <- NULL
  for(x in 1:length(global_offset)){
    sum_gl_off <- c(sum_gl_off,sum(global_offset[1:x]))
  }


  mloc_original <- getMarkerOffsets(cross,global_offset,cmBetween)
  mloc_o <- getMarkerOffsetsFromMap(cur_map,global_offset,cmBetween)

  m_max <- max(mloc_o,mloc_original)
  m_min <- min(mloc_o,mloc_original)

  back <- abs(chromCorMatrix(cross,population,map))
  plot(c(m_min,m_max),c(m_min,m_max),type='n')
  for(i in 1:(length(sum_gl_off)-1)){
    for(j in 1:(length(sum_gl_off)-1)){
        cur_col <- 1-back[i,j]/max(back)
        rect(sum_gl_off[i],sum_gl_off[j],sum_gl_off[i+1],sum_gl_off[j+1],lty=0,col=rgb(cur_col,cur_col,cur_col))
    }
  }
  points(cbind(mloc_o,mloc_o),pch=21,col="red",lwd=4)
  points(cbind(mloc_original,mloc_original),pch=20,col="green",lwd=2)
  abline(v=sum_gl_off[-length(sum_gl_off)],lty=2)
  abline(h=sum_gl_off[-length(sum_gl_off)],lty=2)
  invisible(global_offset)
} 

getChrOffsets.internal <- function(cross, cmBetween){
  offsets <- unlist(lapply(pull.map(cross),max))
  offsets <- offsets+cmBetween
  offsets <-c(0,offsets)
  offsets
}


#From IQTL by Danny Arends, SHOULD NOT MODIFY
getMarkerOffsets <- function(cross, offsets, cmBetween=25){
  if(missing(offsets))offsets <- getChrOffsets(cross,cmBetween)
  cnt <- 1
  myoffsets <- NULL
  for(x in nmar(cross)){
    myoffsets <- c(myoffsets,rep(sum(offsets[1:cnt]),x))
    cnt <- cnt+1
  }

  mlocations <- myoffsets + as.numeric(unlist(pull.map(cross)))
  mlocations
}

getMarkerOffsetsFromMap <- function(map, offsets, cmBetween=25){
  cnt <- 1
  myoffsets <- NULL
  for(x in table(map[,1])){
    myoffsets <- c(myoffsets,rep(sum(offsets[1:cnt]),x))
    cnt <- cnt+1
  }

  mlocations <- myoffsets + as.numeric(map[,2])
  mlocations
}



ascendingMaptoJigSawMap <- function(population,cur_map,verbose=FALSE){
  map <- defaultCheck.internal(map,"map",2,"genetic")
	if(map=="genetic"){
    cur_map <- population$maps$genetic
  }else{
    cur_map <- population$maps$physical
  }
  if(is.null(cur_map)) stop("no ",map," map provided!")
  for(x in unique(cur_map[,1])){
    if(verbose) cat("Processing chromosome:",x,"\n")
    cur_beg <- min(cur_map[which(cur_map[,1]==x),2])
    cur_map[which(cur_map[,1]==x),2] <- cur_map[which(cur_map[,1]==x),2]-cur_beg
  }
  if(map=="genetic"){
     population$maps$genetic <- cur_map
  }else{
    population$maps$physical <- cur_map
  }
  invisible(population)
}

getPopulationOffsets.internal <- function(population, cur_map, cmBetween){
  minima <- NULL
  for(x in unique(cur_map[,1])){
    minima <- c(minima,max(cur_map[which(cur_map[,1]==x),2]))
  }
  minima <- minima + cmBetween
  c(0,minima)
  invisible(minima)
}

chromCorMatrix <- function(cross,population,map=c("genetic","physical")){
    map <- defaultCheck.internal(map,"map",2,"genetic")
	if(map=="genetic"){
    cur_map <- population$maps$genetic
  }else{
    cur_map <- population$maps$physical
  }
  s <- proc.time()
  if(is.null(cur_map)) stop("no ",map," map provided!")  
  result <- matrix(0,nchr(cross),length(unique(cur_map[,1])))
  for(i in 1:nchr(cross)){
      cur_new <- cross$geno[[i]]$data
      for(j in unique(cur_map[,1])){
         cur_ori <- population$offspring$genotypes$real[rownames(cur_map)[which(cur_map[,1]==j)],]
         result[i,j]<- mean(apply(cur_ori,1,function(x){apply(cur_new,2,function(y){cor(x,y,use="pair")})}))
      }
      e <- proc.time()
      cat("Counting correlation matrix:",round(i/nchr(cross)*100),"% done estimated time remaining:",((e-s)[3]/i)*(nchr(cross)-i),"s\n")
  }
  invisible(result)
}