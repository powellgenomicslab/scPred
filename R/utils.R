.intersectMat <- function(ref, new){
  
  # Get gene names from loadings reference
  refRows <- rownames(ref)
  # Get gene names from new dataset
  newRows <- rownames(new)
  
  if(!(all(newRows %in% refRows) & all(refRows %in% newRows))){ # Subset genes if necesary
    newSub <- new[newRows %in% refRows, ] 
    refSub <- ref[refRows %in% newRows, ]
    newSub <- newSub[match(rownames(refSub), rownames(newSub)), ]
  }else if(!all(newRows == refRows)){ # Order new data according to loadings matrix
    newSub <- newData[match(refRows, newRows), ]
    refSub <- ref
  }else{ # Use data directly if genes match and are ordered
    newSub <- new
    refSub <- ref
  }
  
  return(list(ref = refSub, new = newSub))  
  
}


subsetMatrix <- function(x, s, by.col = TRUE, drop = FALSE, verbose = FALSE, ...){

  if(by.col){
    ids <- colnames(x, ...)
    if(is.null(ids) & verbose) stop("No colnames were found!")
    i <- ids %in% s
    if(!any(i) & verbose){message("No matches were found")}
    x <- x[,i, drop = drop]
  }else{
    ids <- rownames(x, ...)
    if(is.null(ids) & verbose) stop("No rownames were found!")
    i <- ids %in% s
    if(!any(i) & verbose){message("No matches were found")}
    x <- x[i, , drop = drop]
  }
  
  x
    
}

scaleDataSeurat <- function(
  data,
  genes.use = NULL,
  center,
  scale,
  scale.max = 10
){
  
  if(is.null(genes.use)){
    genes.use <- rownames(x = data)
  }
  
  genes.use <- rownames(x = data)
  genes.use <- intersect(x = genes.use, y = rownames(x = data))
  data.use <-  data[genes.use, ]
  
  scale.data <- matrix(
    data = NA,
    nrow = length(x = genes.use),
    ncol = ncol(x = data)
  )
  
  
  dimnames(x = scale.data) <- dimnames(x = data.use)
  
  bin.size <- 1000
  max.bin <- floor(length(genes.use)/bin.size) + 1
  message("Scaling data matrix")
  pb <- txtProgressBar(min = 0, max = max.bin, style = 3, file = stderr())
  for (i in 1:max.bin) {
    my.inds <- ((bin.size * (i - 1)):(bin.size * i - 1)) + 1
    my.inds <- my.inds[my.inds <= length(x = genes.use)]

    new.data <- t(
      x = scale(
        x = t(x = as.matrix(x = data.use[genes.use[my.inds], ])),
        center = center[genes.use[my.inds]],
        scale = scale[genes.use[my.inds]]
      )
    )
    
    new.data[new.data > scale.max] <- scale.max
    scale.data[genes.use[my.inds], ] <- new.data
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  scale.data[is.na(scale.data)] <- 0
  
  scale.data
}
