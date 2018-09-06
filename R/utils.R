.intersectLoadings <- function(ref, new){
  
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
