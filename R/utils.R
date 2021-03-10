.make_names <- function(x){
  x <- gsub("\\+", "_plus", x)
  x <- gsub("\\-", "_minus", x)
  x <- make.names(x)
}


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


getPalette <- function(n){
  
  if(n < 6){
    c("#29BF12", "#00A5CF", "#DE1A1A", "#574AE2", "#FFBF00")
  }else if(n < 9){
    c("#558aa6", "#B1740F", "#D5006A", "#08585A", "#FFFD98", "#9449d2", "#BBBE64", "#D7263D")
  }else if(n < 13){
    c("#943CB4", "#194D44", "yellow", "#5B6DC8", "#3CA437", "#6B244C", "#6ACDC5", "#DE1A1A", "#BBB53E", "#2A297A", "#995533", "#D590DA")
  }else{
    stop("Too many classes")
  }
  
}