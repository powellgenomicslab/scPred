#' @title Project new data into training principal components
#' @description Projects a new dataset into the principal components obtained from a training dataset
#' @param newData A matrix object with cells as rows and genes (loci) as columns
#' @param referenceData An \code{eigenPred} object
#' @return A data frame with the projection
#' @keywords test, validation, projection
#' @importFrom methods is
#' @export
#' @author
#' José Alquicira Hernández


projectNewData <- function(newData, referenceData){
  
  if(!is(newData, "matrix")){
    stop("'newData' must be a matrix")
  }
  
  if(!is(referenceData, "scPred")){
    stop("'referenceData' must be of class 'scPred'")
  }
  
  if(referenceData@pseudo){
    newData <- log2(newData + 1)
  }
  
  new <- colnames(newData)
  ref <- rownames(getLoadings(referenceData))
  
  
  if(!(all(new %in% ref) & all(ref %in% new))){ # Subset genes if necesary
    newSub <- newData[,new %in% ref] 
    refSub <- getLoadings(referenceData)[ref %in% new, ]
    newSub <- newSub[, match(rownames(refSub), colnames(newSub))]
  }else if(!all(new == ref)){ # Order new data according to loadings matrix
    newSub <- newData[, match(ref, new)]
    refSub <- getLoadings(referenceData)
  }else{ # Use data directly if genes match and are ordered
    newSub <- newData
    refSub <- getLoadings(referenceData)
  }
  
  # Get new centers and scales for gene features
  newCenter <- referenceData@prcomp$center[names(referenceData@prcomp$center) %in% colnames(newSub)]
  newScale <- referenceData@prcomp$scale[names(referenceData@prcomp$scale) %in% colnames(newSub)]
  
  
  # Perform linear transformation
  newDataProj <- scale(newSub, newCenter, newScale) %*% refSub
  
  as.data.frame(newDataProj)
  
}