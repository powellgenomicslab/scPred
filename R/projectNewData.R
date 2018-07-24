#' @title Project new data onto training principal components
#' @description Projects a new dataset into the principal components obtained from a training dataset
#' @param object An \code{scPred} object
#' @param newData A matrix object with cells as rows and genes (loci) as columns
#' @param informative Perfoms rotation using only informative principal components
#' @return A data frame with the projection
#' @keywords test, validation, projection
#' @importFrom methods is
#' @export
#' @author
#' José Alquicira Hernández
#' @examples 
#' 
#' # Prjects training discriminant principal axes onto test dataset. 
#' ## By settint "informative = FALSE", all principal components
#' ## (including non-informative) are projected.
#' 
#' projection <- projectNewData(object, expTest, informative = FALSE)
#' 


projectNewData <- function(object, newData, informative = TRUE){
  
  if(!is(newData, "matrix")){
    stop("'newData' must be a matrix")
  }
  
  if(!is(object, "scPred")){
    stop("'object' must be of class 'scPred'")
  }
  
  if(object@pseudo){
    newData <- log2(newData + 1)
  }
  
  new <- colnames(newData)
  ref <- rownames(getLoadings(object))
  
  
  if(!(all(new %in% ref) & all(ref %in% new))){ # Subset genes if necesary
    newSub <- newData[,new %in% ref] 
    refSub <- getLoadings(object)[ref %in% new, ]
    newSub <- newSub[, match(rownames(refSub), colnames(newSub))]
  }else if(!all(new == ref)){ # Order new data according to loadings matrix
    newSub <- newData[, match(ref, new)]
    refSub <- getLoadings(object)
  }else{ # Use data directly if genes match and are ordered
    newSub <- newData
    refSub <- getLoadings(object)
  }
  
  # Get new centers and scales for gene features
  newCenter <- object@pca$center[names(object@pca$center) %in% colnames(newSub)]
  newScale <- object@pca$scale[names(object@pca$scale) %in% colnames(newSub)]
  
  if(informative){
  informativePCs  <- object@features %>% 
      lapply("[[", "PC") %>% 
      unlist() %>% 
      as.vector() %>% 
      unique() 
  
  features <- colnames(refSub) %in% informativePCs
  refSub <- refSub[,features]
  }
  
  
  
  # Perform linear transformation
  newDataProj <- scale(newSub, newCenter, newScale) %*% refSub
  
  as.data.frame(newDataProj)
  
}
