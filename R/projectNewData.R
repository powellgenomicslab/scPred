#' @title Project new data onto training principal components
#' @description Projects a new dataset into the principal components obtained from a training dataset
#' @param object An \code{scPred} object
#' @param newData A matrix object with genes as rows and cells as columns
#' @param informative Perfoms rotation using only informative principal components
#' @param seurat Performs Seurat scaling?
#' @return A data frame with the projection
#' @keywords test, validation, projection
#' @importFrom methods is
#' @export
#' @author
#' José Alquicira Hernández
#' @examples 
#' 
#' # Prjects training discriminant principal axes onto test dataset. 
#' ## By setting "informative = FALSE", all principal components
#' ## (including non-informative) are projected.
#' 
#' projection <- projectNewData(object, expTest, informative = FALSE)
#' 


projectNewData <- function(object, newData, informative = TRUE, seurat = FALSE){
  
  if(!is(newData, "matrix")){
    stop("'newData' must be a matrix")
  }
  
  if(!is(object, "scPred")){
    stop("'object' must be of class 'scPred'")
  }
  
  if(object@pseudo){
    newData <- log2(newData + 1)
  }
  
  
  res <- .intersectMat(ref = getLoadings(object), new =  newData)
  refSub <- res$ref
  newSub <- res$new
  rm(res)

  newCenter <- object@svd$center[names(object@svd$center) %in% rownames(newSub)]
  newScale <- object@svd$scale[names(object@svd$scale) %in% rownames(newSub)]
  
  if(informative){
  informativePCs  <- object@features %>% 
      lapply("[[", "PC") %>% 
      unlist() %>% 
      as.vector() %>% 
      unique() 
  
  features <- colnames(refSub) %in% informativePCs
  refSub <- refSub[,features, drop = FALSE]
  }
  
  
  
  # Scale data
  if(seurat){
    newSubScale <- t(scaleDataSeurat(newSub, 
                                   genes.use = rownames(refSub), 
                                   center = newCenter, 
                                   scale = newScale))
  }else{
    newSubScale <- as.matrix(scale(t(newSub), newCenter, newScale))
  }
  
  # Perform linear transformation
  newDataProj <- newSubScale %*% refSub
  
  newDataProj
  
}
