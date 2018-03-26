#' @title Eigendecompose gene expression matrix via SVD
#' @description Performs principal component analysis from a gene expression matrix. Data is centered and scaled by default.
#' @param expdata A matrix object with cells as rows and genes (loci) as columns. Unique row names must be provided
#' @param pseudo Whether to perform a \code{log2(data + 1)} transformation on expression data
#' @param args List with extra arguments provided to \code{prcomp()} function.
#' @param decimals Number of decimals for explained variance results
#' @return A eigenPred object with three filled slots
#' \itemize{
#' \item \code{prcomp}: results from \code{prcomp} function
#' \item \code{expVar}: explained variance by each principal component
#' \item \code{pseudo}: \code{TRUE} if a pseudo-log2 transformation was performed
#' }
#' @keywords singular value decomposition, pca
#' @importFrom methods is
#' @importFrom methods new
#' @export
#' @author José Alquicira Hernández
#' @examples 
#' 
#' # Eigendecompose gene expression matrix
#' res <- eigendecompose(expTrain)

eigenDecompose <- function(expData, pseudo = TRUE, args = NULL, decimals = 2){
  
  # Parameter validations
  
  if(!is(expData, "matrix")){
    stop("Expression data must be a matrix")
  }
  
  
  if(pseudo){
    expData <- log2(expData + 1)
  }else{
    expData <- expData
  }
  
  namesArgs <- names(args)
  
  
  # Call prcomp() function
  
  if(is.null(args)){
    pca <- prcomp(expData, center = TRUE, scale. = TRUE)
  }else{
    args$x <- expData
    pca <- do.call(prcomp, args = args)
  }
  
  # Extract variance
  
  varianceExplained <- round((pca$sdev**2 / sum(pca$sdev**2))*100, decimals)
  names(varianceExplained) <- colnames(pca$x)
  
  
  return(new("eigenPred", prcomp = pca, expVar = varianceExplained, pseudo = pseudo)) 
  
  
}
