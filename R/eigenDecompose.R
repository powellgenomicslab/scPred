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

eigenDecompose <- function(expData, n = 10, pseudo = TRUE){
  
  # Parameter validations
  
  if(!is(expData, "matrix")){
    stop("Expression data must be a matrix")
  }
  
  
  if(pseudo){
    expData <- log2(expData + 1)
  }else{
    expData <- expData
  }
  

  # Remove features with zero variance --------------------------------------
  
  zeroVar <- which(apply(expData, 2, var) == 0)
  
  if(any(zeroVar)){
    expData <- expData[,-zeroVar]
    message("The following genes were removed as their variance is zero across all cells:")
    cat(paste0(names(zeroVar), collapse = "\n"), "\n", sep = "")
  }

  # Call prcomp() function
  message("Performing Lanczos bidiagonalization...")
  pca <- prcomp_irlba(expData, n = n, center = TRUE, scale. = TRUE)
  
  rownames(pca$x) <- rownames(expData)
  rownames(pca$rotation) <- colnames(expData)
  
  expData <- expData + 0
  f <- function(i) sqrt(sum((expData[, i] - pca$center[i])^2)/(nrow(expData) -  1L))
  scale. <- vapply(seq(ncol(expData)), f, pi, USE.NAMES = TRUE)
  names(scale.) <- names(pca$center)
  pca$scale <- scale.
    
  # Extract variance
  varianceExplained <- pca$sdev**2 / sum(pca$sdev**2)*100
  names(varianceExplained) <- colnames(pca$x)
  
  message("DONE!")
  
  return(new("scPred", prcomp = pca, expVar = varianceExplained, pseudo = pseudo)) 
  
  
}
