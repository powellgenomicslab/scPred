#' @title Eigendecompose gene expression matrix
#' @description Performs principal component analysis from a gene expression matrix. Data is centered and scaled by default.
#' @param expData A matrix object with genes as rows and cells as columns. Unique row names for genes must be provided
#' @param n Number of principal components to be computed
#' @param pseudo Whether to perform a \code{log2(data + 1)} transformation on expression data
#' @param returnData Return training data?
#' @param seed Numeric seed for computing the eigenvalues and eigenvectors using the Lanczos algorithm
#' @return A scPred object with three or four filled slots
#' \itemize{
#' \item \code{svd}: results from \code{prcomp_irlba} function
#' \item \code{expVar}: explained variance by each principal component
#' \item \code{pseudo}: \code{TRUE} if a pseudo-log2 transformation was performed
#' \item \code{trainData}: if \code{returnData} is TRUE, the training data is returned
#' }
#' @keywords singular value decomposition, svd
#' @importFrom methods is
#' @importFrom methods new
#' @export
#' @author José Alquicira Hernández
#' @examples 
#' 
#' # Eigendecompose gene expression matrix
#' 
#' # Simulate gene expression data for two groups
#' 
#' class1 <- matrix(rnbinom(10000, 1, 0.1),  ncol = 100)
#' class2 <- matrix(rnbinom(10000, 1, 0.15),  ncol = 100)
#' 
#' # Create gene expression matrix (rows = cells, colums = genes)
#' 
#' expTrain <- cbind(class1, class2)
#' 
#' # Eigendecompose gene expression matrix
#' 
#' object <- eigenDecompose(expTrain, n = 25)
#' plotEigen(object)
#' 

eigenDecompose <- function(expData, n = 10, pseudo = TRUE, returnData = TRUE, seed = 66){
  
  # Parameter validations
  
  if(!is(expData, "matrix")){
    stop("Expression data must be a matrix object")
  }
  
  expData <- t(expData)
  
  if(pseudo){
    expData <- log2(expData + 1)
  }
  

  # Remove features with zero variance --------------------------------------
  
  zeroVar <- apply(expData, 2, var) == 0
  
  if(any(zeroVar)){
    expData <- expData[,!zeroVar]
    message(paste0(sum(zeroVar), " following genes were removed as their variance is zero across all cells:"))
    cat(paste0(names(zeroVar), collapse = "\n"), "\n", sep = "")
  }

  # Call prcomp() function
  message("Performing Lanczos bidiagonalization...")
  set.seed(66)
  svd <- prcomp_irlba(expData, n = n, center = TRUE, scale. = TRUE)
  class(svd)
  
  rownames(svd$x) <- rownames(expData)
  rownames(svd$rotation) <- colnames(expData)
  
  
  expData <- expData + 0
  nCells <- nrow(expData)
  
 
  f <- function(i) sqrt(sum((expData[, i] - svd$center[i])^2)/(nCells -  1L))
  scale. <- vapply(seq(ncol(expData)), f, pi, USE.NAMES = TRUE)
  
  
  names(scale.) <- names(svd$center)
  svd$scale <- scale.
    
  # Extract variance
  varianceExplained <- svd$sdev**2 / sum(svd$sdev**2)*100
  names(varianceExplained) <- colnames(svd$x)
  
  message("DONE!")
  
  svd <- svd[c("x", "rotation", "center", "scale", "sdev")]
  
  if(returnData){
    return(new("scPred", svd = svd, expVar = varianceExplained, pseudo = pseudo, trainData = t(expData)))
    
  }else{
    return(new("scPred", svd = svd, expVar = varianceExplained, pseudo = pseudo)) 
  }
  
}
