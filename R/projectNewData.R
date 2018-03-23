#' @title Project new data into training principal components
#' @description Projects a new dataset into the principal components obtained from a training dataset
#' @param newData A matrix object with cells as rows and genes (loci) as columns
#' @param referenceData An \code{eigenDec} object
#' @return A data frame with the projection
#' @keywords test, validation, projection
#' @export
#' @author
#' José Alquicira Hernández


projectNewData <- function(newData, referenceData){
  
  if(!is(newData, "matrix")){
    stop("'newData' must be a matrix")
  }
  
  if(!is(referenceData, "eigenPred")){
    stop("'referenceData' must be of class 'eigenDec'")
  }

  
  if(referenceData@pseudo){
    newData <- log2(newData + 1)
  }else{
    newData <- newData
  }
  
  newDataProj <- predict(referenceData@prcomp, newdata = newData)
  
  as.data.frame(newDataProj)
  
}
