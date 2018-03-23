#' @title Predict classes of a new dataset using a trained model
#' @description Predicts cell classes for a new dataset based on a training model, a reference \code{eigenPred} object
#' @param trainData An \code{eigenPred} object with metadata and informative features obtained.
#' @param referenceData A matrix object with cells as rows and genes (loci) as columns obtained with \code{projectNewData} function
#' @param trainedModel A \code{train} object returned by \code{trainModel} function
#' @keywords prediction, new, test, validation
#' @export
#' @author
#' José Alquicira Hernández



eigenPredict <- function(trainData, predData, trainedModel, classes = NULL){
  
  if(!is(trainData, "eigenPred")){
    stop("'eigenDec' object must be of class 'eigenDec'")
  }
  
  if(!is(predData, "data.frame")){
    stop("'predData' object must be a dataframe")
  }
  
  if(!is(trainedModel, "train")){
    stop("'trainedModel' object must be of class 'train'")
  }
  
  
  if(length(trainData@features) == 0){
    stop("No informative principal components have been obtained yet.\nSee getInformativePCs() function")
  }
  
  # Subset informative principal components from prediction dataset
  if(length(trainData@features$PC) == 0){
    stop("No significant principal components were found")
  }
  
  featureList <- trainData@features$PC[seq_len(trainedModel$top)]
  
  
  projection <- predData[,featureList]
  prediction <- predict(trainedModel, newdata = as.matrix(projection), type = "prob")
  
  if(!is.null(classes) & is.numeric(classes)){
    prediction <- factor(ifelse(prediction[,1] > classes, names(prediction)[1], names(prediction)[2]),
           levels = names(prediction))
    prediction <- data.frame(class = prediction)
  }

  row.names(prediction) <- row.names(predData)
  
  prediction
}
