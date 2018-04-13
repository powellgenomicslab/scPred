#' @title Predict classes of a new dataset using a trained model
#' @description Predicts cell classes for a new dataset based on a training model, a reference \code{eigenPred} object
#' @param trainData An \code{scPred} object with metadata and informative features obtained.
#' @param referenceData A matrix object with cells as rows and genes (loci) as columns obtained with \code{projectNewData} function
#' @param trainedModel A \code{train} object returned by \code{trainModel} function
#' @keywords prediction, new, test, validation
#' @importFrom methods is
#' @export
#' @author
#' José Alquicira Hernández



eigenPredict <- function(object, newData, threshold = 0.7){
  
  if(!is(object, "scPred")){
    stop("'eigenDec' object must be of class 'eigenDec'")
  }
  
  if(!(is(newData, "matrix") | is(newData, "data.frame"))){
    stop("'predData' object must be a dataframe")
  }
  
  if(length(object@features) == 0){
    stop("No informative principal components have been obtained yet.\nSee getInformativePCs() function")
  }
  
  projection <- projectNewData(newData = newData,
                                referenceData = object)
  
  
  classes <- names(object@features)
  
  res <- lapply(classes, .predictClass, object, projection)
  names(res) <- levels(classes)
  res <- as.data.frame(res)
  row.names(res) <- rownames(projection)
  i <- apply(res, 1, which.max)
  
  prob <- c()
  for(j in seq_len(nrow(res))){
    prob[j] <- res[j,i[j]]
  }
  
  predictions <- data.frame(probability = prob, prePrediction = names(res)[i])
  rownames(predictions) <- rownames(projection)
  
  predictions %>% 
    rownames_to_column("id") %>% 
    mutate(class = ifelse(probability > threshold, as.character(prePrediction), "Unassigned")) %>%
    select(-prePrediction) %>% 
    column_to_rownames("id") -> finalPrediction
  
  finalPrediction
  

}


.predictClass <- function(positiveClass, object, projection){
  # Get features for positive class
  featureList <- object@features[[positiveClass]]
  features <- projection[,as.character(featureList$PC)]
  
  # Perform presictions
  prediction <- predict(object@train[[positiveClass]], 
                        newdata = features, 
                        type = "prob")
  
  prediction[,1, drop = FALSE]
  
}
