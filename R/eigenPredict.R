#' @title Predict classes of a new dataset using a trained model
#' @description Predicts cell classes for a new dataset based on a training model, a reference \code{eigenPred} object
#' @param object An \code{scPred} object with metadata and informative features obtained.
#' @param newData A matrix object with cells as rows and genes (loci) as columns obtained with \code{projectNewData} function
#' @param threshold Threshold used for probabilities to classify cells into classes
#' @return A data frame with prediction probabilities associated to each class and a \code{predClass} column, 
#' indicating the classification based on the provided threshold
#' @keywords prediction, new, test, validation
#' @importFrom methods is
#' @export
#' @author
#' José Alquicira Hernández



scPredict <- function(object, newData, threshold = 0.9){
  
  if(!is(object, "scPred")){
    stop("'object' must be of class 'scPred'")
  }
  
  if(!(is(newData, "matrix") | is(newData, "data.frame"))){
    stop("'predData' object must be a dataframe or a matrix")
  }
  
  if(length(object@features) == 0){
    stop("No informative principal components have been obtained yet.\nSee getInformativePCs() function")
  }
  
  projection <- projectNewData(object = object,
                               newData = newData,
                               informative = TRUE)
  
  
  classes <- names(object@features)
  
  res <- lapply(classes, .predictClass, object, projection)
  names(res) <- levels(classes)
  res <- as.data.frame(res)
  row.names(res) <- rownames(projection)
  
  # if(length(classes) == 1){
  #   classes <- levels(object@metadata[[object@pVar]])
  #   res[[classes[2]]] <- 1 - res[[classes[1]]]
  #   
  #   res$class <- ifelse(res[,1] > threshold, classes[1], 
  #                       ifelse(res[,2] > threshold, classes[2], "unassigned")) %>% 
  #     as.factor()
  #   
  #   res[,2] <- NULL 
  #   names(res) <- c("probability", "class")
  #   return(res)
  # }
  
  i <- apply(res, 1, which.max)
  
  prob <- c()
  for(j in seq_len(nrow(res))){
    prob[j] <- res[j,i[j]]
  }
  
  predictions <- data.frame(res, probability = prob, prePrediction = names(res)[i])
  rownames(predictions) <- rownames(projection)
  
  predictions %>% 
    rownames_to_column("id") %>% 
    mutate(predClass = ifelse(probability > threshold, as.character(prePrediction), "unassigned")) %>%
    select(-probability, -prePrediction) %>% 
    column_to_rownames("id") -> finalPrediction
  
  return(finalPrediction)
  
  
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
