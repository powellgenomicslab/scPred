#' @title Predict cell classes from a new dataset using a trained model
#' @description Predicts cell classes for a new dataset based on a training model, a reference \code{eigenPred} object
#' @param object An \code{scPred} object with metadata and informative features and trained model(s).
#' @param newData A matrix object with cells as rows and genes (loci) as columns
#' @param threshold Threshold used for probabilities to classify cells into classes
#' @param returnProj Set to TRUE to return computed projection
#' @param returnData Returns prediction data (\code{newData})
#' @param informative Set to TRUE to project only informative components
#' @param useProj If a projection matrix is already stored in the \code{scPred} object, perform predictions using this matrix as input
#' @return A data frame with prediction probabilities associated to each class and a \code{predClass} column, 
#' indicating the classification based on the provided threshold
#' @keywords prediction, new, test, validation
#' @importFrom methods is
#' @export
#' @author
#' José Alquicira Hernández
#' @examples 
#' 
#' # Get class probabilities for cells in a new dataset. A 'predClass' columns with 
#' cell classes is returned depending on the threshold parameter.
#' 
#' prediction <- scPredict(object = object, newData = expTest, threshold = 0.9)
#' 



scPredict <- function(object, newData = NULL, threshold = 0.9, 
                      returnProj = TRUE, returnData = FALSE, informative = TRUE,
                      useProj = FALSE){
  
  if(!is(object, "scPred")){
    stop("'object' must be of class 'scPred'")
  }
   
  if(is.null(newData) & (nrow(object@projection) == 0)){ # Neither newData nor projection
    
    stop("No newData or pre-computed projection")
    
  }else if(is.null(newData) & nrow(object@projection)){ # No newData and projection
    
    message("Using projection stored in object as prediction set")
    useProj <- TRUE
    
  }else if(!is.null(newData) & nrow(object@projection)){ # NewData and projection
    
    if(!is(newData, "matrix")){
      stop("'predData' object must be a matrix")
    }
    message("newData provided and projection stored in scPred object. Set 'useProj = TRUE' to override default projection execution")
    
  }
  
  if(length(object@features) == 0){
    stop("No informative principal components have been obtained yet.\nSee getInformativePCs() function")
  }
  
  if(!useProj){
    projection <- projectNewData(object = object,
                                 newData = newData,
                                 informative = informative, 
                                 seurat = if(!is.null(object@svd$seurat)){TRUE}else{FALSE})
  }else{
    projection <- object@projection
  }
  
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
  
  object@predictions <- finalPrediction
  
  if(returnProj & !useProj){
    object@projection <- projection
  }
  
  
  if(returnData & !is.null(newData)){
    if(object@pseudo){
      object@predData <- log2(newData + 1)
    }else{
      object@predData <- newData
    }
  }
  
  return(object)
  
  
}


.predictClass <- function(positiveClass, object, projection){
  # Get features for positive class
  featureList <- object@features[[positiveClass]]
  features <- subsetMatrix(projection, as.character(featureList$PC))

  # Perform presictions
  prediction <- predict(object@train[[positiveClass]], 
                        newdata = features, 
                        type = "prob")
  
  prediction[ , 1, drop = FALSE]
  
}
