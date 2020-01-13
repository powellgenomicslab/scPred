#' @title Predict cell classes from a new dataset using a trained model
#' @description Predicts cell classes for a new dataset based on a training model, a reference \code{eigenPred} object
#' @param object An \code{scPred} object with metadata and informative features and trained model(s).
#' @param newData A matrix object with genes (loci) as rows and cells as columns
#' @param threshold Threshold used for probabilities to classify cells into classes
#' @param returnProj Set to TRUE to return computed projection
#' @param returnData Returns prediction data (\code{newData})
#' @param informative Set to TRUE to project only informative components
#' @param useProj If a projection matrix is already stored in the \code{scPred} object, perform predictions using this matrix as input
#' @return A data frame with prediction probabilities associated to each class and a \code{predClass} column, 
#' indicating the classification based on the provided threshold
#' @keywords prediction, new, test, validation
#' @importFrom methods is
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom dplyr mutate select 
#' @importFrom pbapply pblapply
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

  # Function validations ----------------------------------------------------
  
  # Validate if provided object is an scPred object
  if(!is(object, "scPred")){
    stop("'object' must be of class 'scPred'")
  }
  
  # Validate if scPred models have been trained already
  if(!length(object@train)){
    stop("No models have been trained!")
  }
  
  # Predictions are only possible if new data is provided. The new data can be a gene expression
  # matrix or a loading projection that has been computed independently.
  # Validate if new data is provided. If only a projection is found in the @projection slot,
  # this one is used as the prediction/test data
  if(is.null(newData) & (nrow(object@projection) == 0)){ # Neither newData nor projection
    
    stop("No newData or pre-computed projection")
    
  }else if(is.null(newData) & nrow(object@projection)){ # No newData and projection
    
    message("Using projection stored in object as prediction set")
    useProj <- TRUE
    
  }else if(!is.null(newData) & nrow(object@projection)){ # NewData and projection
    
    if (!(is(newData, "matrix") | is(newData, "Matrix"))){
      stop("'newData' object must be a matrix or seurat object")
    }
    message("newData provided and projection stored in scPred object. Set 'useProj = TRUE' to override default projection execution")
  }
  
  # Evaluate if there are identified features stored in the scPred object
  if(!length(object@features)){
    stop("No informative principal components have been obtained yet.\nSee getInformativePCs() function")
  }
  
  # Convert newData to sparse Matrix object
  if(is(newData, "matrix")){
    newData <- Matrix(newData)
  }
  
  # Data projection ---------------------------------------------------------

  # By default, projection of training loadings is automatically done. This option can be override
  # with the `useProj` argument to use an independent projection stored directly in the object.  
  if(!useProj){
    projection <- projectNewData(object = object,
                                 newData = newData,
                                 informative = informative, 
                                 seurat = if(!is.null(object@svd$seurat)){TRUE}else{FALSE})
  }else{
    projection <- object@projection
  }
  
  # Get cell classes used for training
  classes <- names(object@features)
  

  # Cell class prediction ---------------------------------------------------
  
  # For all cell classes which a training models has been trained for, get
  # the prediction probability for each cell in the `projection` object
  
  message("Predicting cell types")
  res <- pblapply(classes, .predictClass, object, projection)
  names(res) <- levels(classes)
  
  # Gather results in a dataframe
  res <- as.data.frame(res)
  row.names(res) <- rownames(projection)
  
  if(length(classes) == 1){
    
    # If only one model was trained (positive and negative classes present in prediction
    # variable only), get the probability for the negative class using the complement rule
    # P(negative_class) = 1 - positive_class
    allClasses <-  unique(object@metadata[[object@pVar]])
    positiveClass <- classes
    negativeClass <- as.character(allClasses[!allClasses %in% classes])
    res[[negativeClass]] <- 1 - res[[positiveClass]]
    
    # Assign cells according to their associated probabilities
    res$predClass <- ifelse(res[,1] > threshold, positiveClass,
                        ifelse(res[,2] > threshold, negativeClass, "unassigned")) %>%
      as.factor()
    
    # Save results to `@predictions` slot
    object@predictions <- res
    
    
  }else{
    
    # If more than one model was trained (there are 3 classes or more in the prediction variable),
    # obtain the maximum probability for each cell and assign the respective class to that cell.
    
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
    
    # Save results to `@predictions` slot
    object@predictions <- finalPrediction
  }
  
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
