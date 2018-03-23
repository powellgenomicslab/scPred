#' @title Train a prediction model
#' @description Trains a prediction model from an \code{eigenPred} object
#' @param object An \code{eigenPred} object with informative PCs obtained using 
#' the \code{getInformativePCs} function
#' @param positiveClass If a positive class is provided, the best model is selected based on the highest ROC value. Otherwise, accuracy is used
#' @param top Top n significant features from \code{features} slot to be used as predictors
#' @param method Classification model supported via \code{caret} package
#' Default: support vector machine with polynomial kernel
#' @param resampleMethod Resample method used in \code{trainControl} function. Default: K-fold cross validation 
#' @param seed Seed to apply the resample method
#' @param number Number of iterations for resample method. See \code{trainControl} function

#' @param returnData If \code{TRUE}, training data is returned
#' @param savePredictions If \code{TRUE}, predictions for training data are returned
#' @return A \code{train} object with final results. See \code{train} function for details. An aditional value \code{top} is added to the 
#' \code{train} object to track the features used to train the model
#' @keywords train, model
#' @export
#' @author
#' José Alquicira Hernández



trainModel <- function(object,
                       positiveClass = NULL,
                       top = 10,
                       method = "svmPoly",
                       resampleMethod = "cv",
                       seed = NULL,
                       number = 10,
                       returnData = FALSE,
                       savePredictions = FALSE){
  
  # Validate class
  if(!is(object, "eigenPred")){
    stop("object must be 'object class'")
  }
  
  if(nrow(object@metadata) == 0){
    stop("No metadata has been assigned to object")
  }
  
  if(nrow(object@features) == 0){
    stop("No features have been determined. Use 'getInformativePCs' function")
  }
  
  
  # Get features
  if(length(object@features$PC) == 0){
      stop("No significant principal components were found")
  }
  
  if(nrow(object@features) < top){
    warning(sprintf("Only % principal components were determined as significant. Using these as features"))
  }  
  features <- getPCA(object)[,object@features$PC[seq_len(top)]]

  # Get response variable
  response <- object@metadata[[object@pVar]]
  
  if(!is.null(positiveClass)){
    # Get and refactor response variable according to positive class
    # According to twoClassSummary() documentation
    ## "If assumes that the first level of the factor variables corresponds to a relevant result 
    ## but the lev argument can be used to change this."
    
    if(!any(levels(response) == positiveClass)){
      stop(paste0("positiveClass '", positiveClass, "' is not included in prediction variable '", object@pVar, "'"))
    }
    orderedLevels <- c(levels(response)[levels(response) == positiveClass], levels(response)[levels(response) != positiveClass])
    response <- factor(response, levels = orderedLevels)
  }
  
  # Train model
  
  if(!is.null(positiveClass)){
    if(!is.null(seed)) set.seed(seed)
    trCtrl <- trainControl(classProbs = TRUE,
                           method = resampleMethod,
                           number = number,
                           summaryFunction = twoClassSummary,
                           returnData = returnData,
                           savePredictions = savePredictions,
                           allowParallel = FALSE)
    
    fit <- train(x = as.matrix(features), 
                 y = response, 
                 method = method,
                 metric = "ROC",
                 trControl = trCtrl)
  }else{
    if(!is.null(seed)) set.seed(seed)
    trCtrl <- trainControl(method = resampleMethod,
                           number = number,
                           savePredictions = savePredictions,
                           allowParallel = FALSE)
    
    fit <- train(x = as.matrix(features), 
                 y = response, 
                 method = method,
                 trControl = trCtrl,
                 returnData = returnData)
    
  }
  
  fit$top <- top
  
  fit
}
