#' @title Train a prediction model
#' @description Trains a prediction model from an \code{scPred} object
#' @param object An \code{scPred} object with informative PCs obtained using 
#' the \code{getInformativePCs} function
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
#' @importFrom methods is
#' @export
#' @author
#' José Alquicira Hernández



trainModel <- function(object,
                       method = "svmRadial",
                       resampleMethod = "cv",
                       seed = NULL,
                       number = 10,
                       returnData = FALSE,
                       savePredictions = FALSE){
  
  # Validate class
  if(!is(object, "scPred")){
    stop("object must be 'scPred'")
  }
  
  if(nrow(object@metadata) == 0){
    stop("No metadata has been assigned to object")
  }
  
  if(length(object@features) == 0){
    stop("No features have been determined. Use 'getInformativePCs' function")
  }
  
  
  classes <- metadata(object)[[object@pVar]]
  
  if(length(levels(classes)) == 2){
    modelsRes <-  .trainModelByClass(levels(classes)[1],
                                     classes,
                                     object,
                                     method,
                                     resampleMethod,
                                     seed = seed,
                                     number,
                                     returnData,
                                     savePredictions)
    modelsRes <- list(modelsRes)
    names(modelsRes) <- levels(classes)[1]
    
    
  }else{
    modelsRes <- lapply(levels(classes), .trainModelByClass,
                        classes,
                        object,
                        method,
                        resampleMethod,
                        seed,
                        number,
                        returnData,
                        savePredictions)
    names(modelsRes) <- levels(classes)
  }
  
  object@train <- modelsRes
  object
}

.trainModelByClass <- function(positiveClass,
                               classes,
                               object,
                               method,
                               resampleMethod,
                               seed,
                               number,
                               returnData,
                               savePredictions){
  
  if(nrow(object@features[[positiveClass]]) == 0){
    message("No informative principal components were identified for class: ", positiveClass)
  }
  
  
  features <- getPCA(object)[, as.character(object@features[[positiveClass]]$PC)]
  
  
  # Get and refactor response variable according to positive class
  # According to twoClassSummary() documentation
  ## "If assumes that the first level of the factor variables corresponds to a relevant result 
  ## but the lev argument can be used to change this."
  
  i <- classes != positiveClass
  response <- as.character(classes)
  response[i] <- "other"
  response <- factor(response, levels = c(positiveClass, "other"))
  
  
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
  fit
}





