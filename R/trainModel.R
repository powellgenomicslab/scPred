#' @title Train a prediction model
#' @description Trains a prediction model from an \code{scPred} object
#' @param object An \code{scPred} object with informative PCs obtained using 
#' the \code{getInformativePCs} function
#' @param model Classification model supported via \code{caret} package. A list of all models can be found here: 
#' https://topepo.github.io/caret/available-models.html
#' Default: support vector machine with polynomial kernel
#' @param resampleMethod Resample model used in \code{trainControl} function. Default: K-fold cross validation 
#' @param seed Numeric seed for resample model
#' @param number Number of iterations for resample method. See \code{trainControl} function
#' @param returnData If \code{TRUE}, training data is returned
#' @param savePredictions an indicator of how much of the hold-out predictions for each resample should be 
#' saved. Values can be either "all", "final", or "none". A logical value can also be used that convert to
#'  "all" (for true) or "none" (for false). "final" saves the predictions for the optimal tuning parameters.
#' @return A list of  \code{train} objects for each cell class (e.g. cell type). See \code{train} function for details.
#' @keywords train, model
#' @importFrom methods is
#' @export
#' @author
#' José Alquicira Hernández
#' @examples 
#' 
#' # Train a SVM with a Radial kernel
#' ## A numeric seed is provided for the K-fold cross validation
#' ## The metric ROC is used to select the best tuned model. "Accuracy" and "Kappa" may be used too.
#' 
#' object <- trainModel(object = object, 
#'                      seed = 1234, 
#'                      metric = "ROC")
#' 




trainModel <- function(object,
                       model = "svmRadial",
                       resampleMethod = "cv",
                       number = 10,
                       seed = NULL,
                       metric = c("ROC", "Accuracy", "Kappa"),
                       returnData = TRUE,
                       savePredictions = "final"){
  
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
  metric <- match.arg(metric)
  
  if(length(levels(classes)) == 2){
    modelsRes <-  .trainModelByClass(levels(classes)[1],
                                     classes,
                                     object,
                                     model,
                                     resampleMethod,
                                     seed,
                                     metric,
                                     number,
                                     returnData,
                                     savePredictions)
    modelsRes <- list(modelsRes)
    names(modelsRes) <- levels(classes)[1]
    
    
  }else{
    modelsRes <- lapply(levels(classes), .trainModelByClass,
                        classes,
                        object,
                        model,
                        resampleMethod,
                        seed,
                        metric,
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
                               model,
                               resampleMethod,
                               seed,
                               metric,
                               number,
                               returnData,
                               savePredictions){
  
  if(nrow(object@features[[positiveClass]]) == 0){
    message("No informative principal components were identified for class: ", positiveClass)
  }
  
  namesPC <- as.character(object@features[[positiveClass]]$PC)
  features <- subsetMatrix(getPCA(object), namesPC)

  
  # Get and refactor response variable according to positive class
  # According to twoClassSummary() documentation
  ## "If assumes that the first level of the factor variables corresponds to a relevant result 
  ## but the lev argument can be used to change this."
  
  i <- classes != positiveClass
  response <- as.character(classes)
  response[i] <- "other"
  response <- factor(response, levels = c(positiveClass, "other"))
  
  
  if(!is.null(seed)) set.seed(seed)
  
  if(metric == "ROC"){
  trCtrl <- trainControl(classProbs = TRUE,
                         method = resampleMethod,
                         number = number,
                         summaryFunction = twoClassSummary,
                         returnData = returnData,
                         savePredictions = savePredictions,
                         allowParallel = FALSE)
  }else{
    trCtrl <- trainControl(classProbs = TRUE,
                           method = resampleMethod,
                           number = number,
                           returnData = returnData,
                           savePredictions = savePredictions,
                           allowParallel = FALSE)
  }
  
  fit <- train(x = features, 
               y = response, 
               method = model,
               metric = metric,
               trControl = trCtrl)
  fit
}





