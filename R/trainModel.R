#' @title Train a prediction model
#' @description Trains a prediction model from an \code{scPred} object
#' @param object An \code{scPred} object with informative PCs obtained using
#' the \code{getInformativePCs} function
#' @param model Classification model supported via \code{caret} package. A list of all models can be found here:
#' https://topepo.github.io/caret/available-models.html
#' Default: support vector machine with polynomial kernel
#' @param resampleMethod Resample model used in \code{trainControl} function. Default: K-fold cross validation
#' @param number Number of iterations for resample method. See \code{trainControl} function
#' @param seed Numeric seed for resample method
#' @param metric Performance metric to be used to select best model: `ROC` (area under the ROC curve), 
#' `PR` (area under the precision-recall curve), `Accuracy`, and `Kappa`
#' @param imbalance Proportion of cell type compositions to be considered as an imbalance issue. By default,
#' if a cell type is only present in 0.1 (10%) or less from the whole population, scPred will attempt to 
#' train a weighted SVM to account for class imbalance. Set to `1` to avoid this step
#' @param returnData If \code{TRUE}, training data is returned
#' @param savePredictions an indicator of how much of the hold-out predictions for each resample should be
#' saved. Values can be either "all", "final", or "none". A logical value can also be used that convert to
#'  "all" (for true) or "none" (for false). "final" saves the predictions for the optimal tuning parameters.
#' @return A list of  \code{train} objects for each cell class (e.g. cell type). See \code{train} function for details.
#' @keywords train, model
#' @importFrom methods is
#' @importFrom pbapply pblapply
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
                       number = 5,
                       seed = 66,
                       metric = c("ROC", "PR", "Accuracy", "Kappa"),
                       imbalance = 0.1,
                       returnData = FALSE,
                       savePredictions = "final",
                       allowParallel = FALSE){
    
    
    # Validations -------------------------------------------------------------
    
    # Check class
    if(!any(is(object, "scPred") | is(object, "Seurat"))){
        stop("object must be 'scPred' or 'Seurat'")
    }
    
    if(is(object, "scPred")){
        
        # Check metadata
        if(nrow(object@metadata) == 0){
            stop("No metadata has been assigned to object")
        }
        
        # Validate if features have been determined
        if(length(object@features) == 0){
            stop("No features have been determined. Use 'getFeatureSpace()' function")
        }
        
        
        classes <- names(object@features)
        
        if(is.null(classes)){
            stop("Prediction variable is not contained in metadata")
        }
        
    }else{
        
        if(is.null(object@misc$scPred)){
            stop("No features have been determined. Use 'getFeatureSpace()' function")
        }
        classes <- names(object@misc$scPred$features)
        
    }
    
    metric <- match.arg(metric)
    
    
    # Train a prediction model for each class
    
    if(length(classes) == 2){
        modelsRes <-  .trainModel(classes[1],
                                  object,
                                  model,
                                  resampleMethod,
                                  seed,
                                  metric,
                                  imbalance,
                                  number,
                                  returnData,
                                  savePredictions,
                                  allowParallel)
        modelsRes <- list(modelsRes)
        names(modelsRes) <- classes[1]
        
        
    }else{
        modelsRes <- pblapply(classes, .trainModel,
                              object,
                              model,
                              resampleMethod,
                              seed,
                              metric,
                              imbalance,
                              number,
                              returnData,
                              savePredictions,
                              allowParallel)
        names(modelsRes) <- classes
    }
    
    if(inherits(object, "scPred")){
        object@train <- modelsRes
    }else{
        object@misc$scPred$train <- modelsRes
    }
    object
}

.trainModel <- function(positiveClass,
                        object,
                        model,
                        resampleMethod,
                        seed,
                        metric,
                        imbalance,
                        number,
                        returnData,
                        savePredictions, 
                        allowParallel){
    
    
    if(is(object, "scPred")){
        
        if(nrow(object@features[[positiveClass]]) == 0){
            message("No informative principal components were identified for class: ", positiveClass)
        }
        
        namesPC <- as.character(object@features[[positiveClass]]$PC)
        features <- subsetMatrix(getPCA(object), namesPC)
        
        
        # Get and refactor response variable according to positive class
        # According to twoClassSummary() documentation
        ## "If assumes that the first level of the factor variables corresponds to a relevant result
        ## but the lev argument can be used to change this."
        response <-  as.character(object@metadata[[object@pVar]])

        
    }else{
        
        if(nrow(object@misc$scPred$features[[positiveClass]]) == 0){
            message("No informative principal components were identified for class: ", positiveClass)
        }
        
        namesPC <- as.character(object@misc$scPred$features[[positiveClass]]$PC)
        features <- subsetMatrix(Embeddings(object), namesPC)
        response <-  as.character(object[[object@misc$scPred$pVar, drop = TRUE]])
        
    }
    
    classWeight <- table(response)/ length(response)
    
    if(any(classWeight < imbalance)){
        metric <- "PR"
        model <- "svmRadialWeights"
    }
    
    
    i <- response != positiveClass
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
                               allowParallel = allowParallel)
        
    }else if(metric == "PR"){
        trCtrl <- trainControl(classProbs = TRUE,
                               method = resampleMethod,
                               number = number,
                               summaryFunction = prSummary,
                               returnData = returnData,
                               savePredictions = savePredictions,
                               allowParallel = allowParallel)
        metric <- "AUC"
    }else{
        trCtrl <- trainControl(classProbs = TRUE,
                               method = resampleMethod,
                               number = number,
                               returnData = returnData,
                               savePredictions = savePredictions,
                               allowParallel = allowParallel)
    }
    
    
    if(metric == "AUC"){
    fit <- train(x = features,
                 y = response,
                 method = model,
                 metric = metric,
                 trControl = trCtrl,
                 tuneGrid = expand.grid(sigma = .05,
                                        C = c(.25, .5, 1),
                                        Weight = 1:2))
    }else{
        fit <- train(x = features,
                     y = response,
                     method = model,
                     metric = metric,
                     trControl = trCtrl)
    }
    
    fit
}




