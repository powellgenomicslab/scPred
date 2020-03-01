#' @title Train a prediction model
#' @description Trains a prediction model from an \code{scPred} object stored in a \code{Seurat} object
#' @param object An \code{Seurat} object with informative PCs obtained using
#' the \code{getFeatureSpace} function
#' @param model Classification model supported via \code{caret} package. A list of all models can be found here:
#' https://topepo.github.io/caret/available-models.html
#' Default: support vector machine with polynomial kernel
#' @param resampleMethod Resample model used in \code{trainControl} function from \code{caret}. 
#' Default: K-fold cross validation
#' @param number Number of iterations for resample method. See \code{trainControl} function
#' @param seed Numeric seed for resample method. Fixed to ensure reproducibility
#' @param metric Performance metric to be used to select best model: `ROC` (area under the ROC curve), 
#' `PR` (area under the precision-recall curve), `Accuracy`, and `Kappa`
#' @param imbalance Proportion of cell type composition to be considered as an imbalance issue. By default,
#' if a cell type is present only in 0.1 (10%) or less from the whole population, scPred will attempt to 
#' train a weighted SVM to account for class imbalance. Set to `1` to avoid this step
#' @param returnData If \code{TRUE}, training data is returned within \code{scPred} object. 
#' @param savePredictions Specifies the set of hold-out predictions for each resample that should be
#' returned. Values can be either "all", "final", or "none".
#' @return A list of \code{train} objects for each cell class (e.g. cell type). See \code{train} function for details.
#' @keywords train, model
#' @importFrom methods is
#' @importFrom pbapply pblapply
#' @export
#' @author
#' Jose Alquicira Hernandez
#' @examples
#'
#' # Train a SVM with a Radial kernel
#' ## A numeric seed is provided for the K-fold cross validation
#' ## The metric ROC is used to select the best tuned model. "Accuracy" and "Kappa" may be used too.
#'
#' pbmc_small <- getFeatureSpace(pbmc_small, "RNA_snn_res.0.8")
#' pbmc_small <- trainModel(object = pbmc_small,
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
    if(!is(object, "Seurat")){
        stop("object must be 'Seurat'")
    }
    
    if(is.null(object@misc$scPred)){
        stop("No features have been determined. Use 'getFeatureSpace()' function")
    }
    
    classes <- names(object@misc$scPred@features)
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
    
    
    object@misc$scPred@train <- modelsRes
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
    
    
    
    
    if(nrow(object@misc$scPred@features[[positiveClass]]) == 0){
        message("No informative principal components were identified for class: ", positiveClass)
    }
    
    namesPC <- as.character(object@misc$scPred@features[[positiveClass]]$PC)
    features <- subsetMatrix(Embeddings(object, reduction = "pca"), namesPC)
    response <-  as.character(object[[object@misc$scPred@pVar, drop = TRUE]])
    
    
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




