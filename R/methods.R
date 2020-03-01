#' @title show
#' @description Generic display function for \linkS4class{scPred} objects. Displays summary of the
#' object such as number of cells, genes, significant features.
#' @importFrom methods setMethod
#' @importFrom tibble column_to_rownames
#' @importFrom magrittr set_colnames
#' @export

setMethod("show", signature("scPred"), function(object) {
    
    cat("'scPred' object\n")
    
    if(length(object@pVar) != 0){
        cat(sprintf("      Prediction variable = %s\n", object@pVar))
    }
    
    
    
    if(length(object@features) != 0){
        
        cat("\n- Informative PCs per class\n")
        object@features %>%
            sapply(nrow) %>%
            as.data.frame() %>%
            `colnames<-`("Features") %>%
            print()
        
    }
    
    if(length(object@train) != 0){
        
        cat("\n- Training information\n")
        cat(sprintf("      Model: %s\n", object@train[[1]]$modelInfo$label))
        
        
        getMetrics <- function(x, metric){
            
            bestModelIndex <- as.integer(rownames(x$bestTune))
            
            if(metric == "ROC"){
                round(x$results[bestModelIndex, c("ROC", "Sens", "Spec")], 3)
            }else if(metric == "Accuracy"){
                round(x$results[bestModelIndex, c("Accuracy", "Kappa")], 3)
            }else if(metric == "AUC"){
                round(x$results[bestModelIndex, c("AUC", "Precision", "Recall", "F")], 3)
            }
        }
        
        metric <- object@train[[1]]$metric 
        print(t(sapply(object@train, getMetrics, metric)))
        
    }
    
    
})


#' @title Get training probabilities
#' @description Gets training probabilities for each trained model
#' @importFrom methods setMethod
#' @export

setGeneric("getTrainResults", def = function(object) {
    standardGeneric("getTrainResults")
})

#' @title Get training probabilities
#' @description Gets training probabilities for each trained model
#' @importFrom methods setMethod
#' @export

setMethod("getTrainResults", signature("Seurat"), function(object){
    
    if(length(object@misc$scPred@train) == 0){
        stop("No models have been trained")
    }
    
    
    if(ncol(object@misc$scPred@train[[1]]$pred) == 0){
        stop('No training results were calculated. Set savePredictions = "final" and returnData = TRUE')
    }
    
    probs <- lapply(names(object@misc$scPred@train), function(model) extractProb(object@misc$scPred@train[model]))
    names(probs) <- names(object@misc$scPred@train)
    probs
    
})


#' @title Gets contingency table
#' @description Creates a cross table using two columns from the metadata
#' @param object \code{Seurat} object
#' @param true Column name in \code{meta.data} slot that corresponds to the true known classes
#' @param pred Column name in \code{meta.data} slot that corresponds to the predicted classes
#' if they  have been assigned independently from the \code{scPredict()} function
#' @param output Return counts, fraction, or proportions? Default: counts
#' @param digits If proportions are returned, number of digits to round numbers
#' @return A contingency table
#' @export
#' @importFrom dplyr group_by_ summarise
#' @importFrom tidyr spread
#' @importFrom magrittr "%>%"
#' @importFrom tibble column_to_rownames
#' @author Jose Alquicira Hernandez
#'


setGeneric("crossTab", def = function(object, true = NULL, pred = NULL, output = c("counts", "fraction", "prop"), digits = 2) {
    standardGeneric("crossTab")
})

#' @title Gets contingency table
#' @description Creates a cross table using two columns from the metadata
#' @param object \code{Seurat} object
#' @param true Column name in \code{meta.data} slot that corresponds to the true known classes
#' @param pred Column name in \code{meta.data} slot that corresponds to the predicted classes
#' if they  have been assigned independently from the \code{scPredict()} function
#' @param output Return counts, fraction, or proportions? Default: counts
#' @param digits If proportions are returned, number of digits to round numbers
#' @return A contingency table
#' @export
#' @importFrom dplyr group_by_ summarise
#' @importFrom tidyr spread
#' @importFrom magrittr "%>%"
#' @importFrom tibble column_to_rownames
#' @author Jose Alquicira Hernandez
#'

setMethod("crossTab", signature("Seurat"),
          function(object, true = NULL, pred = NULL, output = c("counts", "fraction", "prop"), digits = 2){
              
              if(is.null(true) | is.null(pred)) stop("Provide two column names")
              
              if(!true %in% names(object@meta.data)) stop("'true' column is not present in metadata")
              if(!pred %in% names(object@meta.data)) stop("'pred' column is not present in metadata")
              
              output <- match.arg(output)
              
              cont <- table(object[[pred, drop = TRUE]], object[[true, drop = TRUE]])
              dim_names <- dimnames(cont)
              cont <- as.matrix.data.frame(cont)
              dimnames(cont) <- dim_names
              cont <- as.data.frame(cont)
              sums <- colSums(cont)
              
              if(output == "prop"){
                  mapply(function(x, y) x/y, cont, sums) %>% 
                      as.data.frame() %>% 
                      round(digits) %>% 
                      `dimnames<-`(dim_names) -> cont
              }else if(output == "fraction"){
                  mapply(function(x, y) paste0(x, "/", y), cont, sums) %>% 
                      as.data.frame() %>% 
                      `dimnames<-`(dim_names) -> cont
              }
              
              cont
          })
