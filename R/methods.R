#' @title show
#' @description Displays scPred object info
#' @param string list with output strings
#' @return Prints scPred object
#' 

print.scPred <- function(string) {
    params <- list(
        crayon::bold(string$title),
        crayon::green(string$pVar),
        crayon::underline(string$pVar_value),
        crayon::green(string$features_section),
        crayon::blue(string$features_header),
        string$features,
        string$training_section
    )
    
    for (p in params) {
        cat(p)
    }
}



#' @title show
#' @description Generic display function for \linkS4class{scPred} objects. Displays summary of the
#' object such as number of cells, genes, significant features.
#' @importFrom methods setMethod
#' @importFrom methods show
#' @importFrom tibble column_to_rownames
#' @importFrom magrittr set_colnames
#' @export

setMethod("show", signature("scPred"), function(object) {
    
    training_value <- NA
    
    if(length(object@pVar) != 0){
        pVar <- object@pVar
    }
    
    if(length(object@features) != 0){
        object@features %>%
            sapply(nrow) -> nFeatures
        features <- mapply(function(x,y){ paste0(x, "\t", y)}, nFeatures, names(nFeatures))
        features <- paste0(paste0(features, collapse = "\n"), "\n")
    }
    
    # if(length(object@train) != 0){
    #     
    #     cat("\n- T\n")
    #     cat(sprintf("      Model: %s\n", object@train[[1]]$modelInfo$label))
    #     
    #     
    #     getMetrics <- function(x, metric){
    #         
    #         bestModelIndex <- as.integer(rownames(x$bestTune))
    #         
    #         if(metric == "ROC"){
    #             round(x$results[bestModelIndex, c("ROC", "Sens", "Spec")], 3)
    #         }else if(metric == "Accuracy"){
    #             round(x$results[bestModelIndex, c("Accuracy", "Kappa")], 3)
    #         }else if(metric == "AUC"){
    #             round(x$results[bestModelIndex, c("AUC", "Precision", "Recall", "F")], 3)
    #         }
    #     }
    #     
    #     metric <- object@train[[1]]$metric 
    #     print(t(sapply(object@train, getMetrics, metric)))
    #     
    # }
    
    
    if(is.na(training_value)){
       training_section <- crayon::red(c(cli::symbol$cross, "Training model(s)\n"))
    }
    
    
    string <- list(
        title = "'scPred' object\n",
        pVar = c(cli::symbol$tick, " Prediction variable = "),
        pVar_value =  c(pVar, "\n"),
        features_section = c(cli::symbol$tick, " Discriminant PCs per cell type\n"),
        features_header ="n\tCell Type\n", 
        features = features,
        training_section = training_section
    )
    print.scPred(string)
    
    
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
