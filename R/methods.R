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
    
    if(length(object@pvar) != 0){
        pVar <- object@pvar
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


setOldClass("train")


#' @title Calculate ROC attributes
#' @description Gets ROC object from a cell type
#' @param x Train object generated by caret + scPred
#' @return ROC object
#' @importFrom magrittr "%>%"
#' @importFrom pROC roc
#' @export

setGeneric("calculateROC", def = function(x) {
    standardGeneric("calculateROC")
})


#' @title Calculate ROC attributes
#' @description Gets ROC object from a cell type
#' @param x Train object generated by caret + scPred
#' @return ROC object
#' @importFrom magrittr "%>%"
#' @importFrom pROC roc
#' @export

setMethod("calculateROC", 
          signature("train"),
          function(x){
              data <- x$pred
              levels <- data[, "obs", drop = TRUE] %>% as.character() %>% unique()
              positiveClass <- levels[levels != "other"]
              roc_obj <- roc(data$obs, data[, positiveClass], quiet = TRUE, levels = c("other", positiveClass))
              roc_obj
          })


#' @title Get ROC objects for all cell types
#' @description Gets ROC objects for all cell types from a Seurat/scPred object
#' @param object Seurat object
#' @return Prints scPred object

setGeneric("getROC", def = function(object) {
    standardGeneric("getROC")
})


#' @title Get ROC objects for all cell types
#' @description Gets ROC objects for all cell types from a Seurat/scPred object
#' @param object Seurat object
#' @return Prints scPred object

setMethod("getROC", 
          signature("Seurat"),
          function(object){
              lapply(object@misc$scPred@train, calculateROC)
          })



#' @title Plot training probabilities
#' @description Plots all training probabilities for each cell type
#' @param object Seurat object
#' @param cost The relative cost of a false negative classification (as compared with a false positive classification)
#' @param class_weight Use the proportion of each cell type vs others as the prevalance
#' (n.cases/(n.controls+n.cases)).
#' @param rug Add rug geom to plot
#' @export
#' @importFrom  pROC coords
#' @importFrom magrittr "%>%"
#' @return Plot with the probability distribution for each cell type


setGeneric("plotProbs", def = function(object, cost = 1, class_weight = TRUE, rug = FALSE) {
    standardGeneric("plotProbs")
})




#' @title Plot training probabilities
#' @description Plots all training probabilities for each cell type
#' @param object Seurat object
#' @param cost The relative cost of a false negative classification (as compared with a false positive classification)
#' @param class_weight Use the proportion of each cell type vs others as the prevalance
#' (n.cases/(n.controls+n.cases)).
#' @param rug Add rug geom to plot
#' @export
#' @importFrom  pROC coords
#' @importFrom magrittr "%>%"
#' @return Plot with the probability distribution for each cell type

setMethod("plotProbs", 
          signature("Seurat"),
          function(object, cost = 1, class_weight = TRUE, rug = FALSE){
              
              
              roc_all <- getROC(object)
              
              plot_train <- function(roc_obj, cell_type, cost, class_weight, rug){
                  if(class_weight){
                      weights <- table(roc_obj$original.response) / length(roc_obj$original.response)
                      positiveWeight <- weights[1]
                  }else{
                      p
                  }
                  cell_type <- gsub("\\.", " ", cell_type)
                  x <- coords(roc_obj, "best", "threshold", transpose = TRUE, best.weights = c(cost, positiveWeight))[1]
                  dat <- roc_obj[c("original.predictor", "original.response")] %>% as.data.frame()
                  
                  ggplot(dat, aes(original.predictor, fill = original.response)) +
                      #geom_density(alpha = 0.7) +
                      geom_histogram(alpha = 0.7, bins = 30) +
                      geom_vline(xintercept = x, lty = "dotted", color = "red") +
                      scale_fill_manual(values = c("#472F7D", "#36B779")) +
                      scale_color_manual(values = c("#472F7D", "#36B779")) +
                      xlab(paste0("P(", cell_type, ")")) +
                      ggtitle(cell_type) +
                      theme_classic() +
                      theme(legend.position = "none", legend.title = element_blank()) +
                      theme(plot.title = element_text(hjust = 0.5, face = "bold")) -> p
                  
                  if(rug) p <- p + geom_rug(aes(color = original.response), alpha = 0.7)
                  
                  p
                  
                  
                  
              }
              
              
              mapply(plot_train, roc_all, names(roc_all), cost = cost, class_weight = class_weight, rug = rug, SIMPLIFY = FALSE) %>% 
                  Reduce(`+`, .)
              
              
          }
)


#' @title Plot ROC curve
#' @description Plots a ROC curve for each cell type using the training probabilities
#' @param object Seurat object
#' @param cost The relative cost of a false negative classification (as compared with a false positive classification)
#' @param class_weight Use the proportion of each cell type vs others as the prevalance
#' (n.cases/(n.controls+n.cases)).
#' @importFrom  pROC coords
#' @importFrom magrittr "%>%"
#' @export
#' @return Plot with the ROC curve for each cell type

setGeneric("plotROC", def = function(object, cost = 1, class_weight = TRUE, rug = FALSE) {
    standardGeneric("plotROC")
})


#' @title Plot ROC curve
#' @description Plots a ROC curve for each cell type using the training probabilities
#' @param object Seurat object
#' @param cost The relative cost of a false negative classification (as compared with a false positive classification)
#' @param class_weight Use the proportion of each cell type vs others as the prevalance
#' (n.cases/(n.controls+n.cases)).
#' @importFrom  pROC coords ggroc
#' @importFrom magrittr "%>%"
#' @export
#' @return Plot with the ROC curve for each cell type

setMethod("plotROC", 
          signature("Seurat"),
          function(object, cost = 1, class_weight = TRUE){
              
              roc_all <- getROC(object)
              
              plot_train <- function(roc_obj, cell_type, class_weight, cost){
                  if(class_weight){
                      weights <- table(roc_obj$original.response) / length(roc_obj$original.response)
                      positiveWeight <- weights[1]
                  }else{
                      positiveWeight <- 0.5
                  }
                  cell_type <- gsub("\\.", " ", cell_type)
                  x <- coords(roc_obj, "best", "threshold", transpose = TRUE, best.weights = c(cost, positiveWeight))[1]
                  dat <- roc_obj[c("original.predictor", "original.response")] %>% as.data.frame()
                  
                  p2 <- ggroc(roc_obj, color = "darkblue")
                  thr <- p2$data[p2$data$threshold == x,]
                  p2 +
                      ggtitle(cell_type) +
                      geom_point(aes(specificity, sensitivity), size = 0.1) +
                      geom_abline(intercept = c(1, 0), slope = 1, color = "red", lty = "dotted") +
                      geom_point(aes(specificity, sensitivity), color = "darkred", data = thr) +
                      theme_classic()
                  
                  
              }
              
              mapply(plot_train, roc_all, names(roc_all), cost = cost, class_weight = class_weight, SIMPLIFY = FALSE) %>% 
                  Reduce(`+`, .)
              
          }
)



#' @title Get probability thresholds
#' @description Calculates probability thresholds based on median absolute deviations
#' @param object Seurat object
#' @param tolerance Probability tolerance value for threshold. This value is substracted or added from
#' the computed probability thresholds for the cell type of interest and the remaining cell types respectively.
#' Useful when applying using a new dataset where computed probabilities may be slightly lower
#' @param nmads Number of median absolute deviations to use
#' @importFrom magrittr "%>%"
#' @export
#' @return A data frame of 2 x n cell types. Positive thresholds refer to the cell type of interest and
#' negative to other cell types in a binary classification setting 


setGeneric("getThresholds", def = function(object, tolerance = 0.05, nmads = 3) {
    standardGeneric("getThresholds")
})


#' @title Get probability thresholds
#' @description Calculates probability thresholds based on median absolute deviations
#' @param object Seurat object
#' @param tolerance Probability tolerance value for threshold. This value is substracted or added from
#' the computed probability thresholds for the cell type of interest and the remaining cell types respectively.
#' Useful when applying using a new dataset where computed probabilities may be slightly lower
#' @param nmads Number of median absolute deviations to use
#' @importFrom magrittr "%>%"
#' @export
#' @return A data frame of 2 x n cell types. Positive thresholds refer to the cell type of interest and
#' negative to other cell types in a binary classification setting 

setMethod("getThresholds", 
          signature("Seurat"), 
          function(object, tolerance = 0.05, nmads = 3){
              
              if(!"scPred" %in% names(object@misc)) stop("No feature space has been determined!")
              if(!length(object@misc$scPred@train)) stop("No models have been trained!")
              if(tolerance > 0.5) stop("Tolerance must be lower than 0.5")
              
              findThreshold <- function(values, 
                                        nmads, 
                                        type = c("lower", "upper"), 
                                        na.rm = FALSE) {
                  med_val <- stats::median(values, na.rm = na.rm)
                  mad_val <- stats::mad(values, center = med_val, na.rm = na.rm)
                  if (type == "lower"){
                      limit <- med_val - nmads * mad_val
                  } else if (type == "upper") {
                      limit <- med_val + nmads * mad_val
                  }
                  return(limit)
              }
              
              
              cell_types <- names(object@misc$scPred@train)
              
              calculate_threshold <- function(x){
                  object@misc$scPred@train[[x]]$pred %>% 
                      split(.$obs) -> classes 
                  
                  other <- classes$other[[x]]
                  positiveClass <- classes[[x]][[x]]
                  
                  i <- findThreshold(positiveClass, type = "lower", nmads = nmads)
                  j <- findThreshold(other, type = "upper", nmads = nmads)
                  
                  c(positive = i, negative = j)
              }
              
              
              if(length(cell_types) == 1){
                  classes <- unique(object@meta.data[, "scPred_response"]) %>% as.character()
                  negativeClass <- classes[classes != cell_types]
                  res <- calculate_threshold(cell_types) %>%
                      t() %>% 
                      as.data.frame()
                  res <- rbind(res, rev(as.numeric(res[1,])))
                  rownames(res) <- c(cell_types, negativeClass)
              }else{
                  res <- lapply(cell_types, calculate_threshold) %>% 
                      Reduce(rbind, .) %>% 
                      as.data.frame() %>% 
                      `rownames<-`(cell_types)
              }
              
              if(tolerance){
                  res$positive <- res$positive - tolerance
                  res$negative <- res$negative + tolerance
              }
              
              res %>% 
                  apply(1, function(x){
                      if(x[2] > x[1]){
                          rep(mean(x), 2)
                      }else{
                          x
                      }
                  }) %>% 
                  as.data.frame() %>% 
                  `rownames<-`(c("positive", "negative"))
              
          })


#' @title Get training probabilities
#' @description Gets training probabilities for each cell type
#' @param object Seurat object
#' @importFrom magrittr "%>%"
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr full_join
#' @export
#' @return A data frame with all cell-type probabilities associated to each cell

setGeneric("getProbabilities", def = function(object) {
    standardGeneric("getProbabilities")
})



#' @title Get training probabilities
#' @description Gets training probabilities for each cell type
#' @param object Seurat object
#' @importFrom magrittr "%>%"
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr full_join
#' @export
#' @return A data frame with all cell-type probabilities associated to each cell

setMethod("getProbabilities", 
          signature("Seurat"), 
          function(object){
    
    if(!"scPred" %in% names(object@misc)) stop("No feature space has been determined!")
    if(!length(object@misc$scPred@train)) stop("No models have been trained!")
    
    cell_types <- names(object@misc$scPred@train)
    pVar <- object@misc$scPred@pVar
    
    get_props <- function(x){
        object@misc$scPred@train[[x]]$pred[c(x, "rowIndex")]
    }
    
    probs <- lapply(cell_types, get_props)
    
    mapply(function(x, x_name){
        res <- data.frame(x, barcode = rownames(object@meta.data)[x$rowIndex])
        res$rowIndex <- NULL
        res
    }, probs, cell_types, SIMPLIFY = FALSE)  %>% 
        Reduce(function(x, y) full_join(x, y, by = "barcode"), .) %>%
        column_to_rownames("barcode") -> probs
    
    probs <- probs[match(colnames(object), rownames(probs)), ]
    probs
})
