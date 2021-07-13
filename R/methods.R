#' @title show
#' @description Displays scPred object info
#' @param string list with output strings
#' @return Prints scPred object
#' 

print.scPred <- function(string) {
    params <- list(
        title = crayon::bold("'scPred' object\n"),
        prediction_title = crayon::green(c(cli::symbol$tick, " Prediction variable = ")),
        prediction_variable = crayon::underline(string$pVar_value),
        feature_section = crayon::green(c(cli::symbol$tick, " Discriminant features per cell type\n")),
        summary = string$training_section,
        summary = crayon::blue("Summary\n\n")
    )
    
    for (p in params) {
        cat(p)
    }
    
    cat(string$display, sep = "\n")
    
    
}



#' @title show
#' @description Generic display function for \linkS4class{scPred} objects. Displays summary of the
#' object such as number of cells, genes, significant features.
#' @importFrom methods setMethod
#' @importFrom methods show
#' @importFrom knitr kable
#' @export

setMethod("show", signature("scPred"), function(object) {
    
    # Extract prediction variable
    if(length(object@pvar) != 0){
        pVar <- object@pvar
    }
    
    # Extract number of cells per cell type
    n <- get_metadata(object)[,"pvar"] 
    n <- table(n)
    n <- as.data.frame(n)
    names(n) <- c("Cell type", "n")
    
    # Extract number of features per cell type
    if(length(object@features) != 0){
        features <- sapply(object@features, nrow)
        features <- as.data.frame(features)
        names(features) <- "Features"
        features$`Cell type` <- rownames(features)
        rownames(features) <- NULL
    }
    
    training_value <- NULL
    # Extract performance metrics from each classifier
    if(length(object@train) != 0){
        
        classifiers <- get_classifiers(object)
        
        training_value <- lapply(classifiers, function(x, metric){
            bestModelIndex <- as.integer(rownames(x$bestTune))
            metric <- x$metric
            method <- x$method
            if(metric == "ROC"){
                perf <- round(x$results[bestModelIndex, c("ROC", "Sens", "Spec")], 3)
            }else if(metric == "Accuracy"){
                perf <- round(x$results[bestModelIndex, c("Accuracy", "Kappa")], 3)
            }else if(metric == "AUC"){
                perf <- round(x$results[bestModelIndex, c("AUC", "Precision", "Recall", "F")], 3)
            }
            cbind(Method = method, perf)
        })
        
        training_value <- do.call(rbind, training_value)
        training_value$`Cell type` <- rownames(training_value)
        
    }
    
    
    
    display <- merge(n, features, by = "Cell type")
    
    if(is.null(training_value)){
        training_section <- crayon::red(c(cli::symbol$cross, " Training model(s)\n"))
    }else{
        training_section <- crayon::green(c(cli::symbol$tick, " Training model(s)\n"))
        display <- merge(display, training_value, by = "Cell type")
    }
    
    
    string <- list(
        pVar_value = c(pVar, "\n"),
        training_section = training_section,
        display = kable(display)
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
                  cont <- as.data.frame(mapply(function(x, y) x/y, cont, sums))
                  cont <- round(cont, digits)
                  dimnames(cont) <- dim_names   
                   
              }else if(output == "fraction"){
                  cont <- as.data.frame(mapply(function(x, y) paste0(x, "/", y), cont, sums))
                  dimnames(cont) <- dim_names
              }
              
              cont
          })


setOldClass("train")


#' @title Get scPred object
#' @description Accessor function to retrieve scPred models from Seurat objects
#' @param object Seurat object
#' @return scPred object
#' @export

setMethod("get_scpred", 
          signature("Seurat"),
          function(object){
              scpred <- object@misc$scPred
              if(is.null(scpred))
                  stop("No scPred model is stored in this Seurat object")
              else
                  scpred
          })



#' @title Get classification models
#' @description Accessor function to retrieve scPred models from Seurat objects
#' @param object \code{Seurat} object
#' @return A list of \code{train} objects
#' @export

setMethod("get_classifiers", 
          signature("Seurat"),
          function(object){
              scpred <-  get_scpred(object)
              get_classifiers(scpred)
          })


#' @title Get classification models
#' @description Accessor function to retrieve trained models from scPred objects
#' @param object \code{scPred} object
#' @importFrom  dplyr distinct
#' @return A list of \code{train} objects
#' @export

setMethod("get_classifiers", 
          signature("scPred"),
          function(object){
              models <- object@train
              if(!length(models)) stop("No models have been trained!")
              models
          })



#' @title Plot training probabilities
#' @description Plots all training probabilities for each cell type
#' @param object Seurat or scPred object
#' @param size Point size for each cell
#' @importFrom ggplot2 ggplot aes xlab ylab scale_color_manual facet_wrap theme element_text element_blank element_rect element_line
#' @importFrom ggbeeswarm geom_quasirandom
#' @return Plot with the probability distribution for each cell type

.plot_probabilities <- function(object, size = 0.8){
    
    
    dat <- cbind(get_probabilities(object), 
          get_metadata(object)) 
    
    dat$response <- NULL
    
    dat <- reshape(dat,
            v.names = "value", 
            timevar = "class",
            times = unique(as.character(dat$pvar)), 
            direction = "long", 
            varying = unique(as.character(dat$pvar)))

    
    dat$response <- ifelse(as.character(dat$pvar) == dat$class, "Positive", "Negative")
    dat$response <- factor(dat$response, c("Positive", "Negative"))
    
    
    ggplot(dat, aes(response, value,  color = response)) +
        geom_quasirandom(method = "smiley", size = size) +
        xlab("Response") +
        ylab("Probability") +
        scale_color_manual(values = c("#277C8E", "#50164B")) +
        facet_wrap(~class) +
        theme(text = element_text(size = 14),
              panel.background = element_blank(),
              axis.text.x = element_text(color = "black"),
              axis.text.y = element_text(color = "black"),
              axis.line = element_line(size = 0.25),
              strip.background = element_rect(color = "black", fill =  "#ffe5cc"),
              panel.border = element_rect(fill = NA, colour = "grey20")) +
        theme(legend.position = "none")
    
}

#' @title Plot training probabilities
#' @description Plots all training probabilities for each cell type
#' @param object Seurat object
#' @param size Point size for each cell
#' @return Plot with the probability distribution for each cell type

setMethod("plot_probabilities", 
          signature("Seurat"), .plot_probabilities)

#' @title Plot training probabilities
#' @description Plots all training probabilities for each cell type
#' @param object scPred object
#' @param size Point size for each cell
#' @return Plot with the probability distribution for each cell type

setMethod("plot_probabilities", 
          signature("scPred"), .plot_probabilities)



#' @title Get training probabilities
#' @description Gets training probabilities for each cell type
#' @param object Seurat object
#' @export
#' @return A data frame with all cell-type probabilities associated to each cell

setMethod("get_probabilities", 
          signature("Seurat"), 
          function(object){
              get_probabilities(get_scpred(object))
          })


#' @title Get training probabilities
#' @description Gets training probabilities for each cell type
#' @param object scPred object
#' @export
#' @return A data frame with all cell-type probabilities associated to each cell

setMethod("get_probabilities", 
          signature("scPred"), 
          function(object){
              
              models <- get_classifiers(object)
              
              probs <- lapply(models, function(x){
                  i <- x$levels[x$levels != "other"]
                  x$pred[c(i, "rowIndex")]
                  
              })
              
              barcodes <- rownames(get_metadata(object))
              
              probs <- mapply(function(x, x_name){
                  res <- data.frame(x, barcode = barcodes[x$rowIndex])
                  res$rowIndex <- NULL
                  names(res)[1] <- x_name
                  res
              }, probs, names(probs), SIMPLIFY = FALSE)
              
              probs <- Reduce(function(x, y) merge(x, y, by = "barcode"), probs)
              rownames(probs) <- probs$barcode
              probs$barcode <- NULL
              
              probs <- probs[match(barcodes, rownames(probs)), ]
              probs
          })


#' @title Get metadata from scPred object
#' @description Accessor function to retrieve metadata from scPred object
#' @return A dataframe including the cell barcodes and prediction variable 
#' (cell type labels)
#' @export

setMethod("get_metadata", 
          signature("scPred"), 
          function(object){
              object@metadata
          })


#' @title Get metadata from scPred object
#' @description Accessor function to retrieve metadata from scPred object
#' @return A dataframe including the cell barcodes and prediction variable 
#' (cell type labels)
#' @export

setMethod("get_metadata", 
          signature("Seurat"), 
          function(object){
              get_metadata(get_scpred(object))
          })