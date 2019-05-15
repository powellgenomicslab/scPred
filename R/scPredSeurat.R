#' @title Perform predictions via Seurat integration
#' @description Integrates data using Seurat and performs predictions based on informative features
#' identified by \code{scPred}
#' @param reference A \code{Seurat} object used as reference to perform predictions
#' @param new A \code{Seurat} object from which cells will be classified
#' @param threshold Threshold used for probabilities to classify cells into classes
#' @param weight Whether output probabilities are averaged across models for every cell.
#' @param ... Extra arguments passed to \code{FindTransferAnchors()}
#' @return A matrix:.
#' \itemize{
#' \item Probabilities of each cell type for all cells in the \code{new} object
#' \item \code{prediction} Cell type predicted based on probability threshold
#' \item \code{generic} Cell type predicted based on maximum probability acroos all models
#' }
#' @importFrom magrittr "%>%"
#' @export
#' @author
#' José Alquicira Hernández
#' 
#' @examples 
#' 
#' res <- scPredSeurat(reference = integrated, new = query)
#' 

scPredSeurat <- function(reference, new, threshold = NULL, weight = TRUE, ...){
  
  # Validate input parameters
  
  if(!is(reference, "Seurat")){
    stop("Invalid class for 'reference': must be 'Seurat'")
  }
  if(!is(new, "Seurat")){
    stop("Invalid class for 'new': must be 'Seurat'")
  }
  

  # "Align" data using Seurat
  anchors <- Seurat::FindTransferAnchors(reference = reference, 
                                 query = new, 
                                 dims = seq_len(ncol(reference@reductions$pca@cell.embeddings)), 
                                 l2.norm = FALSE, 
                                 ...)
  # Extract cell embeddings
  allEmbeddings <- anchors@object.list[[1]]@reductions$pcaproject@cell.embeddings
  
  # Subset test loading projection
  i <- grepl("_query$", rownames(allEmbeddings))
  testEmbeddings <- allEmbeddings[i, ]  
  
  # Classify cells using all trained models 
  cellTypeModelNames <- names(reference@misc$scPred$features)
  res <- sapply(cellTypeModelNames, .predictCellClass, model, reference, testEmbeddings)
  
  # Gather results
  res <- Reduce(cbind, res)
  colnames(res) <- cellTypeModelNames
  rownames(res) <- colnames(new)
  
  # Make sum of probabilities across models equal to zero
  if(weight){
    res <- res/rowSums(res)
  }
  
  
  # Extract maximum probability for each class
  i <- apply(res, 1, which.max)
  resSplit <- split(res, row(res))
  maxProps <- mapply(function(x,y) x[y], resSplit, i)
  
  # Store classification based on maximum probability
  genericPred <- cellTypeModelNames[i]
  
  # Classify cells according to probability threshold
  
  if(is.null(threshold)){
    classThreshold <- sapply(cellTypeModelNames, .getThreshold, reference)
    pred <- data.frame(maxProps, genericPred)
    pred <- ifelse(pred$maxProps > classThreshold[pred$genericPred], 
                   as.character(pred$genericPred), 
                   "unassigned")
    
  }else{
    
  pred <- ifelse(maxProps > threshold, genericPred, "unassigned")
  }
  
  names(pred) <- colnames(new)
  res <- as.data.frame(res)
  
  # Format results
  res$generic <- genericPred
  res$prediction <- pred
  
  # Return results
  res
   
}




.predictCellClass <-  function(cellType, model, reference, testEmbeddings){
  
  # Extract features for a given cell type
  as.character(reference@misc$scPred$features[[cellType]]$PC) -> features
  
  # Format test cell embeddings
  testEmbeddings %>% 
    colnames() %>% 
    gsub("Project", "", .) %>% 
    `colnames<-`(testEmbeddings, .) -> testEmbeddings
  
  # Extract cell type model
  model <- reference@misc$scPred$train[[cellType]]
  
  # Perform predictions based on informative PCs
  prediction <- predict(model, 
                        newdata = subsetMatrix(testEmbeddings, features), 
                        type = "prob")
  
  # Add cell names to results
  rownames(prediction) <- rownames(testEmbeddings)
  
  # Return positive-class probability
  prediction[,1, drop = FALSE]
  
}


#' @title Get thresholds from training models
#' @description Calculates the threshold to classify a cell type.
#' @param cellType Cell type of interest
#' @param reference A \code{Seurat} object used as reference to perform predictions
#' @importFrom magrittr "%>%"
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr pull

.getThreshold <- function(cellType, reference){
  props <- as.data.frame(reference@misc$scPred$train[[cellType]]$pred)
  bestTune <- reference@misc$scPred$train[[cellType]]$bestTune
  
  mapply(bestTune, names(bestTune), 
         FUN = function(par, name) paste0(name, " == ", par)) %>% 
    paste0(collapse = " & ") -> subsetRule 
  
  
  props %>% 
    filter(!!rlang::parse_expr(subsetRule)) %>% 
    select(other, obs) %>% 
    group_by(obs) %>% 
    summarize(median = median(other)) %>% 
    pull(median) %>% 
    mean() 
}
