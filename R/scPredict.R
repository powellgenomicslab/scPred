#' @title Predict cell classes from a new dataset using a trained model
#' @description Predicts cell classes for a new dataset based on trained model(s)
#' @param new A seurat object containing cells to be classified
#' @param reference A \code{Seurat} object with trained model(s) using \code{scPred} or an \code{scPred} object
#' @param threshold Threshold used for probabilities to classify cells into classes. All cells below 
#' this threshold value will be labels as "unassigned". In the case of binary classification (two cell tyoes),
#' a threshold of \code{0.5} will force all cells to be classified to any of the two cell types. For multi-class
#' classification, if there's no probability higher than the threshold associated to a cell type, this will
#' be labelled as "unassigned"
#' @param max.iter.harmony Maximum number of rounds to run Harmony. One round of Harmony involves one clustering and one correction step.
#' @param recompute_alignment Recompute alignment? Useful if \code{scPredict()} has already been run
#' @param seed Numeric seed for harmony 
#' @return A Seurat object with additional metadata columns with prediction probabilities associated to each class, a \code{prediction} column, 
#' indicating the classification based on the provided threshold and a \code{generic_class} column without "unassigned" labels. Additionally,
#' two new reductions are returned: 
#' \itemize{
#' \item \code{scpred}: Aligned data using harmony
#' \item \code{scpred_projection}: Raw projection using reference loadings
#' }
#' @keywords prediction, new, test, validation
#' @importFrom methods is
#' @importFrom stats predict
#' @importFrom SeuratObject Embeddings AddMetaData
#' @importFrom pbapply pblapply
#' @export
#' @author
#' José Alquicira Hernández


scPredict <- function(new,
                      reference, 
                      threshold = 0.55, 
                      max.iter.harmony = 20,
                      recompute_alignment = TRUE,
                      seed = 66){
  
  # Function validations ----------------------------------------------------
  
  # Validate if provided object is an scPred object
  if(!(is(reference, "Seurat") | is(reference, "scPred"))) stop("'object' must be of class 'scPred' or 'Seurat'")
  
  if(is(reference, "Seurat")){
    spmodel <- reference@misc$scPred
  }else{
    spmodel <- reference
  }
  
  if(is.null(spmodel)) stop("No feature space has been determined!")
  if(!length(spmodel@train)) stop("No models have been trained!")
  if(!is(new, "Seurat")) stop("New data must be a Seurat object")
  

  # Project query data ------------------------------------------------------

  new <- project_query(new, 
                reference = spmodel, 
                max.iter.harmony = max.iter.harmony, 
                recompute_alignment = recompute_alignment,
                seed = seed)
  
  new_embeddings_aligned <- Embeddings(new[["scpred"]])
  colnames(new_embeddings_aligned) <- colnames(spmodel@cell_embeddings)
  
  
  
  # Classify cells using all trained models 
  cellTypeModelNames <- names(spmodel@features)
  .predictCellClass <-  function(cellType, spmodel, testEmbeddings){
    
    # Extract features for a given cell type
    as.character(spmodel@features[[cellType]]$feature) -> features
    
    # Extract cell type model
    model <- spmodel@train[[cellType]]
    
    # Perform predictions based on informative PCs
    prediction <- predict(model, 
                          newdata = scPred:::subsetMatrix(testEmbeddings, features), 
                          type = "prob")
    
    # Add cell names to results
    rownames(prediction) <- rownames(testEmbeddings)
    
    # Return positive-class probability
    prediction[,1, drop = FALSE]
    
  }
  
  
  cat(crayon::green(cli::symbol$record, " Classifying cells...\n"))
  res <- sapply(cellTypeModelNames, .predictCellClass, spmodel, new_embeddings_aligned)
  
  # Gather results
  res <- as.data.frame(res)
  colnames(res) <- cellTypeModelNames
  rownames(res) <- colnames(new)
  
  classes <- cellTypeModelNames
  #plot(res$Lymphoid, col = as.factor(test$CellType))
  # If there is only 2 classes, compute complementary probability for negative class
  if(length(cellTypeModelNames) == 1){
    metadata <- get_metadata(spmodel)
    cellClasses <- levels(metadata$pvar)
    res_comp <- 1 - res[,1]
    negClass <- cellClasses[cellClasses != names(res)]
    res[[negClass]] <- res_comp
    
  }
  
  # Extract maximum probability for each class
  max_props <- as.data.frame(t(apply(res, 1, function(x) c(index = which.max(x),  max = x[which.max(x)]))))
  names(max_props) <- c("index", "max")
  
  
  # Store classification based on maximum probability
  max_props$generic_class <- names(res)[max_props$index]
  res <- cbind(res, max_props)
  
  # Classify cells according to probability threshold
    
  pred <- ifelse(res$max > threshold, res$generic_class, "unassigned")
  
  names(pred) <- colnames(new)
  
  # Format results
  res$prediction <- pred
  res$index <- NULL
  res$no_rejection <- res$generic_class
  res$generic_class <- NULL
  
  names(res) <- .make_names(paste0("scpred_", names(res)))
  
  
  # Return results
  new <- AddMetaData(new, res)
  
  cat(crayon::green("DONE!\n"))
  
  new
  
  
}
