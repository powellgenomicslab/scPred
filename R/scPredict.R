#' @title Predict cell classes from a new dataset using a trained model
#' @description Predicts cell classes for a new dataset based on trained model(s)
#' @param new A seurat object containing cells to be classified
#' @param reference An \code{Seurat} object with trained model(s) using \code{scPred} or an \code{scPred} object
#' @param threshold Threshold used for probabilities to classify cells into classes. All cells below 
#' this threshold value will be labels as "unassigned". In the case of binary classification (two cell tyoes),
#' a threshold of \code{0.5} will force all cells to be classified to any of the two cell types. For multi-class
#' classification, if there's no probability higher than the threshold associated to a cell type, this will
#' be labelled as "unassigned"
#' @param max.iter.harmony Maximum number of rounds to run Harmony. One round of Harmony involves one clustering and one correction step.
#' @param recompute_alignment Recompute alignment? Useful if \code{scPredict()} has already been run
#' @param seed Numeric seed for harmony 
#' @return A Seurat object with addtional metadata columns with prediction probabilities associated to each class, a \code{prediction} column, 
#' indicating the classification based on the provided threshold and a \code{generic_class} column without "unassigned" labels.
#' @keywords prediction, new, test, validation
#' @importFrom methods is
#' @importFrom stats predict
#' @importFrom SeuratObject Embeddings Stdev AddMetaData DefaultAssay CreateDimReducObject VariableFeatures Loadings GetAssayData
#' @importFrom pbapply pblapply
#' @importFrom harmony HarmonyMatrix
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
  
  
  
  # Dataset alignment -------------------------------------------------------
  
  if("scpred" %in% names(new@reductions)){
    if(recompute_alignment){
      alignment <- TRUE
      crayon::yellow(cli::symbol$figure_dash, "Data has already been aligned to a reference.\n") %>% 
        cat(sep = "")
      crayon::yellow(cli::symbol$sup_plus, "Skip data alignment using `recompute.alignment = FALSE`.\n") %>% 
        cat(sep = "")
    } 
    else {
      alignment <- FALSE
    }
    
  }else{
    alignment <- TRUE
  }
  
  if(alignment){
    
    cat(crayon::green(cli::symbol$record, " Matching reference with new dataset...\n"))
    
    # Subset data
    ref_loadings <- spmodel@feature_loadings
    ref_embeddings <- spmodel@cell_embeddings
    new_features <- rownames(new)
    
    # Get genes
    reference_features <- rownames(ref_loadings)
    
    
    # Get intersection between reference and new datasets
    shared_features <- intersect(reference_features, new_features)
    cat(crayon::cyan("\t", cli::symbol$line, paste(length(reference_features), "features present in reference loadings\n")))
    cat(crayon::cyan("\t", cli::symbol$line, paste(length(shared_features), "features shared between reference and new dataset\n")))
    cat(crayon::cyan("\t", cli::symbol$line, paste0(round(length(shared_features)/length(reference_features) * 100, 2), 
                                                    "% of features in the reference are present in new dataset\n")))
    
    
    # Subset shared genes from reference
    ref_loadings <- ref_loadings[shared_features, ]
    
    
    
      new_data <- GetAssayData(new, "data")[shared_features,]
      means <- spmodel@scaling$means
      stdevs  <- spmodel@scaling$stdevs
      new_data <- Matrix::t(new_data)
      names(means) <- rownames(spmodel@scaling) -> names(stdevs)
      
      # Subset means and standard deviations
      means <- means[shared_features]
      stdevs <- stdevs[shared_features]
      #all(colnames(new_data) == names(means))

      i <- stdevs == 0
      
      if(any(i)){
        warning(paste0(sum(i), " features have zero variance but are present in the feature loadings. \nDid you subset or integrated this data before?"))
        cat(crayon::yellow("Removing zero-variance genes from projection\n"))
        
        new_data <- new_data[,!i]
        ref_loadings <- ref_loadings[!i,]
        means <- means[!i]
        stdevs <- stdevs[!i]
        #all(colnames(new_data) == rownames(ref_loadings))
      }
      
      
      scaled_data <- scale(new_data, means, stdevs)
      
      
      new_embeddings <- scaled_data %*% ref_loadings
      
    
    
    dataset <- factor(c(rep("reference", nrow(ref_embeddings)), rep("new", nrow(new_embeddings))), 
                      levels = c("reference", "new"))
    
    
    rownames(ref_embeddings) <- paste0("ref_", rownames(ref_embeddings))
    rownames(new_embeddings) <- paste0("new_", rownames(new_embeddings))
    
    
    eigenspace <- as.data.frame(rbind(ref_embeddings, new_embeddings))
    meta_data <- data.frame(rownames(eigenspace), dataset = dataset)
    
    cat(crayon::green(cli::symbol$record, " Aligning new data to reference...\n"))
    
    set.seed(seed)
    harmony_embeddings <- HarmonyMatrix(eigenspace, 
                                        meta_data, 
                                        'dataset', 
                                        do_pca = FALSE, 
                                        reference_values = "reference",
                                        max.iter.harmony = max.iter.harmony)
    
    new_embeddings_aligned <- harmony_embeddings[dataset == "new", , drop = FALSE]
    
  }else{
    new_embeddings_aligned <- Embeddings(new, reduction = "scpred")
    colnames(new_embeddings_aligned) <- gsub("scpred_", spmodel@reduction_key, colnames(new_embeddings_aligned))
  }
  
  
  
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
  res <- as.data.frame(Reduce(cbind, res))
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
  rownames(new_embeddings_aligned) <- gsub("^new_", "", rownames(new_embeddings_aligned))
  new@reductions[["scpred"]] <- CreateDimReducObject(embeddings = new_embeddings_aligned, 
                                                     key = "scpred_",
                                                     assay = DefaultAssay(object = new))
  if(recompute_alignment){
    
    rownames(new_embeddings) <- gsub("^new_", "", rownames(new_embeddings))
    new@reductions[["scpred_projection"]] <- CreateDimReducObject(embeddings = new_embeddings, 
                                                           key = "Projection_",
                                                           assay = DefaultAssay(object = new))
  }
  
  cat(crayon::green("DONE!\n"))
  
  new
  
  
}
