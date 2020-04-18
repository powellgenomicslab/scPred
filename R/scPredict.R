#' @title Predict cell classes from a new dataset using a trained model
#' @description Predicts cell classes for a new dataset based on trained model(s)
#' @param reference An \code{Seurat} object with trained model(s) using \code{scPred}.
#' @param new A seurat object containing cells to be classified
#' @param threshold Threshold used for probabilities to classify cells into classes
#' @param weight If TRUE, probabilities from binary classifiers are weighted to sum 1.
#' @param max.iter.harmony Maximum number of rounds to run Harmony. One round of Harmony involves one clustering and one correction step.
#' @return A Seurat object with addtional metadata columns with prediction probabilities associated to each class, a \code{prediction} column, 
#' indicating the classification based on the provided threshold and a \code{generic_class} column without "unassigned" labels.
#' @keywords prediction, new, test, validation
#' @importFrom methods is
#' @importFrom stats predict
#' @importFrom Seurat Embeddings Stdev ScaleData AddMetaData DefaultAssay CreateDimReducObject VariableFeatures Loadings GetAssayData
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom dplyr mutate select distinct mutate_all
#' @importFrom pbapply pblapply
#' @importFrom harmony HarmonyMatrix
#' @export
#' @author
#' José Alquicira Hernández
#' @examples 
#' 
#' \dontrun{
#' new <- scPredict(reference, new)
#' }



scPredict <- function(reference, 
                      new, 
                      threshold = NULL, 
                      weight = TRUE, 
                      max.iter.harmony = 20,
                      ...){
  
  # Function validations ----------------------------------------------------
  
  # Validate if provided object is an scPred object
  if(!is(reference, "Seurat")) stop("'object' must be of class 'scPred' or 'Seurat'")
  if(!"scPred" %in% names(reference@misc)) stop("No feature space has been determined!")
  if(!length(reference@misc$scPred@train)) stop("No models have been trained!")
  if(!is(new, "Seurat")) stop("Both reference and new data must be Seurat objects")
  
  
  cat(crayon::green(cli::symbol$record, " Matching reference with new dataset...\n"))
  
  # Subset data
  ref_loadings <- Loadings(reference, reduction = "pca")
  ref_embeddings <- Embeddings(reference, reduction = "pca")
  new_genes <- rownames(new)
  
  # Get genes
  reference_genes <- rownames(ref_loadings)
  
  
  # Get intersection between reference and new datasets
  shared_genes <- intersect(reference_genes, new_genes)
  cat(crayon::cyan("\t", cli::symbol$line, paste(length(reference_genes), "genes present in reference loadings\n")))
  cat(crayon::cyan("\t", cli::symbol$line, paste(length(shared_genes), "genes shared between reference and new dataset\n")))
  cat(crayon::cyan("\t", cli::symbol$line, paste0(round(length(shared_genes)/length(reference_genes) * 100, 2), 
                 "% of genes in the reference are present in new dataset\n")))
  
  
  # Subset shared genes from reference
  ref_loadings <- ref_loadings[shared_genes, ]
  
  
  # Scale new data
  new <- ScaleData(new, features = reference_genes, ...)
  #new <- ScaleData(new, features = reference_genes)
  
  
  ## Subset shared genes from new dataset
  
  new_data <- GetAssayData(new, "scale.data")
  new_data <- new_data[shared_genes, ]
  new_data <- new_data[match(shared_genes, rownames(new_data)), ]
  
  # all(rownames(new_data) == shared_genes)
  
  new_embeddings <- t(new_data) %*% ref_loadings
  
  
  dataset <- factor(c(rep("reference", nrow(ref_embeddings)), rep("new", nrow(new_embeddings))), 
                    levels = c("reference", "new"))
  
  rownames(ref_embeddings) <- paste0("ref_", rownames(ref_embeddings))
  rownames(new_embeddings) <- paste0("new_", rownames(new_embeddings))
  
  
  eigenspace <- as.data.frame(rbind(ref_embeddings, new_embeddings))
  meta_data <- data.frame(rownames(eigenspace), dataset = dataset)
  
  cat(crayon::green(cli::symbol$record, " Aligning new data to reference...\n"))
  
  harmony_embeddings <- HarmonyMatrix(eigenspace, 
                                      meta_data, 
                                      'dataset', 
                                      do_pca = FALSE, 
                                      reference_values = "reference",
                                      max.iter.harmony = max.iter.harmony)
  
  new_embeddings_aligned <- harmony_embeddings[dataset == "new", ]
  
  # Classify cells using all trained models 
  cellTypeModelNames <- names(reference@misc$scPred@features)
  .predictCellClass <-  function(cellType, reference, testEmbeddings){
    
    # Extract features for a given cell type
    as.character(reference@misc$scPred@features[[cellType]]$PC) -> features
    
    # Format test cell embeddings
    namesEmbeddings <- gsub("Project", "", colnames(testEmbeddings))
    colnames(testEmbeddings) <- namesEmbeddings
    
    # Extract cell type model
    model <- reference@misc$scPred@train[[cellType]]
    
    # Perform predictions based on informative PCs
    prediction <- predict(model, 
                          newdata = subsetMatrix(testEmbeddings, features), 
                          type = "prob")
    
    # Add cell names to results
    rownames(prediction) <- rownames(testEmbeddings)
    
    # Return positive-class probability
    prediction[,1, drop = FALSE]
    
  }
  
  cat(crayon::green(cli::symbol$record, " Classifying cells...\n"))
  res <- sapply(cellTypeModelNames, .predictCellClass, reference, new_embeddings_aligned)
  
  # Gather results
  res <- Reduce(cbind, res) %>% as.data.frame()
  colnames(res) <- cellTypeModelNames
  rownames(res) <- colnames(new)
  
  classes <- cellTypeModelNames
  #plot(res$Lymphoid, col = as.factor(test$CellType))
  # If there is only 2 classes, compute complementary probability for negative class
  if(length(cellTypeModelNames) == 1){
    
    cellClasses <- levels(reference$scPred_response)
    res_comp <- 1 - res[,1]
    negClass <- cellClasses[cellClasses != names(res)]
    res[[negClass]] <- res_comp
    
  }else if(weight){
    # Make sum of probabilities across models equal to zero
    res <- res/rowSums(res)
  }
  
  # Extract maximum probability for each class
  max_props <- as.data.frame(t(apply(res, 1, function(x) c(index = which.max(x),  max = x[which.max(x)]))))
  names(max_props) <- c("index", "max")
  
  
  # Store classification based on maximum probability
  max_props$generic_class <- names(res)[max_props$index]
  res <- cbind(res, max_props)
  
  
  # Relabel cells according to original values
  pVar <- reference@misc$scPred@pVar
  dict <- reference@meta.data[c(pVar, "scPred_response")] %>% 
    distinct() %>% 
    mutate_all(as.character)
  
  
  res$no_rejection <- dict[match(res$generic_class, dict[,2]), 1]
  

  
  # Classify cells according to probability threshold
  
  if(is.null(threshold)){
    classThreshold <- sapply(cellTypeModelNames, .getThreshold, reference)
    
    if(length(classThreshold) == 1){
      cellClasses <- levels(reference$scPred_response)
      classThreshold <- c(classThreshold, classThreshold)
      names(classThreshold) <- cellClasses
    }
    
    dict <- cbind(dict, classThreshold)
    dict[[res$generic_class %in% dict["scPred_response"] , "threshold"]]
    
    pred <- ifelse(res$max >  dict[dict$scPred_response %in% res$generic_class, "classThreshold"], 
                   as.character(res$no_rejection), 
                   "unassigned")
    
  }else{
    
    pred <- ifelse(res$max > threshold, res$no_rejection , "unassigned")
  }
  
  names(pred) <- colnames(new)
  
  # Format results
  res$prediction <- pred
  res$index <- NULL
  res$generic_class <- NULL
  
  names(res) <- paste0("scPred_", names(res))
  
  cat(crayon::green("DONE!\n"))
  
  # Return results
  res <- AddMetaData(new, res)
  rownames(new_embeddings_aligned) <- gsub("^new_", "", rownames(new_embeddings_aligned))
  res@reductions[["scPred"]] <- CreateDimReducObject(embeddings = new_embeddings_aligned, 
                                                     key = "harmony_",
                                                     assay = DefaultAssay(object = new))
  res
  
  
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
  props <- as.data.frame(reference@misc$scPred@train[[cellType]]$pred)
  bestTune <- reference@misc$scPred@train[[cellType]]$bestTune
  
  mapply(bestTune, names(bestTune), 
         FUN = function(par, name) paste0(name, " == ", par)) %>% 
    paste0(collapse = " & ") -> subsetRule 
  
  
  props %>% 
    #filter(!!rlang::parse_expr(subsetRule)) %>% 
    select(other, obs) %>% 
    group_by(obs) %>% 
    summarize(median = median(other)) %>% 
    pull(median) %>% 
    mean() 
}




