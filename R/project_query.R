#' @title Project query data onto PCA reference
#' @description Projects data from Seurat object onto a pre-computed PCA
#' @param new A Seurat object containing cells to be projected
#' @param reference A \code{Seurat} object after running \code{getFeatureSpace()}
#' @param max.iter.harmony Maximum number of rounds to run Harmony. One round of Harmony involves one clustering and one correction step.
#' @param recompute_alignment Recompute alignment? Useful if \code{scPredict()} has already been run
#' @param seed Numeric seed for harmony 
#' @param ... Further arguments to \code{HarmonyMatrix()}
#' @return A Seurat object with two new reductions: 
#' \itemize{
#' \item \code{scpred}: Aligned data using harmony
#' \item \code{scpred_projection}: Raw projection using reference loadings
#' }
#' @keywords prediction, new, test, validation
#' @importFrom methods is
#' @importFrom SeuratObject Embeddings Stdev DefaultAssay CreateDimReducObject VariableFeatures Loadings GetAssayData
#' @importFrom harmony HarmonyMatrix
#' @export
#' @author
#' José Alquicira Hernández


project_query <- function(new, 
                          reference, 
                          max.iter.harmony = 20,
                          recompute_alignment = TRUE,
                          seed = 66, 
                          ...){


# Validate if provided object is an scPred object
if(!(is(reference, "Seurat") | is(reference, "scPred"))) stop("'object' must be of class 'scPred' or 'Seurat'")

if(is(reference, "Seurat")){
  spmodel <- reference@misc$scPred
}else{
  spmodel <- reference
}

if(is.null(spmodel)) stop("No feature space has been determined!")
if(!is(new, "Seurat")) stop("New data must be a Seurat object")


# Dataset alignment -------------------------------------------------------

if("scpred" %in% names(new@reductions)){
  if(recompute_alignment){
    alignment <- TRUE
    cat(crayon::yellow(cli::symbol$figure_dash, "Data has already been aligned to a reference.\n"), sep = "")
    cat(crayon::yellow(cli::symbol$sup_plus, "Skip data alignment using `recompute.alignment = FALSE`.\n"),  sep = "")
      
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
  
  
  
  new_data <- GetAssayData(new, layer="data")[shared_features,]
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
                                      max.iter.harmony = max.iter.harmony, 
                                      ...)
  
  new_embeddings_aligned <- harmony_embeddings[dataset == "new", , drop = FALSE]
  
}else{
  new_embeddings_aligned <- Embeddings(new, reduction = "scpred")
  colnames(new_embeddings_aligned) <- gsub("scpred_", spmodel@reduction_key, colnames(new_embeddings_aligned))
}

  
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

  new
  
}
