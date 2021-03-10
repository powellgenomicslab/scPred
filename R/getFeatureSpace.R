#' @title Get discriminant feature space
#' @description Given a prediction variable, finds a feature set of class-informative principal components that
#' explain variance differences between cell types. 
#' @param object A \code{seurat} object
#' @param pvar Column in \code{meta.data} slot containing the cell-type labels of each single cell
#' @param correction Multiple testing correction method used. Default: false discovery rate. See \code{p.adjust} function 
#' @param sig Significance level to determine principal components explaining class identity
#' @param reduction Name of reduction in Seurat objet to be used to determine the feature space. Default: "pca"
#' @return An \code{Seurat} object along with a \code{scPred} object stored in the \code{@misc} slot 
#' containing a data.frame of significant features with the following columns:
#' \itemize{
#' \item PC: Principal component
#' \item pValue: p-value obtained from Wilcoxon rank sum test
#' \item pValueAdj: Adjusted p-value according to \code{correction} parameter
#' \item expVar: Explained variance for each principal component
#' \item cumExpVar: All principal components are ranked according to their significance and variance explained. 
#' This column contains the cumulative variance explained across the ranked principal components
#' }
#' @keywords informative, significant, features
#' @importFrom methods is
#' @importFrom tidyr gather
#' @importFrom magrittr "%>%"
#' @importFrom dplyr mutate arrange filter distinct
#' @importFrom pbapply pblapply
#' 
#' @export
#' @author
#' Jose Alquicira Hernandez
#' 
#' @examples 
#' 
#' library(scPred)
#' pbmc_small <- getFeatureSpace(pbmc_small, "RNA_snn_res.0.8")
#' 




getFeatureSpace <- function(object, pvar, correction = "fdr", sig = 1, reduction = "pca"){
  
  
  # Validations -------------------------------------------------------------
  
  if(!is(object, "Seurat")){
    stop("Invalid class for object: must be 'Seurat'")
  }
  
  if(!any(correction %in% stats::p.adjust.methods)){
    stop("Invalid multiple testing correction method. See ?p.adjust function")
  }
  
  if(!pvar %in% names(object@meta.data)){
    stop("Prediction variable is not stored in metadata slot")
  }
  
  classes <- object[[pvar, drop = TRUE]]
  
  if(!is.factor(classes)){
    classes <- as.factor(classes)
  }
  
  assay <- DefaultAssay(object)
  # Filter principal components by variance ---------------------------------
  
  # Check if a PCA has been computed
  if(!(reduction %in% Reductions(object))){
    stop("No ",reduction, " reduction has been computet yet. See RunPCA() function?")
  }else{
    reduction_data <- Reductions(object, slot = reduction)
    if(reduction_data@assay.used != assay) 
      stop("No ", 
           reduction, 
           " reduction associated with default assay: ", 
           assay, "\nChange default assay or compute a new reduction")
  }
  
  # Check if available was normalized
  
  cellEmbeddings <- Embeddings(reduction_data)
  loadings <- Loadings(reduction_data)
  reduction_key <- reduction_data@key
  
  # Store original labels in metadata slot
  
  spmodel <- new("scPred", metadata = data.frame(pvar = classes))
  
  # Validate response variable values
  original_classes <- classes
  uniqueClasses <- unique(classes)
  isValidName <- uniqueClasses == make.names(uniqueClasses)
  
  if(!all(isValidName)){
    classes <- .make_names(classes)
    classes <- factor(classes, levels = unique(classes))
    names(classes) <- Cells(object)
  }
  
  spmodel@metadata$response <- classes

  # Get means and sds -------------------------------------------------------

  features <- rownames(loadings)
  
  data <- GetAssayData(object, "data", assay = assay)[features,]
  means <- Matrix::rowMeans(data)
  
  rowVar <- function(x, ...) {
    sqrt(Matrix::rowSums((x - means)^2, ...)/(ncol(x) - 1))
  }
  
  stdevs  <- rowVar(data)
  
  i <- stdevs == 0
  
  if(any(i)){
    warning(paste0(sum(i), " genes have zero variance but are present in the gene loadings. \nDid you subset or integrated this data before?"))
    cat(crayon::yellow("Removing zero-variance genes from loadings\n"))
    
    loadings <- loadings[!i,]
    means <- means[!i]
    stdevs <- stdevs[!i]
  }
  
  
  spmodel@scaling <- data.frame(means, stdevs)
  
  
  # Select informative principal components
  # If only 2 classes are present in prediction variable, train one model for the positive class
  # The positive class will be the first level of the factor variable
  
  cat(crayon::green(cli::symbol$record, " Extracting feature space for each cell type...\n"))
  if(length(levels(classes)) == 2){
    
    message("First factor level in '", pvar, "' metadata column considered as positive class:")
    message(levels(original_classes)[1])
    res <- .getFeatures(make.names(levels(original_classes)[1]), classes, cellEmbeddings, correction, sig)
    res <- list(res)
    names(res) <- levels(original_classes)[1]
    
  }else{
    
    res <- pblapply(levels(classes), .getFeatures, classes, cellEmbeddings, correction, sig)
    dict <- data.frame(classes, original_classes) %>% 
      distinct()
    
    i <- match(levels(classes), dict$classes)
    names(res) <- as.character(dict$original_classes[i])
    
  }
  
  
  nFeatures <- unlist(lapply(res, nrow))
  
  noFeatures <- nFeatures == 0
  
  if(any(noFeatures)){
    
    warning("\nWarning: No features were found for classes:\n",
            paste0(names(res)[noFeatures], collapse = "\n"), "\n")
    res[[names(res)[noFeatures]]] <- NULL
    
  }
  
  
  
  # Create scPred object
  
  spmodel@pvar <- pvar
  spmodel@features <- res
  spmodel@cell_embeddings <- cellEmbeddings
  spmodel@feature_loadings <- loadings
  spmodel@reduction <- reduction
  spmodel@reduction_key <- reduction_key
  
  object@misc$scPred <- spmodel
  
  cat(crayon::green("DONE!\n"))
  
  object
  
  
}

.getFeatures <- function(positiveClass, classes, cellEmbeddings, correction, sig){
  
  # Set non-positive classes to "other"
  i <- classes != positiveClass
  newClasses <- as.character(classes)
  newClasses[i] <- "other"
  newClasses <- factor(newClasses, levels = c(positiveClass, "other"))
  
  # Get indices for positive and negative class cells
  positiveCells <- newClasses == positiveClass
  negativeCells <- newClasses == "other"
  
  # Get informative features
  apply(cellEmbeddings, 2, function(d) stats::wilcox.test(d[positiveCells], d[negativeCells])) %>%
    lapply('[[', "p.value") %>% # Extract p-values
    as.data.frame() %>% 
    gather(key = "feature", value = "pValue") %>%
    mutate(pValueAdj = stats::p.adjust(pValue, method = correction, n = nrow(.))) %>% # Perform multiple test correction
    arrange(pValueAdj) %>% 
    filter(pValueAdj < sig) -> sigDims
  
  sigDims
}

