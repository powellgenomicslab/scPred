#' @title Get discriminant feature space
#' @description Given a prediction variable, finds a feature set of class-informative principal components that
#' explain variance differences between cell types. 
#' @param object A \code{seurat} object
#' @param pVar Column in \code{meta.data} slot containing the cell-type labels of each single cell
#' @param correction Multiple testing correction method used. Default: false discovery rate. See \code{p.adjust} function 
#' @param sig Significance level to determine principal components explaining class identity
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
#' @importFrom dplyr mutate arrange filter
#' @importFrom pbapply pblapply
#' @importFrom robustbase covMcd 
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




getFeatureSpace <- function(object, pVar, correction = "fdr", sig = 1, reduction = "pca"){
  
  
  # Validations -------------------------------------------------------------
  
  if(!is(object, "Seurat")){
    stop("Invalid class for object: must be 'Seurat'")
  }
  
  if(!any(correction %in% stats::p.adjust.methods)){
    stop("Invalid multiple testing correction method. See ?p.adjust function")
  }
  
  if(!pVar %in% names(object@meta.data)){
    stop("Prediction variable is not stored in metadata slot")
  }
  
  classes <- object[[pVar, drop = TRUE]]
  
  if(!is.factor(classes)){
    classes <- as.factor(classes)
  }
  
  # Filter principal components by variance ---------------------------------
  
  # Check if a PCA has been computed
  if(!(reduction %in% names(object@reductions))){
    stop("No ",reduction, " reduction has been computet yet. See RunPCA() function?")
  }
  
  # Check if available was normalized
  
  assay <- DefaultAssay(object)
  cellEmbeddings <- Embeddings(object, reduction = reduction)
  
  
  # Validate response variable values
  uniqueClasses <- unique(classes)
  isValidName <- uniqueClasses == make.names(uniqueClasses)
  
  if(!all(isValidName)){
    classes <- make.names(classes)
    classes <- factor(classes, levels = unique(classes))
    names(classes) <- Cells(object)
  }
  
  object@meta.data[["scPred_response"]] <- classes
  
  # Select informative principal components
  # If only 2 classes are present in prediction variable, train one model for the positive class
  # The positive class will be the first level of the factor variable
  
  cat(crayon::green(cli::symbol$record, " Extracting feature space for each cell type...\n"))
  if(length(levels(classes)) == 2){
    
    message("First factor level in '", pVar, "' metadata column considered as positive class:")
    message(levels(classes)[1])
    res <- .getFeatures(levels(classes)[1], classes, cellEmbeddings, correction, sig)
    res <- list(res)
    names(res) <- levels(classes)[1]
    
  }else{
    
    res <- pblapply(levels(classes), .getFeatures, classes, cellEmbeddings, correction, sig)
    names(res) <- levels(classes)
    
  }
  
  
  nFeatures <- unlist(lapply(res, nrow))
  
  noFeatures <- nFeatures == 0
  
  if(any(noFeatures)){
    
    warning("\nWarning: No features were found for classes:\n",
            paste0(names(res)[noFeatures], collapse = "\n"), "\n")
    res[[names(res)[noFeatures]]] <- NULL
    
  }
  
  cat(crayon::green("DONE!\n"))
  
  
  
  # Create scPred object
  
  object@misc$scPred <- new("scPred", 
                            pVar = pVar,
                            reduction = reduction,
                            features = res)
  
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

