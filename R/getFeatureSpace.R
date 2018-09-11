#' @title Get informative principal components
#' @description Given a prediction variable, finds a feature set of class-informative principal components. 
#' A Wilcoxon rank sum test is used to determine a difference between the score distributions of cell classes from the prediction variable.
#' @param object An \code{scPred} object
#' @param pVar Prediction variable corresponding to a column in \code{metadata} slot
#' @param varLim Threshold to filter principal components based on variance explained.
#' @param correction Multiple testing correction method. Default: false discovery rate. See \code{p.adjust} function 
#' @param sig Significance level to determine principal components explaining class identity
#' @return A \code{scPred} object with two additional filled slots:
#' \itemize{
#' \item \code{features}: A data frame with significant principal components the following information:
#' \itemize{
#' \item PC: Principal component
#' \item pValue: p-value obtained from Mann-Whitney test
#' \item pValueAdj: Adjusted p-value according to \code{correction} parameter
#' \item expVar: Explained variance by the principal component
#' \item cumExpVar: All principal components are ranked accoriding to their frequency of ocurrence and their variance explained. 
#' This column contains the cumulative variance explained across the ranked principal components
#' }
#' \item \code{pVar}: Column name from metadata to use as the variable to predict using
#' the informative principal components. Informative principal components are selected based on this variable.
#' }
#' @keywords informative, significant, features
#' @importFrom methods is
#' @export
#' @author
#' José Alquicira Hernández
#' 
#' @examples 
#' 
#' # Assign cell information to scPred object
#' # Cell information must be a data.frame with rownames as cell ids matching the eigendecomposed 
#' gene expression matrix rownames.
#' 
#' metadata(object) <- cellInfo
#' 
#' # Get feature space for column "cellType" in metadata slot
#' 
#' object <- getFeatureSpace(object = object, pVar = "cellType")
#' 




getFeatureSpace <- function(object, pVar, varLim = 0.01, correction = "fdr", sig = 0.05){
  
  
  # Validations -------------------------------------------------------------
  
  if(!is(object, "scPred") & !is(object, "seurat")){
    stop("Invalid class for object: must be 'scPred' or 'seurat'")
  }
  
  if(!any(correction %in% stats::p.adjust.methods)){
    stop("Invalid multiple testing correction method. See ?p.adjust function")
  }
  
  if(is(object, "scPred")){
    classes <- metadata(object)[[pVar]]
  }else{
    classes <- object@meta.data[[pVar]]
  }
  
  if(is.null(classes)){
    stop("Prediction variable is not stored in metadata slot")
  }
  
  if(!is.factor(classes)){
    stop("Prediction variable must be a factor object")
  }else if(!all(levels(classes) %in% unique(classes))){
    stop("Not all levels are included in prediction variable")
  }else if(length(levels(classes)) == 1){
    stop("No training is possible with only one classification class. Check prediction variable")
  }
  
  
  # Filter principal components by variance ---------------------------------
  
  if(is(object, "scPred")){ # scPred object
    
    # Get PCA
    i <- object@expVar > varLim
    pca <- getPCA(object)[,i]
    
    # Get variance explained
    expVar <- object@expVar
    
  }else{ # seurat object
    
    # Get PCA
    i <- object@dr[["pca"]]@sdev > varLim
    pca <- object@dr[["pca"]]@cell.embeddings[,i]
    
    # Get variance explained
    expVar <- object@dr[["pca"]]@sdev**2/sum(object@dr[["pca"]]@sdev**2)
    names(expVar) <- colnames(pca)
    
    # Create svd slot for scPred object
    svd <- list(x = pca, 
                rotation = object@dr$pca@gene.loadings, 
                sdev = object@dr$pca@sdev, 
                center = rowMeans(object@scale.data), 
                scale = apply(object@scale.data, 1, sd))
    
    # Create scPred object
    object <- new("scPred", 
                  svd = svd, 
                  metadata = object@meta.data,
                  expVar = expVar, 
                  pseudo = FALSE, 
                  trainData = object@scale.data)
  }
  
  
  if(length(levels(classes)) == 2){
    message("First factor level in '", pVar, "' metadata column considered as positive class")
    res <- .getPcByClass(levels(classes)[1], object, classes, pca, correction, sig)
    res <- list(res)
    names(res) <- levels(classes)[1]
  }else{
    res <- lapply(levels(classes), .getPcByClass, object, classes, pca, correction, sig)
    names(res) <- levels(classes)
  }
  
  object@features <- res
  object@pVar <- pVar
  
  message("DONE!")
  object
  
}

.getPcByClass <- function(positiveClass, object, classes, pca, correction, sig){
  
  i <- classes != positiveClass
  newClasses <- as.character(classes)
  newClasses[i] <- "other"
  newClasses <- factor(newClasses, levels = c(positiveClass, "other"))
  
  
  apply(pca, 2, function(pc) wilcox.test(pc[newClasses == positiveClass], pc[newClasses == "other"])) %>% 
    lapply('[[', "p.value") %>% 
    as.data.frame() %>% 
    gather(key = "PC", value = "pValue") %>% 
    mutate(pValueAdj = p.adjust(pValue, method = correction, n = nrow(.))) %>% 
    arrange(pValueAdj) %>% 
    filter(pValueAdj < sig) %>% 
    mutate(expVar = object@expVar[match(PC, names(object@expVar))]) %>% 
    mutate(PC = factor(PC, levels = PC), cumExpVar = cumsum(expVar)) -> sigPCs
  
  sigPCs
}

