#' @title show
#' @description Generic display function for \linkS4class{scPred} objects. Displays summary of the
#' object such as number of cells, genes, significant features.
#' @importFrom methods setMethod
#' @export

setMethod("show", signature("scPred"), function(object) {
  
  cat("'scPred' object\n")
  
  cat("- Expression data\n")
  nCells <- nrow(object@pca$x)
  nPCs <- ncol(object@pca$x)
  
  cat(sprintf("      Cells =  %i\n", nCells))
  cat(sprintf("      Genes =  %i\n", nrow(getLoadings(object))))
  cat(sprintf("      PCs =  %i\n", nPCs))
  
  
  if(length(object@metadata) > 0){
    cat("- Metadata information\n")
    cat(sprintf("      %s\n", paste0(colnames(object@metadata), collapse = ", ")))
  }
  
  if(length(object@pVar) != 0 & length(object@metadata) > 0){
    cat("- Prediction\n")
    cat(sprintf("      Variable = %s\n", object@pVar))
    
  }
  
  if(length(object@features) != 0){
    cat("- Informative PCs per class\n")
    for(i in seq_len(length(object@features))){
      cat(sprintf("      %s = %i\n", names(object@features)[i], nrow(object@features[[i]])))
    }
    
  }
  
  if(length(object@train) != 0){
    cat("- Training\n")
    cat(sprintf("      Model: %s\n", object@train[[1]]$modelInfo$label))
    for(i in seq_len(length(object@train))){
      cat(sprintf("      Class -> %s\n", names(object@train)[i]))
      bestModelIndex <- as.integer(rownames(object@train[[i]]$bestTune))
      metrics <- round(object@train[[i]]$results[bestModelIndex,c("ROC", "Sens", "Spec")], 3)
      cat(sprintf("      AUROC = %s, Sensitivity = %s, Specificity = %s\n", 
                  metrics$ROC, metrics$Sens, metrics$Spec))
    }
  }
  
  
})


#' @title Get metadata
#' @description Gets metadata from \code{scPred} object
#' @importFrom methods setMethod
#' @export


setGeneric("metadata", function(object) standardGeneric("metadata"))
setMethod("metadata", "scPred", function(object) {
  out <- object@metadata
  return(out)
})


#' @title Set metadata
#' @description Sets metadata to a \code{scPred} object. Metadata must be a dataframe
#' \itemize{
#' \item row names: ids matching the column names of the expression matrix
#' \item columns: associated metadata such as cell type, conditions, sample or, batch. 
#' } 
#' @importFrom methods setMethod
#' @importFrom methods loadMethod
#' @export

setGeneric("metadata<-", function(object, value) standardGeneric("metadata<-"))
setMethod("metadata<-", signature("scPred"), 
          function(object, value) {
            if(!all(rownames(getPCA(object)) == row.names(value))){
              stop("Cells Ids do not match cell IDs in metadata")
            }
            object@metadata <- value
            object
          }
)


#' @title Get principal components
#' @description Gets matrix with principal components from a \code{scPred} object
#' @importFrom methods setMethod
#' @export

setGeneric("getPCA", def = function(object) {
  standardGeneric("getPCA")
})


#' @title Get principal components
#' @description Gets matrix with principal components from a \code{scPred} object
#' @importFrom methods setMethod
#' @export

setMethod("getPCA", signature("scPred"), function(object) {
  return(as.data.frame(object@pca$x))
})


#' @title Get loadings matrix
#' @description Gets rotation matrix (right singular vectors) from a \code{scPred} object
#' @importFrom methods setMethod
#' @export

setGeneric("getLoadings", def = function(object) {
  standardGeneric("getLoadings")
})

#' @title Get loadings matrix
#' @description Gets rotation matrix (right singular vectors) from a \code{scPred} object
#' @importFrom methods setMethod
#' @export

setMethod("getLoadings", signature("scPred"), function(object) {
  return(object@pca$rotation)
})


#' @title Plot PCA
#' @description Plots PCA with raw results or group by a variable
#' @importFrom methods setMethod
#' @export

setGeneric("plotEigen", def = function(object, group = NULL, pc = c(1,2), geom = c("points", "density_2d", "both"), marginal = TRUE) {
  standardGeneric("plotEigen")
})

#' @title Plot PCA
#' @description Plots PCA 
#' @importFrom methods setMethod
#' @export


setMethod("plotEigen", signature("scPred"), function(object, group = NULL, pc = c(1,2), geom = c("both", "points", "density_2d"), marginal = TRUE){
  
  geom <- match.arg(geom)
  pca <- getPCA(object)[,pc] 
  
  
  if(!is.null(group)){
    
    if(ncol(object@metadata) == 0){
      stop("No metadata has been assigned to 'scPred' object")
    }
    
    if(!any(group %in% names(object@metadata))){
      stop("'group' variable is not included in metadata")
    }
    
    if(!is(object@metadata[[group]], "factor")){
      stop("'group' variable must be a factor variable")
    }
    
    metadata <- object@metadata[group]
    pca <- cbind(pca, metadata)
    p <- ggplot(pca, aes_string(x = names(pca)[1], y = names(pca)[2], color = group))
    
  }else{
    p <- ggplot(pca, aes_string(x = names(pca)[1], y = names(pca)[2]))
  }
  
  if(geom == "points" | geom == "both"){
    p <- p + geom_point()
  }
  if(geom == "density_2d" | geom == "both"){
    p <- p +  geom_density_2d() 
  }
  p <- p + 
    scale_color_brewer(palette = "Set1") +
    theme_bw()
  if(marginal){
    if(!is.null(group)){
    ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE)
    }else{
      ggMarginal(p, type = "density")
    }
  }else{
    p
  }
  
})

