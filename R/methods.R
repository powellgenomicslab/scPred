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
      
      if(object@train[[i]]$metric == "ROC"){
        metrics <- round(object@train[[i]]$results[bestModelIndex,c("ROC", "Sens", "Spec")], 3)
        cat(sprintf("      AUROC = %s, Sensitivity = %s, Specificity = %s\n", 
                    metrics$ROC, metrics$Sens, metrics$Spec))
      }else{
        metrics <- round(object@train[[i]]$results[bestModelIndex,c("Accuracy", "Kappa")], 3)
        cat(sprintf("      Accuracy = %s, Kappa = %s\n", 
                    metrics$Accuracy, metrics$Kappa))
        
        
      }
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
  
  
  # Validations -------------------------------------------------------------
  
  # Validate class object
  if(!(is.numeric(pc) | is.integer(pc))){
    stop("'pc' parameter values must be numeric or integer")
  }
  
  # Validate dimensions to be plotted
  if(length(pc) != 2){
    stop("Only two principal components can be plotted")
  }
  if(!all(pc %in% seq_len(ncol(object@pca$x)))){
    stop(paste0("Principal component(s) does not exist. Min 1, Max ", ncol(object@pca$x)))
  }
  
  # Validate valid geom
  geom <- match.arg(geom)
  
  
  # Main function -----------------------------------------------------------
  
  # Obtain scores
  pca <- getPCA(object)[,pc] 
  
  if(!is.null(group)){
    
    # Validate if provided grouping variable is valid
    
    ## Is metadata data.frame available?
    if(ncol(object@metadata) == 0){
      stop("No metadata has been assigned to 'scPred' object")
    }
    
    ## Present in metadata?
    if(!any(group %in% names(object@metadata))){
      stop("'group' variable is not included in metadata")
    }
    
    ## Is grouping variable a factor?
    if(!is(object@metadata[[group]], "factor")){
      stop("'group' variable must be a factor variable")
    }
    
    ## Obtain grouping variable and integrate to scores data.frame
    metadata <- object@metadata[group]
    pca <- cbind(pca, metadata)
    
    # Set plotting elements
    p <- ggplot(pca, aes_string(x = names(pca)[1], y = names(pca)[2], color = group))
    
  }else{
    # Set plotting elements if no grouping variable is provided
    p <- ggplot(pca, aes_string(x = names(pca)[1], y = names(pca)[2]))
  }
  
  # Add layers depending of geom_ provided
  if(geom == "points" | geom == "both"){
    p <- p + geom_point()
  }
  if(geom == "density_2d" | geom == "both"){
    p <- p +  geom_density_2d() 
  }
  
  # Add color palette and swith theme
  p <- p + 
    scale_color_brewer(palette = "Set1") +
    theme_bw()
  
  # Add marginal density plots if specified
  if(marginal){
    # Plot grouped marginal plot if grouping variable was provided
    if(!is.null(group)){
      ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE)
    }else{
      ggMarginal(p, type = "density")
    }
    
  }else{
    p
  }
  
})


#' @title Get training probabilities
#' @description Gets training probabilities for each trained model
#' @importFrom methods setMethod
#' @export

setGeneric("getTrainResults", def = function(object) {
  standardGeneric("getTrainResults")
})

#' @title Get training probabilities
#' @description Gets training probabilities for each trained model
#' @importFrom methods setMethod
#' @export

setMethod("getTrainResults", signature("scPred"), function(object){
  
  if(length(object@train) == 0){
    stop("No model has been trained")
  }
  
  
  if(ncol(object@train[[1]]$pred) == 0){
    stop('No training result was calculated. Set savePredictions = "final" and returnData = TRUE')
  }
  
  probs <- lapply(names(object@train), function(model) extractProb(object@train[model]))
  names(probs) <- names(object@train)
  probs
  
})


#' @title Plot loadings
#' @description Plot loading values for any given principal component in a \code{scPred} object
#' @param object A \code{scPred} object
#' @param pc The number of the principal component to be plotted
#' @param n Top `n` variable genes to plot. Notice that the number of genes plotted is n*2 as both 
#' negative and positive loadings are considered
#' @importFrom methods setMethod
#' @export

setGeneric("plotLoadings", def = function(object, pc = 1, n = 10) {
  standardGeneric("plotLoadings")
})

#' @title Plot loadings
#' @description Plot loading values for any given principal component in a \code{scPred} object
#' @param object A \code{scPred} object
#' @param pc The number of the principal component to be plotted
#' @param n Top `n` variable genes to plot. Notice that the number of plotted genes is n*2 as both 
#' negative and positive loadings are considered
#' @return A ggplot2 object
#' @importFrom methods setMethod
#' @export

setMethod("plotLoadings", signature("scPred"), function(object, pc = 1, n = 10){
  
  # Validations -------------------------------------------------------------
  
  # Validate class object for `pc`
  if(!(is.numeric(pc) | is.integer(pc))){
    stop("`pc` parameter value must be numeric or integer")
  }
  
  # Validate class object for `n`
  if(!(is.numeric(n) | is.integer(n))){
    stop("`n` parameter values must be numeric or integer")
  }
  
  # Validate that only one principal component to be plotted was provided
  if(length(pc) != 1){
    stop("Only one principal component can be plotted. Provide a single number to `pc` parameter")
  }
  
  # Validate that only one principal component to be plotted was provided
  if(length(n) != 1){
    stop("`n` must be a scalar integer")
  }
  
  # Validate the principal component provided is valid
  if(!pc %in% seq_len(ncol(object@pca$rotation))){
    stop(paste0("Principal component does not exist. Min 1, Max ", ncol(object@pca$rotation)))
  }
  
  # Validate that number of genes to be plotted is valid
  if(n > nrow(object@pca$rotation) | n < 1){
    stop(paste0("Only ", nrow(object@pca$rotation), 
                " genes are included in the loadings matrix. Check provided number to `n` parameter"))
  }
  
  
  # Main function -----------------------------------------------------------
  
  # Create column variable label for principal component
  pc <- paste0("PC", pc)
  
  # Obtain loadings, selects and orders genes according to their loading values
  object %>% 
    getLoadings() %>% 
    as.data.frame() %>% 
    rownames_to_column("gene") %>% 
    select(gene, !! sym(pc)) %>% 
    arrange(!! sym(pc)) -> scores
  
  # Get more variable genes (absolute value)
  top_positive <- top_n(scores, n = n, !! sym(pc))
  top_negative <- top_n(scores, n = -n, !! sym(pc))
  
  # Merge variable genes and sets factor variable to set the order of the genes when plotted
  rbind(top_negative, top_positive) %>% 
    mutate(gene = factor(gene, levels = gene)) %>% 
    mutate(direction = as.factor(c(rep("negative", n), rep("positive", n)))) -> top
  
  # Plot "lollipop" graph
  top %>%  
    ggplot() +
    aes_string(x = "gene", y = pc, color = "direction") +
    xlab("Genes") +
    geom_point() +
    geom_segment(aes_string(xend = "gene", yend = mean(top[[pc]]))) +
    coord_flip() +
    scale_color_brewer(palette = "Set1") +
    theme_bw() +
    theme(legend.position="none")
  
  
})

