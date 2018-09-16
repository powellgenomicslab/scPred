#' @title show
#' @description Generic display function for \linkS4class{scPred} objects. Displays summary of the
#' object such as number of cells, genes, significant features.
#' @importFrom methods setMethod
#' @export

setMethod("show", signature("scPred"), function(object) {
  
  cat("'scPred' object\n")
  
  cat("\n- Expression data\n")
  nCells <- nrow(object@svd$x)
  nPCs <- ncol(object@svd$x)
  
  cat(sprintf("      Cell embeddings =  %i\n", nCells))
  cat(sprintf("      Gene loadings =  %i\n", nrow(getLoadings(object))))
  cat(sprintf("      PCs =  %i\n", nPCs))
  
  
  if(length(object@metadata) > 0){
    cat("\n- Metadata information\n")
    cat(sprintf("      %s\n", paste0(colnames(object@metadata), collapse = ", ")))
    
    if(length(object@pVar) != 0){
      cat(sprintf("      Prediction variable = %s\n", object@pVar))
      object@metadata[[object@pVar]] %>% 
        table() %>% 
        as.data.frame() %>% 
        column_to_rownames(".") %>% 
        set_colnames("Counts") %>% 
        print()
    }
  }
  

  
  if(length(object@features) != 0){
    
    cat("\n- Informative PCs per class\n")
    object@features %>% 
    sapply(nrow) %>% 
    as.data.frame() %>% 
    set_colnames("Features") %>% 
    print()
    
  }
  
  if(length(object@train) != 0){
    
    cat("\n- Training information\n")
    cat(sprintf("      Model: %s\n", object@train[[1]]$modelInfo$label))
    
    data.frame(names(object@train))
    
    getMetrics <- function(x, roc = TRUE){
      if(roc){
        round(x$results[bestModelIndex,c("ROC", "Sens", "Spec")], 3)
      }else{
        round(x$results[bestModelIndex,c("Accuracy", "Kappa")], 3)
      }
    }
    
    if(object@train[[1]]$metric == "ROC"){
      print(t(sapply(object@train, getMetrics)))
    }else{
      print(t(sapply(object@train, getMetrics, roc = FALSE)))
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
  return(object@svd$x)
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
  return(object@svd$rotation)
})


#' @title Plot PCA
#' @description Plots PCA with raw results or group by a variable
#' @param object \code{scPred} object
#' @param group Groupping variable included in metadata
#' @param pc Vector of length 2 with the dimensions to plot
#' @param predGroup Equivalent of \code{group} parameter for prediction dataset. Projection and predMeta slots must be filled
#' @param geom ggplot geom
#' @importFrom methods setMethod
#' @export

setGeneric("plotEigen", def = function(object, 
                                       group = NULL, 
                                       pc = c(1,2), 
                                       predGroup = NULL, 
                                       geom = c("points", "density_2d", "both")) {
  standardGeneric("plotEigen")
})

#' @title Plot PCA
#' @description Plots PCA 
#' @param object \code{scPred} object
#' @param group Groupping variable included in metadata
#' @param pc Vector of length 2 with the dimensions to plot
#' @param predGroup Equivalent of \code{group} parameter for prediction dataset. Projection and predMeta slots must be filled
#' @param geom ggplot geom
#' @importFrom methods setMethod
#' @export


setMethod("plotEigen", signature("scPred"), function(object, 
                                                     group = NULL, 
                                                     pc = c(1,2), 
                                                     predGroup = NULL,
                                                     geom = c("both", "points", "density_2d")){
  
  geom <- match.arg(geom)
  
  namesPC <- paste0("PC", pc)
  pca <- as.data.frame(subsetMatrix(getPCA(object), namesPC))
  
  pca$dataset <- "Train"
  
  # Check if a grouppping variable is provided
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
    
  }
  
  
  
  
  
  if(length(object@projection)){
    
    if(any(!namesPC %in% colnames(object@projection))){
      message("Performing projection of non-informative principal components...")
      pcaPred <- projectNewData(object, object@predData, informative = FALSE)
    }else{
      pcaPred <-  object@projection
    }
    
    pcaPred <- as.data.frame(subsetMatrix(pcaPred, namesPC))
    pcaPred$dataset <- "Prediction"
    
    if(!is.null(group)){
      if(!is.null(predGroup)){
        pcaPred[group] <- predGroup
        
      }else{
        pcaPred[group] <- "Unknown" 
      }
    }
    
    
    pcaAll <- rbind(pca, pcaPred)
    pcaAll$dataset <- factor(pcaAll$dataset, levels = c("Train", "Prediction"))
    
    
  }else{
    pcaAll <- pca
  }
  
  
  
  if(!is.null(group)){
    p <- ggplot(pcaAll, aes_string(x = namesPC[1], y = namesPC[2], color = group))
    
  }else{
    p <- ggplot(pcaAll, aes_string(x = namesPC[1], y = namesPC[2]))
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
  
  if(any("Prediction" == pcaAll$dataset)){
    p <- p + facet_wrap(~dataset)
  }
  
  
  p
  
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
  if(!pc %in% seq_len(ncol(object@svd$rotation))){
    stop(paste0("Principal component does not exist. Min 1, Max ", ncol(object@svd$rotation)))
  }
  
  # Validate that number of genes to be plotted is valid
  if(n > nrow(object@svd$rotation) | n < 1){
    stop(paste0("Only ", nrow(object@svd$rotation), 
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
    theme(legend.position = "none")
  
  
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
    stop("No models have been trained")
  }
  
  
  if(ncol(object@train[[1]]$pred) == 0){
    stop('No training results were calculated. Set savePredictions = "final" and returnData = TRUE')
  }
  
  probs <- lapply(names(object@train), function(model) extractProb(object@train[model]))
  names(probs) <- names(object@train)
  probs
  
})


#' @title Get predictions
#' @description Gets prediction probabilities for each cell class
#' @importFrom methods setMethod
#' @export

setGeneric("getPredictions", def = function(object) {
  standardGeneric("getPredictions")
})

#' @title Get training probabilities
#' @description Gets prediction probabilities for each cell class
#' @importFrom methods setMethod
#' @export

setMethod("getPredictions", signature("scPred"), function(object){
  
  if(length(object@train) == 0){
    stop("No predictions have been performed")
  }
  
  return(object@predictions)
  
})

#' @title Plot gene expression data
#' @description Plots a PCA and projection with the gene expression values for a given gene
#' @param object \code{scPred} object
#' @param gene Gene id. Must match one of the row names in the gene expression matrix used as input
#' @param pc Vector of length 2 with the two dimensions to be plotted
#' @param low Color for low gene expression
#' @param high Color for high gene expression
#' @importFrom methods setMethod
#' @export

setGeneric("plotExp", def = function(object, gene, pc = c(1,2), low = "gray", high = "red") {
  standardGeneric("plotExp")
})

#' @title Plot gene expression data
#' @description Plots a PCA and projection with the gene expression values for a given gene
#' @param object \code{scPred} object
#' @param gene Gene id. Must match one of the row names in the gene expression matrix used as input
#' @param pc Vector of length 2 with the two dimensions to be plotted
#' @param low Color for low gene expression
#' @param high Color for high gene expression
#' @importFrom methods setMethod
#' @export

setMethod("plotExp", signature("scPred"), function(object, gene, pc = c(1,2), low = "gray", high = "red"){
  
  iTrain <- which(rownames(object@trainData) == gene)
  
  if(length(iTrain) == 0){
    stop("Gene not found in training data")
  }
  
  trainGenes <- t(object@trainData[iTrain, , drop = FALSE])
  
  namesPC <- paste0("PC", pc)
  pca <- as.data.frame(subsetMatrix(getPCA(object), namesPC))
  
  pca$gene <- trainGenes
  pca$dataset <- "Train"
  
  if(length(object@projection) & length(object@predData)){
    
    if(any(!namesPC %in% names(object@projection))){
      message("Performing projection of non-informative principal components...")
      projection <- projectNewData(object, object@predData, informative = FALSE)
      object@projection <- projection
    }
    
    iPred <- which(rownames(object@predData) == gene)
    
    if(length(iPred) == 0){
      stop("Gene not found in prediction data")
    }
    
    
    predGenes <- t(object@predData[iPred,])
    pcaPred <- data.frame(object@projection[namesPC], gene = predGenes)
    pcaPred$dataset <- "Prediction"
    
    pcaAll <- rbind(pca, pcaPred)
    pcaAll$dataset <- factor(pcaAll$dataset, levels = c("Train", "Prediction"))
    
    
  }else{
    pcaAll <- pca
  }
  
  
  
  ggplot(pcaAll, aes_string(x = namesPC[1], y =  namesPC[2])) +
    geom_point(aes(color = gene)) +
    scale_color_gradient(low = low, high = high) +
    ggtitle(gene) +
    theme_bw() -> p
  
  if(any("Prediction" == pcaAll$dataset)){
    p <- p + facet_wrap(~dataset)
  }
  
  p
  
})


