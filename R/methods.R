#' @title show
#' @description Generic display function for \linkS4class{scPred} objects. Displays summary of the
#' object such as number of cells, genes, significant features.
#' @importFrom methods setMethod
#' @importFrom tibble column_to_rownames
#' @importFrom magrittr set_colnames
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
        set_colnames("n") %>% 
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
      
<<<<<<< HEAD
      if(object@train[[i]]$metric == "ROC"){
        metrics <- round(object@train[[i]]$results[bestModelIndex,c("ROC", "Sens", "Spec")], 3)
        cat(sprintf("      AUROC = %s, Sensitivity = %s, Specificity = %s\n", 
                    metrics$ROC, metrics$Sens, metrics$Spec))
=======
      bestModelIndex <- as.integer(rownames(x$bestTune))
      
      if(roc){
        round(x$results[bestModelIndex,c("ROC", "Sens", "Spec")], 3)
>>>>>>> development
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
                                       pc1 = 1,
                                       pc2 = 2,
                                       predGroup = NULL, 
                                       geom = c("points", "density_2d", "both"),
                                       plotPred = TRUE) {
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
                                                     pc1 = 1, 
                                                     pc2 = 2,
                                                     predGroup = NULL,
                                                     geom = c("both", "points", "density_2d"),
                                                     plotPred = TRUE){
  
  
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
  
  pc <- c(pc1, pc2)
  
<<<<<<< HEAD
  # Main function -----------------------------------------------------------
  
  # Obtain scores
  pca <- getPCA(object)[,pc] 
  
=======
  namesPC <- paste0("PC", pc)
  
  pca <- getPCA(object)
  allpcs <- colnames(pca)                    
  
  pcsInEigen <- namesPC %in% allpcs
  if(!all(pcsInEigen)){
    label <- paste0(namesPC[!pcsInEigen], collapse = " and ")
    stop(label, " not present in SVD results")
  }
                      
  pca <- as.data.frame(subsetMatrix(pca, namesPC))
  
  pca
  
  
  
  pca$dataset <- "Train"
  
  # Check if a grouppping variable is provided
>>>>>>> development
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
      metadata <- object@metadata[[group]]
    }else{
      metadata <- as.factor(as.character(object@metadata[[group]]))
    }
    
<<<<<<< HEAD
    ## Obtain grouping variable and integrate to scores data.frame
    metadata <- object@metadata[group]
    pca <- cbind(pca, metadata)
    
    # Set plotting elements
    p <- ggplot(pca, aes_string(x = names(pca)[1], y = names(pca)[2], color = group))
    
  }else{
    # Set plotting elements if no grouping variable is provided
    p <- ggplot(pca, aes_string(x = names(pca)[1], y = names(pca)[2]))
=======
    pca <- cbind(pca, metadata)
    names(pca)[4] <- group 
  }
  
  
  
  
  
  if(length(object@projection)){
    
    if(any(!namesPC %in% colnames(object@projection))){
      message("Performing projection of non-informative principal components...")
      pcaPred <- projectNewData(object, object@predData, informative = FALSE)
    }else{
      pcaPred <- object@projection
    }
    
    pcaPred <- as.data.frame(subsetMatrix(pcaPred, namesPC))
    pcaPred$dataset <- "Prediction"
    
    if(!is.null(group)){
      if(!is.null(predGroup)){
        predVar <- object@predMeta[[predGroup]]
        pcaPred[group] <- factor(predVar, levels = unique(predVar))
        
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
>>>>>>> development
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
  
<<<<<<< HEAD
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
=======
  if(any("Prediction" == pcaAll$dataset) & plotPred){
    p <- p + facet_wrap(~dataset)
>>>>>>> development
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
    stop("No model has been trained")
  }
  
  
  if(ncol(object@train[[1]]$pred) == 0){
    stop('No training result was calculated. Set savePredictions = "final" and returnData = TRUE')
  }
  
  probs <- lapply(names(object@train), function(model) extractProb(object@train[model]))
  names(probs) <- names(object@train)
  probs
  
})

<<<<<<< HEAD

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
=======

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




#' @title Get accuracy
#' @description Returns the accuracy per group (recall) given a known "true" for the prediction dataset
#' @param object \code{scPred} object
#' @param var Variable in \code{predMeta} slot containing them true classes. True classes will be compared to
#' the classifications provided by scPred in the \code{predictions} slot
#' @importFrom methods setMethod
#' @export

setGeneric("getAccuracy", def = function(object, var) {
  standardGeneric("getAccuracy")
})

#' @title Get accuracy
#' @description Returns the accuracy per group (recall) given a known "true" for the prediction dataset
#' @param object \code{scPred} object
#' @param var Variable in \code{predMeta} slot containing them true classes. True classes will be compared to
#' the classifications provided by scPred in the \code{predictions} slot
#' @importFrom methods setMethod
#' @export

setMethod("getAccuracy", signature("scPred"), function(object, var){
  
  if(length(object@predMeta) == 0){
    stop("No metadata for prediction dataset has been stored")
  }
  
  if(!var %in% names(object@predMeta)){
    stop("Variable not present in metadata")
  }
  
  response <- make.names(as.character(object@predMeta[[var]]))
  predictions <- object@predictions["predClass"]
  
  truePred <- cbind(predictions, response)
  
  # Get accuracy
  
  res <- as.data.frame(table(truePred$response))
  names(res) <- c("true", "total")
  
  truePred %>% 
    set_colnames(c("prediction", "true")) %>% 
    mutate(result = if_else(prediction == true, "correct", "incorrect")) %>% 
    group_by(true, result) %>% 
    summarise(n = n()) %>% 
    filter(result == "correct") %>% 
    select(-result) -> counts
  
  
  left_join(res, counts, by = "true") %>% 
    select(true, n, total) %>% 
    mutate(n = if_else(is.na(n), 0, as.numeric(n))) %>% 
    mutate(accuracy = n/total) %>% 
    column_to_rownames("true")
  
})

#' @title Gets contingency table
#' @description If the \code{predMeta} slot contains a column with the true classes of the cells,
#' builds a contingency table by using this column as reference and comparing it to the predicted classes
#' obtained with \code{scPredict}. If an independent column with predicted classes is in the prediction
#' metadata, this column instead of the default classes assigned by \code{scPred} can be provided using the
#' \code{pred} parameter.
#' @param object \code{scPred} object
#' @param true Column name in \code{predMeta} slot that corresponds to the true known classes 
#' @param pred Column name in \code{predMeta} slot that corresponds to the predicted classes
#' if they  have been assigned independently from the \code{scPredict()} function
#' @param fill Value to fill contingency table ff unique cell classes from the true and the 
#' predicted columns do not match.
#' @param prop Return proportions or counts? Default: proportions
#' @param digits If proportions are returned, number of digits to round numbers
#' @return A contingency table 
#' @export
#' @importFrom dplyr group_by_ summarise
#' @importFrom tidyr spread
#' @importFrom magrittr "%>%"
#' @importFrom tibble column_to_rownames
#' @author José Alquicira Hernández
#' 


setGeneric("crossTab", def = function(object, true, pred = NULL, fill = 0, prop = TRUE, digits = 2) {
  standardGeneric("crossTab")
})

#' @title Gets contingency table
#' @description If the \code{predMeta} slot contains a column with the true classes of the cells,
#' builds a contingency table by using this column as reference and comparing it to the predicted classes
#' obtained with \code{scPredict}. If an independent column with predicted classes is in the prediction
#' metadata, this column instead of the default classes assigned by \code{scPred} can be provided using the
#' \code{pred} parameter.
#' @param object \code{scPred} object
#' @param true Column name in \code{predMeta} slot that corresponds to the true known classes 
#' @param pred Column name in \code{predMeta} slot that corresponds to the predicted classes
#' if they  have been assigned independently from the \code{scPredict()} function
#' @param fill Value to fill contingency table ff unique cell classes from the true and the 
#' predicted columns do not match.
#' @param prop Return proportions or counts? Default: proportions
#' @param digits If proportions are returned, number of digits to round numbers
#' @return A contingency table 
#' @export
#' @importFrom dplyr group_by_ summarise
#' @importFrom tidyr spread
#' @importFrom magrittr "%>%"
#' @importFrom tibble column_to_rownames
#' @author José Alquicira Hernández
#' 

setMethod("crossTab", signature("scPred"), 
          function(object, true, pred = NULL, fill = 0, prop = TRUE, digits = 2){
  
  res <- processPreds(object = object, true = true, pred = pred)
  predictions <- res$predictions
  true <- res$true
  pred <- res$pred
  
  predictions %>% 
    group_by_(pred, true) %>% 
    summarise(n = n()) %>% 
    spread(key = true, value = "n", fill = fill) %>% 
    as.data.frame() %>% 
    column_to_rownames(pred) -> x
  
  if(prop){
    row_names <- rownames(x)
    x <- mapply(function(x,d){x/d}, x, colSums(x))
    rownames(x) <- row_names
    x %>% 
      round(digits) %>% 
      as.data.frame() -> x
    
  }
  x
})



#' @title Plot prediction probabilities
#' @description Plots the probability distributions according to a grouping variable. 
#' @param object \code{scPred} object
#' @param facet Column name in \code{predMeta} slot
#' if they have been assigned independently from the \code{scPredict()} function
#' @return A distribution probability plot divided into facets representing the groups contained in the
#' provided facet column name and colors to the corresponding probabilities for each class. If no grouping
#' variable is provided to the \code{facet} parameter, distributions are divided according to the predicted classes
#' by default
#' @export
#' @author José Alquicira Hernández
#' 


setGeneric("plotPredProbs", def = function(object, facet = NULL){
  standardGeneric("plotPredProbs")
})



#' @title Plot prediction probabilities
#' @description Plots the probability distributions according to a grouping variable. 
#' @param object \code{scPred} object
#' @param facet Column name in \code{predMeta} slot
#' if they have been assigned independently from the \code{scPredict()} function
#' @return A distribution probability plot divided into facets representing the groups contained in the
#' provided facet column name and colors to the corresponding probabilities for each class. If no grouping
#' variable is provided to the \code{facet} parameter, distributions are divided according to the predicted classes
#' by default
#' @export
#' @author José Alquicira Hernández
#' 


setMethod("plotPredProbs", signature("scPred"), function(object, facet = NULL) {
  
  if(nrow(object@predMeta) == 0){
    stop("No prediction metadata has been assigned to scPred object")
  }
  
  if(is.null(facet)){
    predictions <- getPredictions(object)
    facet <- "predClass"
    
  }else{
    res  <- processPreds(object, true = facet)
    predictions <- res$predictions
    fill <- res$true
  }
  
  n <- length(object@train)
  
  
  predictions %>% 
    gather(key = "Class", value = "prob", seq_len(n)) %>% 
    ggplot() +
    aes_string(x = "prob", fill = "Class") +
    geom_histogram(color = "black") +
    facet_wrap(as.formula(paste("~", facet)), scales = "free") +
    scale_fill_manual(values = getPalette(n)) +
    theme_bw() +
    xlab("Probability") +
    ylab("Number of cells")
>>>>>>> development
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

