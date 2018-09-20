#' @title Align train and test datasets
#' @description Aligns train and test datasets in a low dimensional space using the dynamic time warping alogorithm implemented in Seurat



alignDataset <- function(object, prediction){
  z -> object
  b -> prediction

  # Validations -------------------------------------------------------------

  ## Training object type
  if(!is(object, "scPred")){
    stop("Invalid class for 'training' object: must be 'scPred'")
  }
  ## Prediction object type
  if(!is(prediction, "seurat")){
    stop("Invalid class for 'prediction' object: must be 'seurat'")
  }
  
  
  ## Check if a PCA has been computed
  if(!("pca" %in% names(prediction@dr))){
    stop("No PCA has been computet yet for prediction dataset. See RunPCA() function")
  }
  
  

  # Set up cell embeddings and training loadings ----------------------------


  # Perform projection
  res <- projectNewData(object, newData = as.matrix(prediction@data), seurat = TRUE, returnLoadings = TRUE)
  
  # Imtersect prediction data and loadings based on genes
  predictionData <- as.matrix(prediction@data[match(rownames(res$loadings), rownames(prediction@data)),])
  
  
  # Intersect training and prediction data data by genes
  trainingData <- as.matrix(object@trainData[match(rownames(predictionData), rownames(object@trainData)),])
  

  # Aggregate train and prediction data
  trainPred <- cbind(trainingData,  predictionData)
  
  # Create seurat object with aggregated data
  align <- CreateSeuratObject(trainPred)
  align@raw.data <- NULL
  
  # Add dataset origin info to metadata
  
  batch <- as.factor(c(rep("train", ncol(trainingData)), rep("pred", ncol(predictionData))))
  align@meta.data$dataset_dummy <- batch
  
  
  # Set scores/cell embeddings for training and prediction 
  eigenSpace <- subsetMatrix(object@svd$x, colnames(res$proj))
  embeddings <- rbind(eigenSpace, res$proj)
  
  
  # Are genes in loadings and data intersection the same and ordered?
  # all(rownames(res$loadings) == rownames(trainPred))
  

  # Create “dim.reduction” object -------------------------------------------
  new.type <- "pca.scpred"
  align <- Seurat:::SetDimReduction(object = align, 
                                     reduction.type = new.type, 
                                     slot = "cell.embeddings", 
                                     new.data = embeddings)
  align <- Seurat:::SetDimReduction(object = align, 
                                    reduction.type = new.type, 
                                    slot = "gene.loadings", 
                                    new.data = res$loadings)
  
  align <- Seurat:::SetDimReduction(object = align, reduction.type = new.type, 
                                     slot = "key", new.data = "PC")
  align@data <- trainPred
  align@calc.params$NormalizeData <- list(assay.type = "RNA", 
                                           normalization.method = "LogNormalize", 
                                           scale.factor = 10000)
  
  # Align features

  align <- AlignSubspaceSeurat(align, 
                               reduction.type = "pca.scpred", 
                               grouping.var = "dataset_dummy",
                               dims.align = seq_len(ncol(embeddings)))

  
  # Extract projection
  
  i <- align@meta.data$dataset_dummy == "pred"
  alignedProjection <- align@dr$pca.scpred.aligned@cell.embeddings[i,]
  dimensionNames <- colnames(embeddings)
  colnames(alignedProjection) <- dimensionNames
  
  alignedProjection
  
}
