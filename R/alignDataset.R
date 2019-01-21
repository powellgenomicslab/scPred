#' @title Align train and test datasets
#' @description Aligns train and test datasets in a low dimensional space using the 
#' dynamic time warping algorithm implemented in Seurat
#' @param object \code{scPred} object
#' @param prediction A seurat object containing the prediction dataset. \code{data} slot must be filled with the 
#' after running NormalizeData() function
#' @return A matrix with the aligned projection
#' @author José Alquicira Hernández
#' @importFrom Seurat SetDimReduction
#' @export

alignDataset <- function(object, prediction){


  # Validations -------------------------------------------------------------

  ## Training object type
  if(!is(object, "scPred")){
    stop("Invalid class for 'training' object: must be 'scPred'")
  }
  ## Prediction object type
  if(!is(prediction, "seurat")){
    stop("Invalid class for 'prediction' object: must be 'seurat'")
  }
  

  # Set up cell embeddings and training loadings ----------------------------


  # Perform projection
  ## We perform a Seurat normalization and data scaling before the projection and return both
  ## the loadings and the projection
  
  res <- projectNewData(object, 
                        newData = as.matrix(prediction@data), 
                        seurat = TRUE, 
                        returnLoadings = TRUE)
  
  
  # Subset/intersect prediction data and training loadings based on gene labels
  predictionData <- as.matrix(prediction@data[match(rownames(res$loadings), rownames(prediction@data)),])
  # Subset/intersect prediction and training data based on gene labels
  trainingData <- as.matrix(object@trainData[match(rownames(predictionData), rownames(object@trainData)),])
  ## Training and prediction data contain the same genes in the same order

  
  # Merge train and prediction data
  trainPred <- cbind(trainingData,  predictionData)
  
  # Create seurat object with merged data
  align <- CreateSeuratObject(trainPred)
  align@raw.data <- NULL
  
  # Add dataset origin info to metadata
  batch <- as.factor(c(rep("train", ncol(trainingData)), 
                       rep("pred", ncol(predictionData))))
  align@meta.data$dataset_dummy <- batch
  
  
  # Extract significant PCS from scores/cell embeddings for training 
  # and prediction datasets
  eigenSpace <- subsetMatrix(object@svd$x, colnames(res$proj))
  # Merge all cell embeddings
  embeddings <- rbind(eigenSpace, res$proj)
  
  # Are genes in loadings and data intersection the same and ordered?
  # all(rownames(res$loadings) == rownames(trainPred))
  

  # Create “dim.reduction” object -------------------------------------------
  new.type <- "pca.scpred"
  
  # Set dimension reduction slot
  align <- SetDimReduction(object = align, 
                                     reduction.type = new.type, 
                                     slot = "cell.embeddings", 
                                     new.data = embeddings)
  align <- SetDimReduction(object = align, 
                                    reduction.type = new.type, 
                                    slot = "gene.loadings", 
                                    new.data = res$loadings)
  
  align <- SetDimReduction(object = align, reduction.type = new.type, 
                                     slot = "key", new.data = "PC")
  
  # Add merged data to align object
  align@data <- trainPred
  
  # Add normalization parameters
  align@calc.params$NormalizeData <- list(assay.type = "RNA", 
                                           normalization.method = "LogNormalize", 
                                           scale.factor = 10000)
  
  # Align features
  align <- .alignSubspaceSeurat(align,
                               dims.align = seq_len(ncol(embeddings)))

  
  # Extract projection
  i <- align@meta.data$dataset_dummy == "pred"
  alignedProjection <- align@dr$pca.scpred.aligned@cell.embeddings[i,]
  dimensionNames <- colnames(embeddings)
  colnames(alignedProjection) <- dimensionNames
  
  alignedProjection
  
}
