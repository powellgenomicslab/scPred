#' @title Gets contingency table
#' @description Creates a cross table using two columns from the metadata
#' @param object \code{Seurat} object
#' @param true Column name in \code{meta.data} slot that corresponds to the true known classes
#' @param pred Column name in \code{meta.data} slot that corresponds to the predicted classes
#' if they  have been assigned independently from the \code{scPredict()} function
#' @param output Return counts, fraction, or proportions? Default: counts
#' @param digits If proportions are returned, number of digits to round numbers
#' @return A contingency table
#' @export
#' @author Jose Alquicira Hernandez
#'


setGeneric("crossTab", def = function(object, true = NULL, pred = NULL, output = c("counts", "fraction", "prop"), digits = 2) {
  standardGeneric("crossTab")
})


#' @title Get scPred object
#' @description Accessor function to retrieve scPred models from Seurat objects
#' @param object Seurat object
#' @return scPred object
#' @export


setGeneric("get_scpred", def = function(object) {
  standardGeneric("get_scpred")
})

#' @title Get classification models
#' @description Accessor function to retrieve classification models from Seurat 
#' and scPred objects
#' @param object \code{Seurat} or scPred object
#' @return A list of \code{train} objects


setGeneric("get_classifiers", def = function(object) {
  standardGeneric("get_classifiers")
})

#' @title Plot training probabilities
#' @description Plots all training probabilities for each cell type
#' @param object Seurat or scPred object
#' @param size Point size for each cell
#' @export
#' @return Plot with the probability distribution for each cell type


setGeneric("plot_probabilities", def = function(object, size = 0.8) {
  standardGeneric("plot_probabilities")
})


#' @title Get training probabilities
#' @description Gets training probabilities for each cell type
#' @param object Seurat or scPred object
#' @export
#' @return A data frame with all cell-type probabilities associated to each cell

setGeneric("get_probabilities", def = function(object) {
  standardGeneric("get_probabilities")
})



#' @title Get metadata from scPred object
#' @description Accessor function to retrieve metadata from scPred object
#' @return A dataframe including the cell barcodes and prediction variable 
#' (cell type labels)
#' @export


setGeneric("get_metadata", def = function(object) {
  standardGeneric("get_metadata")
})