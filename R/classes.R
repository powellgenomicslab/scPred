#' @title Definition of 'scPred' class
#' @description An S4 class to containing features, dimensionality reduction information, and trained models
#' @slot pvar Column name from metadata to use as the variable to predict using
#' @slot metadata A data frame to store metadata including prediction variable (pvar)
#' @slot features A data frame with the following information:
#' \itemize{
#' \item feature: Eigenvector (e.g. principal component)
#' \item pValue: Significance value from a wilcoxon test 
#' \item pValueAdj: Adjusted p-value for multiple testing
#' }
#' @slot loadings Gene loadings
#' @slot scaling Means and standard deviation to center and standardize data
#' @slot reduction Dimensionality reduction name
#' @slot reduction_key Dimensionality reduction name key
#' @slot train A list with all trained models using the \code{caret} package. Each model correspond to a cell type
#' @slot mist A list to store extra information and for developing testing
#' @name scPred
#' @rdname scPred
#' @aliases scPred-class
#' @exportClass scPred
#' 


setClass("scPred", representation(pvar = "character",
                                  metadata = "data.frame",
                                  features = "list",
                                  cell_embeddings = "matrix",
                                  feature_loadings = "matrix",
                                  scaling = "data.frame",
                                  reduction = "character",
                                  reduction_key = "character",
                                  train = "list",
                                  misc = "list"),
         prototype(pvar = character(),
                   metadata = data.frame(),
                   features = list(),
                   cell_embeddings = matrix(),
                   feature_loadings = matrix(),
                   scaling = data.frame(),
                   reduction = character(),
                   reduction_key = character(),
                   train = list(),
                   misc = list())
         )






