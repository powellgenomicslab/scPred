#' @title Definition of 'scPred' class
#' @description An S4 class to contain principal component analysis of a gene expression matrix, metadata, training, and
#' prediction information.
#' @slot pVar Column name from metadata to use as the variable to predict using
#' @slot reduction Dimensionality reduction
#' the informative principal components
#' @slot features A data frame with the following information:
#' \itemize{
#' \item PC: Principal component
#' \item Freq: Frequency of occurencxe of the principal component over a number of random samples from the PCA matrix
#' \item expVar: Explained variance by the principal component
#' \item cumExpVar: All principal components are ranked accoriding to their frequency of ocurrence and their variance explained. 
#' This column contains the cumulative variance explained across the ranked principal components
#' }
#' @slot train A list with all trained models using the \code{caret} package. Each model correspond to a cell type
#' @name scPred
#' @rdname scPred
#' @aliases scPred-class
#' @exportClass scPred
#' 


setClass("scPred", representation(pVar = "character",
                                  features = "list",
                                  reduction = "character",
                                  train = "list"),
         prototype(pVar = character(),
                   reduction = character(),
                   features = data.frame(),
                   train = list()))






