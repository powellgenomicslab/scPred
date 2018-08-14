setOldClass("irlba_prcomp")

#' @title Definition of 'eigenPred' class
#' @description An S4 class to contain principal component analysis of a gene expression matrix, metadata, training, and
#' prediction information.
#' @slot pca Principal component analysis performed with \code{prcomp_irlba()} function
#' @slot metadata A dataframe with:
#' \itemize{
#' \item row names: ids matching the column names of the gene expression matrix
#' \item columns: associated metadata such as cell type, conditions, sample or batch. 
#' }
#' @slot expVar Explained variance by each principoal component 
#' @slot pVar Column name from metadata to use as the variable to predict using
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
#' @slot pseudo TRUE of a \code{log2(data + 1)} transformation was performed before performing the PCA 
#' @name scPred
#' @rdname scPred
#' @aliases scPred-class
#' @exportClass scPred
#' 


setClass("scPred", representation(pca = "irlba_prcomp",
                                    metadata = "data.frame",
                                    trainData = "matrix",
                                    predData = "matrix",
                                    expVar = "numeric",
                                    pVar = "character",
                                    features = "list",
                                    train = "list",
                                    pseudo = "logical"),
         prototype(metadata = data.frame(),
                   features = data.frame()))






