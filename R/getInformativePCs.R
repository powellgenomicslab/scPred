#' @title Get informative principal components
#' @description Given a prediction variable, finds a set of class-informative principal components that are significant across n subsets of cells from the 
#' PCA matrix. A Wilcoxon rank sum test is used to determine a difference between the value distributions of classes from the prediction variable.
#' @param object An \code{eigenPred} object
#' @param pVar Prediction variable corresponding to a column in \code{metadata} slot
#' @param varLim Threshold to filter principal components based on variance explained.
#' @param correction Multiple testing correction method. Default: false discovery rate. See \code{p.adjust} function 
#' @param sig Significance level to determine principal components explaining class identity
#' @return A \code{eigenPred} object with two additional filled slots:
#' \itemize{
#' \item \code{features}: A data frame with significant principal components the following information:
#' \itemize{
#' \item PC: Principal component
#' \item pValue: p-value obtained from Mann-Whitney test
#' \item pValueAdj: Adjusted p-value according to \code{correction} parameter
#' \item expVar: Explained variance by the principal component
#' \item cumExpVar: All principal components are ranked accoriding to their frequency of ocurrence and their variance explained. 
#' This column contains the cumulative variance explained across the ranked principal components
#' }
#' \item \code{pVar}: Column name from metadata to use as the variable to predict using
#' the informative principal components. Informative principal components are selected based on this variable.
#' }
#' @keywords informative, significant, features
#' @importFrom methods is
#' @export
#' @author
#' José Alquicira Hernández




getInformativePCs <- function(object, pVar, varLim = 0.01, correction = "fdr", sig = 0.1){
  if(!is(object, "eigenPred")){
    stop("Invalid class for object: must be 'eigenPred'")
  }
  
  if(!any(correction %in% stats::p.adjust.methods)){
    stop("Invalid multiple testing correction method. See ?p.adjust function")
  }
  
  i <- object@expVar > varLim
  
  pca <- getPCA(object)[,i]
  classes <- metadata(object)[[pVar]]
  
  testDiff <- function(pc){
    wilcox.test(pc[classes == levels(classes)[1]], pc[classes == levels(classes)[2]])
  }
  
  lapply(pca,testDiff) %>% 
    lapply('[[', "p.value") %>% 
    as.data.frame() %>% 
    gather(key = "PC", value = "pValue") %>% 
    mutate(pValueAdj = p.adjust(pValue, method = correction, n = nrow(.))) %>% 
    arrange(pValueAdj) %>% 
    filter(pValueAdj < sig) %>% 
    mutate(expVar = object@expVar[match(PC, names(object@expVar))]) %>% 
    mutate(PC = factor(PC, levels = PC), cumExpVar = cumsum(expVar)) -> sigPCs
  
  
  
  object@features <- sigPCs
  object@pVar <- pVar
  
  message("DONE!")
  object
  
}


