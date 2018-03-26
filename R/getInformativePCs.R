#' @title Get informative principal components
#' @description Given a prediction variable, finds a set of class-informative principal components that are significant across n subsets of cells from the 
#' PCA matrix. A Wilcoxon rank sum test is used to determine a difference between the value distributions of classes from the prediction variable.
#' @param object An \code{eigenPred} object
#' @param pVar Prediction variable corresponding to a column in \code{metadata} slot
#' @param rep Number of random samples to assess stable principal components
#' @param seed Seed to create random samples
#' @param varLim Threshold to filter principal components based on variance explained.
#' @param correction Multiple testing correction method. Default: false discovery rate. See \code{p.adjust} function 
#' @param sig Significance level to determine principal components explaining class identity
#' @return A \code{eigenPred} object with three additional filled slots:
#' \itemize{
#' \item \code{features}: A data frame with significant principal components the following information:
#' \itemize{
#' \item PC: Principal component
#' \item Freq: Frequency of occurence of the principal component over a number of random samples from the PCA matrix
#' \item expVar: Explained variance by the principal component
#' \item cumExpVar: All principal components are ranked accoriding to their frequency of ocurrence and their variance explained. 
#' This column contains the cumulative variance explained across the ranked principal components
#' }
#' \item \code{pVar}: Column name from metadata to use as the variable to predict using
#' the informative principal components. Informative principal components are selected based on this variable.
#' \item \code{rep}: Number of random samples to determine stable principal components
#' }
#' @keywords informative, significant, features
#' @importFrom methods is
#' @export
#' @author
#' José Alquicira Hernández




getInformativePCs <- function(object, pVar, rep = 100, seed = NULL, varLim = 0.01, correction = "fdr", sig = 0.05){
  if(!is(object, "eigenPred")){
    stop("Invalid class for object: must be 'eigenPred'")
  }
  
  if(!any(correction %in% stats::p.adjust.methods)){
    stop("Invalid multiple testing correction method. See ?p.adjust function")
  }
  
  if(!is.null(seed)){
    set.seed(seed)   
    seeds <- sample(1:1000000L, rep)
  }else{
    seeds <- sample(1:1000000L, rep)
  }
  
  res <- list()
  
  message("Identifying stable features...")
  pb <- txtProgressBar(min = 0, max = length(seeds), style = 3)
  
  for(i in seq_len(length(seeds))){
    res[[i]] <- .getSignificantPCs(object,
                                   pVar = pVar,
                                   seed = seeds[i],
                                   varLim = varLim,
                                   correction = correction,
                                   sig = sig)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  message("Gathering results...")
  
  res %>% 
    lapply("[", "PC") %>% 
    unlist() %>% 
    table() %>% 
    as.data.frame() %>% 
    set_names(c("PC", "Freq")) %>%
    mutate(expVar = object@expVar[match(PC, names(object@expVar))]) %>% 
    arrange(-Freq, -expVar) %>% 
    mutate(PC = factor(PC, levels = PC), cumExpVar = cumsum(expVar)) -> sigPCs
  
  
  object@features <- sigPCs
  object@pVar <- pVar
  object@rep <- rep
  
  message("DONE!")
  object
  
}



.getSignificantPCs <- function(object, pVar, seed, varLim = 0.01, correction = "fdr", sig = 0.05){
  nCells <- nrow(object@metadata)
  set.seed(seed)
  i <- createDataPartition(object@metadata[[pVar]], p = 0.5, list = FALSE, times = 1)
  
  # Filter PCs by variance explained
  pcFilter <- object@expVar > varLim
  pcs <- names(pcFilter)[pcFilter]
  
  
  subpca <- getPCA(object)[i, pcs]
  submeta <- metadata(object)[i, pVar]
  
  lapply(subpca, function(pc){
    wilcox.test(pc[submeta == levels(submeta)[1]],
                pc[submeta == levels(submeta)[2]])}) %>% 
    lapply('[[', "p.value") %>% 
    as.data.frame() %>% 
    gather(key = "PC", value = "pValue") %>% 
    mutate(pValueAdj = p.adjust(pValue, method = correction, n = nrow(.))) %>% 
    arrange(pValueAdj) %>% 
    filter(pValueAdj < sig) -> pcaSig
  
  pcaSig
  
}
