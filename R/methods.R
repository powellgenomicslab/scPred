#' @title show
#' @description Generic display function for \linkS4class{eigenPred} objects. Displays summary of the
#' object such as number of cells, genes, significant features.
#' @importFrom methods setMethod
#' @export

setMethod("show", signature("eigenPred"), function(object) {
  
  cat("'eigenPred' object\n")
  
  nCells <- nrow(object@prcomp$x)
  nPCs <- ncol(object@prcomp$x)
  
  cat(sprintf("  Number of cells: %i\n", nCells))
  cat(sprintf("  Number of genes: %i\n", nrow(getLoadings(object))))
  cat(sprintf("  Number of principal components: %i\n", nPCs))
  cat("  Maximum variance explained by a single PC:", max(object@expVar), "\n")
  
  if(nrow(object@features) != 0){
    cat("  ---\n")
    cat(sprintf("  Number of significant features: %i\n", nrow(object@features)))
  }
  
  if(length(object@pVar) != 0 & nrow(object@metadata) > 0){
    cat(sprintf("  Prediction variable: %s\n", object@pVar))
    metadata(object) %>% 
      select(one_of(object@pVar)) %>% 
      table() -> freq 
      cat("                       |_", names(freq)[1], ":", freq[1], "\n")
      cat("                       |_", names(freq)[2], ":", freq[2], "\n")

  }
  
})


#' @title Get metadata
#' @description Gets metadata from \code{eigenPred} object
#' @importFrom methods setMethod
#' @export


setGeneric("metadata", function(object) standardGeneric("metadata"))
setMethod("metadata", "eigenPred", function(object) {
  out <- object@metadata
  return(out)
})


#' @title Set metadata
#' @description Sets metadata to a \code{eigenPred} object. Metadata must be a dataframe
#' \itemize{
#' \item row names: ids matching the column names of the expression matrix
#' \item columns: associated metadata such as cell type, conditions, sample or, batch. 
#' } 
#' @importFrom methods setMethod
#' @export

setGeneric("metadata<-", function(object, value) standardGeneric("metadata<-"))
setMethod("metadata<-", signature("eigenPred"), 
          function(object, value) {
            if(!all(rownames(getPCA(object)) == row.names(value))){
              stop("Cells Ids do not match cell IDs in metadata")
            }
            object@metadata <- value
            object
          }
)



#' @title Get principal components
#' @description Gets matrix with principal components from a \code{eigenPred} object
#' @importFrom methods setMethod
#' @export

setGeneric("getPCA", def = function(object) {
  standardGeneric("getPCA")
})


#' @title Get principal components
#' @description Gets matrix with principal components from a \code{eigenPred} object
#' @importFrom methods setMethod
#' @export

setMethod("getPCA", signature("eigenPred"), function(object) {
  return(as.data.frame(object@prcomp$x))
})


#' @title Get loadings matrix
#' @description Gets rotation matrix (right singular vectors) from a \code{eigenPred} object
#' @importFrom methods setMethod
#' @export

setGeneric("getLoadings", def = function(object) {
  standardGeneric("getLoadings")
})

#' @title Get loadings matrix
#' @description Gets rotation matrix (right singular vectors) from a \code{eigenPred} object
#' @importFrom methods setMethod
#' @export

setMethod("getLoadings", signature("eigenPred"), function(object) {
  return(object@prcomp$rotation)
})


#' @title Plot PCA
#' @description Plots PCA with raw results or group by a variable
#' @importFrom methods setMethod
#' @export

setGeneric("plotEigen", def = function(object, group = NULL, pc = c(1,2), geom = c("points", "density_2d", "both")) {
  standardGeneric("plotEigen")
})

#' @title Plot PCA
#' @description Plots PCA 
#' @importFrom methods setMethod
#' @export


setMethod("plotEigen", signature("eigenPred"), function(object, group = NULL, pc = c(1,2), geom = c("both", "points", "density_2d")){
  
  geom <- match.arg(geom)
  pca <- getPCA(object)[,pc] 
  
  
  if(!is.null(group)){
    
    if(ncol(object@metadata) == 0){
      stop("No metadata has been assigned to 'eigenPred' object")
    }
    
    if(!any(group %in% names(object@metadata))){
      stop("'group' variable is not included in metadata")
    }
    
    if(!is(object@metadata[[group]], "factor")){
      stop("'group' variable must be a factor variable")
    }
    
    metadata <- object@metadata[group]
    pca <- cbind(pca, metadata)
    p <- ggplot(pca, aes_string(x = names(pca)[1], y = names(pca)[2], color = group))
    
  }else{
    p <- ggplot(pca, aes_string(x = names(pca)[1], y = names(pca)[2]))
  }
  
  if(geom == "points" | geom == "both"){
    p <- p + geom_point()
  }
  if(geom == "density_2d" | geom == "both"){
    p <- p +  geom_density_2d() 
  }
  p + 
    scale_color_brewer(palette = "Set1") +
    theme_bw()
})



#' @title Generate diagnostic plot
#' @description Plots a graph showing the frequency of each principal component across a number of random samples and the cumulative variance
#' explained by them. This plot may be used to select stable principal components to be used to train a prediction model.
#' @importFrom methods setMethod
#' @export

setGeneric("plotFeatures", def = function(object) {
  standardGeneric("plotFeatures")
})

#' @title Generate diagnostic plot
#' @description Plots a graph showing the frequency of each principal component across a number of random samples and the cumulative variance
#' explained by them. This plot may be used to select stable principal components to be used to train a prediction model.
#' @importFrom methods setMethod
#' @export


setMethod("plotFeatures", signature("eigenPred"), function(object){
  object@features %>% 
    mutate(cumExpVarScaled = ( (cumExpVar - min(cumExpVar)) / (max(cumExpVar) - min(cumExpVar)) ) * (max(Freq) - min(Freq)) + min(Freq)) -> pcFreq
  
  ggplot(pcFreq, aes(PC, Freq)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_line(aes(as.numeric(PC), cumExpVarScaled)) +
    geom_point(aes(as.numeric(PC), cumExpVarScaled), shape = 1) +
    scale_y_continuous(sec.axis = sec_axis(~( (. - min(pcFreq$Freq)) / (max(pcFreq$Freq) - min(pcFreq$Freq)) ) * (max(pcFreq$cumExpVar) - min(pcFreq$cumExpVar)) + min(pcFreq$cumExpVar), 
                                           name = "Cumulative explained variance")) +
    theme_bw() +
    xlab("Significant principal components") +
    ylab(paste("Frequency of occurence in", object@rep, "interations")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) -> p
  
  p
  
  
})
