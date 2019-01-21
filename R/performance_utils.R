#' @title Prepare and tidy prediction results
#' @description Formats prediction results with true values (if known) contained in metadata. Optionally,
#'  adds an extra column for new prediction labels
#' @param object \code{scPred} object
#' @param true Column name in \code{predMeta} slot that corresponds to the true known classes 
#' @param pred Column name in \code{predMeta} slot that corresponds to the predicted classes
#' if they have been assigned independently from the \code{scPredict()} function
#' @return A list with 
#' \itemize{
#' \item \code{predictions}. A dataframe with two extra columns: true and prediction columns from the prediction metadata
#' and/or from class labels assigned by \code{scPred}
#' \item \code{true}. Renamed true name column
#' \item \code{pred}. Renamed pred name column
#' \item 
#' }
#' @author José Alquicira Hernández
#' 

processPreds <- function(object, true, pred = NULL){
  
  # Get predictions
  predictions <- getPredictions(object)
  
  # Check if there is metadata where the true and/or the predicted classes are
  # stored
  
  if(!length(object@predMeta)){
    stop("No metadata associated with prediction data has been assigned")
  }
  
  # Check if true column name is in metadata
  if(!true %in% names(object@predMeta)){
    stop(true, " variable is not present in prediction metadata")
  }
  
  
  # Check if pred column name was provided. If so, assign it to the predictions 
  # dataframe.
  
  if(is.null(pred)){
    pred <- "predClass"
  }else if(!pred %in% names(object@predMeta)){
    stop(pred, " prediction variable is not present in prediction metadata")
  }else{
    predictions[[pred]] <- object@predMeta[[pred]]
  }
  
  
  if(pred == true){
    trueLab <- paste0(true, ".")
    predictions[[trueLab]] <- object@predMeta[[true]]
    true <- trueLab
  }else{
    predictions[[true]] <- object@predMeta[[true]]
  }
  
  list(predictions = predictions, true = true, pred = pred)
}