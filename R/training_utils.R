#' @title Get training predictions
#' @description For all trained models, retrieves the predictions and associated probabilities
#' @param object An \code{scPred} or \code{seurat} object
#' @return An list with the prediction results from each model
#' @export
#' @author José Alquicira Hernández

getTrainPred <- function(object){
  
  if(!is(object, "scPred")){
    stop("'object' must be of class 'scPred'")
  }
  
  if(!length(object@train)){
    stop("No training model found!")
  }
  
  
  preds <- lapply(object@train, "[[", "pred")
  ids <- lapply(object@train, "[[", "trainingData") %>% 
    lapply(rownames)
  
  mapply(function(x, y){ rownames(x) <- y; x}, preds, ids, SIMPLIFY = FALSE)
  
}


#' @title Plot training probabilities for all models
#' @description Plots the probability distributions for each model.
#' @param object An \code{scPred} or \code{seurat} object
#' @param ... Extra arguments for \code{plot_grid} function
#' @return A plot with the probability distributions. Each panel corresponds to a trained model (cell class)
#' and the positive class to the class of interest
#' @importFrom cowplot plot_grid get_legend
#' @export
#' @author José Alquicira Hernández

plotTrainProbs <- function(object, ...){
  
  preds <- getTrainPred(object)
  namesClasses <- names(preds)
  
  plotProbs <- function(data, class){
    
    ggplot(data, aes_string(class, fill = "obs")) +
      geom_histogram(color = "black", size = 0.1) +
      xlab(paste0("P(", class,")")) +
      ylab("Number of cells") +
      labs(fill = "Prediction") +
      scale_fill_manual(values = c("#1874CD", "#EE2C2C")) +
      theme_classic()
    
  }
  
  plots <- mapply(plotProbs, preds, namesClasses, SIMPLIFY = FALSE) -> plot_res
  n <- length(plots)
  
  if(n > 1){
    legend <- get_legend(plots[[1]] +
                           scale_fill_manual(labels = c("+ class", "- class"),
                                             values = c("#1874CD", "#EE2C2C"))
    )
    
    plots <- lapply(plots, function(p) p + theme(legend.position = "none"))
    plots[2:n] <- lapply(plots[2:n], function(p) p + ylab(""))
    plots[2:n] <- lapply(plots[2:n], function(p) p + ylab(""))
    plots <- plot_grid(plotlist = plots, ...)
    plots <- plot_grid(plots, legend, rel_widths = c(n, 1/n), rel_heights = c(n, 1/n), ...)
    
    
  }
  
  plots
  
}