#' @title Align low-dimensional space using Seurat algorithm
#' @description Uses the manifold-alignment Seurat algorithm to align the training eigenspace and the prediction projection
#' See ?AlignSubspace() for more details. Note: this helper function is a modified version from Seurat.
#' @author José Alquicira Hernández
#' @importFrom  pbapply pbsapply


.alignSubspaceSeurat <- function (object, reduction.type = "pca.scpred", grouping.var = "dataset_dummy", dims.align, 
          num.possible.genes = 2000, num.genes = 30, show.plots = FALSE, 
          verbose = TRUE, ...) 
{
  parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("AlignSubspace"))]
  object <- Seurat:::SetCalcParams(object = object, calculation = paste0("AlignSubspace.", 
                                                                reduction.type), ... = parameters.to.store)
  ident.orig <- object@ident
  object <- SetAllIdent(object = object, id = grouping.var)
  levels.split <- names(x = sort(x = table(object@ident), decreasing = T))
  num.groups <- length(levels.split)
  objects <- list()
  for (i in 1:num.groups) {
    objects[[i]] <- SubsetData(object = object, ident.use = levels.split[i])
  }
  object@ident <- ident.orig
  cc.loadings <- list()
  scaled.data <- list()
  cc.embeds <- list()
  for (i in 1:num.groups) {
    if (verbose) {
      cat(paste0("Rescaling group ", i, "\n"), file = stderr())
    }
    objects[[i]] <- ScaleData(object = objects[[i]], display.progress = verbose, 
                              ...)
    objects[[i]]@scale.data[is.na(x = objects[[i]]@scale.data)] <- 0
    objects[[i]] <- ProjectDim(object = objects[[i]], reduction.type = reduction.type, 
                               do.print = FALSE)
    cc.loadings[[i]] <- GetGeneLoadings(object = objects[[i]], 
                                        reduction.type = reduction.type, use.full = TRUE)
    cc.embeds[[i]] <- GetCellEmbeddings(object = objects[[i]], 
                                        reduction.type = reduction.type)
    scaled.data[[i]] <- objects[[i]]@scale.data
  }
  cc.embeds.all <- GetCellEmbeddings(object = object, reduction.type = reduction.type, 
                                     dims.use = dims.align)
  colnames(cc.embeds.all) <- paste0("A", colnames(x = cc.embeds.all))
  cc.embeds.orig <- cc.embeds.all
  for (cc.use in dims.align) {
    for (g in 2:num.groups) {
      if (verbose) {
        cat(paste0("Aligning dimension ", cc.use, "\n"), 
            file = stderr())
      }
      genes.rank <- data.frame(rank(x = abs(x = cc.loadings[[1]][, 
                                                                 cc.use])), rank(x = abs(x = cc.loadings[[g]][, 
                                                                                                              cc.use])), cc.loadings[[1]][, cc.use], cc.loadings[[g]][, 
                                                                                                                                                                      cc.use])
      genes.rank$min <- apply(X = genes.rank[, 1:2], MARGIN = 1, 
                              FUN = min)
      genes.rank <- genes.rank[order(genes.rank$min, decreasing = TRUE), 
                               ]
      genes.top <- rownames(x = genes.rank)[1:min(num.possible.genes, 
                                                  nrow(genes.rank))]
      bicors <- list()
      for (i in c(1, g)) {
        cc.vals <- cc.embeds[[i]][, cc.use]
        if (verbose) {
          bicors[[i]] <- pbsapply(X = genes.top, FUN = function(x) {
            return(Seurat:::BiweightMidcor(x = cc.vals, y = scaled.data[[i]][x, 
                                                                    ]))
          })
        }
        else {
          bicors[[i]] <- sapply(X = genes.top, FUN = function(x) {
            return(Seurat:::BiweightMidcor(x = cc.vals, y = scaled.data[[i]][x, 
                                                                    ]))
          })
        }
      }
      genes.rank <- data.frame(rank(x = abs(x = bicors[[1]])), 
                               rank(x = abs(x = bicors[[g]])), bicors[[1]], 
                               bicors[[g]])
      genes.rank$min <- apply(X = abs(x = genes.rank[, 
                                                     1:2]), MARGIN = 1, FUN = min)
      genes.rank <- genes.rank[sign(genes.rank[, 3]) == 
                                 sign(genes.rank[, 4]), ]
      genes.rank <- genes.rank[order(genes.rank$min, decreasing = TRUE), 
                               ]
      genes.use <- rownames(x = genes.rank)[1:min(num.genes, 
                                                  nrow(genes.rank))]
      if (length(genes.use) == 0) {
        stop("Can't align group ", g, " for dimension ", 
             cc.use)
      }
      metagenes <- list()
      multvar.data <- list()
      for (i in c(1, g)) {
        scaled.use <- sweep(x = scaled.data[[i]][genes.use, 
                                                 ], MARGIN = 1, STATS = sign(x = genes.rank[genes.use, 
                                                                                            which(c(1, g) == i) + 2]), FUN = "*")
        scaled.use <- scaled.use[, names(x = sort(x = cc.embeds[[i]][, 
                                                                     cc.use]))]
        metagenes[[i]] <- (cc.loadings[[i]][genes.use, 
                                            cc.use] %*% scaled.data[[i]][genes.use, ])[1, 
                                                                                       colnames(x = scaled.use)]
      }
      mean.difference <- mean(x = Seurat:::ReferenceRange(x = metagenes[[g]])) - 
        mean(x = Seurat:::ReferenceRange(x = metagenes[[1]]))
      align.1 <- Seurat:::ReferenceRange(x = metagenes[[g]])
      align.2 <- Seurat:::ReferenceRange(x = metagenes[[1]])
      a1q <- sapply(X = seq(from = 0, to = 1, by = 0.001), 
                    FUN = function(x) {
                      return(quantile(x = align.1, probs = x))
                    })
      a2q <- sapply(X = seq(from = 0, to = 1, by = 0.001), 
                    FUN = function(x) {
                      quantile(x = align.2, probs = x)
                    })
      iqr <- (a1q - a2q)[100:900]
      iqr.x <- which.min(x = abs(x = iqr))
      iqrmin <- iqr[iqr.x]
      if (show.plots) {
        print(iqrmin)
      }
      align.2 <- align.2 + iqrmin
      alignment <- dtw::dtw(x = align.1, y = align.2, keep.internals = TRUE)
      alignment.map <- data.frame(alignment$index1, alignment$index2)
      alignment.map$cc_data1 <- sort(cc.embeds[[g]][, cc.use])[alignment$index1]
      alignment.map$cc_data2 <- sort(cc.embeds[[1]][, cc.use])[alignment$index2]
      alignment.map.orig <- alignment.map
      alignment.map$dups <- duplicated(x = alignment.map$alignment.index1) | 
        duplicated(x = alignment.map$alignment.index1, 
                   fromLast = TRUE)
      alignment.map <- alignment.map %>% group_by(alignment.index1) %>% 
        mutate(cc_data1_mapped = ifelse(dups, mean(cc_data2), 
                                        cc_data2))
      alignment.map <- alignment.map[!duplicated(x = alignment.map$alignment.index1), 
                                     ]
      cc.embeds.all[names(x = sort(x = cc.embeds[[g]][, 
                                                      cc.use])), cc.use] <- alignment.map$cc_data1_mapped
      if (show.plots) {
        par(mfrow = c(3, 2))
        plot(x = Seurat:::ReferenceRange(x = metagenes[[1]]), 
             main = cc.use)
        plot(x = Seurat:::ReferenceRange(x = metagenes[[g]]))
        plot(x = Seurat:::ReferenceRange(x = metagenes[[1]])[(alignment.map.orig$alignment.index2)], 
             pch = 16)
        points(x = Seurat:::ReferenceRange(metagenes[[g]])[(alignment.map.orig$alignment.index1)], 
               col = "red", pch = 16, cex = 0.4)
        plot(x = density(x = alignment.map$cc_data1_mapped))
        lines(x = density(x = sort(x = cc.embeds[[1]][, 
                                                      cc.use])), col = "red")
        plot(x = alignment.map.orig$cc_data1)
        points(x = alignment.map.orig$cc_data2, col = "red")
      }
    }
  }
  new.type <- paste0(reduction.type, ".aligned")
  new.key <- paste0("A", GetDimReduction(object = object, reduction.type = reduction.type, 
                                         slot = "key"))
  object <- Seurat:::SetDimReduction(object = object, reduction.type = new.type, 
                            slot = "cell.embeddings", new.data = cc.embeds.all)
  object <- Seurat:::SetDimReduction(object = object, reduction.type = new.type, 
                            slot = "key", new.data = new.key)
  return(object)
}