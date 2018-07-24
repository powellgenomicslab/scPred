## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load_libraries,  message=FALSE, results='hide', warning=FALSE-------
library("scPred")
library("tidyverse")

## ----download_data, eval=FALSE-------------------------------------------
#  download.file("https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/pollen.rds", destfile = "~/Downloads/pollen.rds")

## ----read_data, message=FALSE, results='hide', warning=FALSE-------------
require(SingleCellExperiment)
pollen <- readRDS("~/Downloads/pollen.rds")

## ----get_data------------------------------------------------------------
pollen_counts <- normcounts(pollen)
pollen_cpm  <- apply(pollen_counts, 2, function(x) (x/sum(x))*1000000)
pollen_metadata <- as.data.frame(colData(pollen))

## ----get_counts_per_class------------------------------------------------
table(pollen_metadata$cell_type2)

## ----create_partitions---------------------------------------------------
set.seed(1234)
i <- createDataPartition(pollen_metadata$cell_type2, p = 0.70, list = FALSE)
train_data <- t(pollen_cpm[, i])
test_data <- t(pollen_cpm[, -i])

train_info <- pollen_metadata[i, , drop = FALSE]
test_info <- pollen_metadata[-i, , drop = FALSE]

## ----eigendecomposition, results=FALSE, message=FALSE, warning=FALSE-----
set.seed(1234)
scp <- eigenDecompose(train_data, n = 10)

## ----metadata_assignment-------------------------------------------------
scPred::metadata(scp) <- train_info

## ----get_feature_spaces--------------------------------------------------
scp <- getFeatureSpace(scp, pVar = "cell_type2")

## ----show_features-------------------------------------------------------
scp@features

## ----plot_pca, fig.height=5, fig.width=6---------------------------------
plotEigen(scp, group = "cell_type2")

## ----train_models--------------------------------------------------------
scp <- trainModel(scp, seed = 66)

## ----show_summary--------------------------------------------------------
scp

## ----get_train_results, fig.height=5, fig.width=6------------------------
res <- getTrainResults(scp)

## ----plot_train_probabilities,  fig.height=5, fig.width=7----------------
mapply(function(x,y){dplyr::rename(x, probability = !! enquo(y))}, res, names(res), SIMPLIFY = FALSE) %>% 
  Reduce(rbind, .) -> train_probabilities

train_probabilities %>% 
    select(object, obs, probability, "other") %>% 
    ggplot() +
    aes(probability, fill = obs) +
    geom_histogram(bins = 30, color = "black") +
    geom_vline(xintercept = 0.9, color = "red") +
  facet_wrap(~object) +
    theme_bw()


## ----prediction----------------------------------------------------------
predictions <- scPredict(scp, newData = test_data, threshold = 0.9)

## ----print_predictions---------------------------------------------------
predictions

## ----plot_prediction_porbs, fig.height=5, fig.width=7--------------------
predictions %>% 
  mutate(true = test_info$cell_type2) %>% 
  gather(key = "model", value = "probability", 1:4) %>% 
  ggplot() +
  aes(probability, fill = true) +
  geom_histogram(color = "black") +
  geom_vline(xintercept = 0.9, color = "red") +
  facet_wrap(~model) +
  theme_bw() +
  scale_fill_viridis_d()

## ----session_info---------------------------------------------------
options(width = 70)
devtools::session_info(include_base = TRUE)

