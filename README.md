[![DOI](https://zenodo.org/badge/206724827.svg)](https://zenodo.org/badge/latestdoi/206724827)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://zenodo.org/badge/latestdoi/206724827)

# scPred: accurate supervised method for cell-type classification from single-cell RNA-seq data


`scPred` is a general method to predict cell types based on variance structure decomposition.
It selects the most cell type-informative principal components from a dataset and trains a prediction model for each cell type. The principal training axes are projected onto the test dataset to obtain the PCs scores for the test dataset and the trained model(s) is/are used to classify single cells.

For more details see our paper in **Genome Biology**:

[scPred: accurate supervised method for cell-type classification from single-cell RNA-seq data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1862-5)

This [introduction to scPred](https://joseah.github.io/post/introduction-to-scpred/) shows a basic workflow for cell type prediction.

You can install scPred via devtools as follows:

```r
devtools::install_github("powellgenomicslab/scPred")
```
