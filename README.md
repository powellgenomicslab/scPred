# scPred: Single cell prediction using singular value decomposition and machine learning classification


`scPred` is a general method to predict cell types based on variance structure decomposition.
It selects the most cell type-informative principal components from a dataset and trains a prediction model for each cell type. The principal training axes are projected onto the test dataset to obtain the PCs scores for the test dataset and the trained model(s) is/are used to classify single cells.

For more details see our pre-print on **bioRxiv**:

[scPred: Single cell prediction using singular value decomposition and machine learning classification](https://www.biorxiv.org/content/early/2018/07/15/369538)

This [introduction to scPred](https://joseah.github.io/scPred_docs/) shows a basic workflow for cell type prediction.

# Authors

- [José Alquicira-Hernández](https://graduate-school.uq.edu.au/profile/530/jose)
- [Quan Nguyen](https://imb.uq.edu.au/profile/1672/quan-nguyen)
- [Joseph E Powell](https://qbi.uq.edu.au/profile/469/joseph-powell)
