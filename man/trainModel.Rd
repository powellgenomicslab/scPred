% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trainModel.R
\name{trainModel}
\alias{trainModel}
\title{Train a prediction model}
\usage{
trainModel(object, model = "svmRadial", resampleMethod = "cv",
  seed = NULL, number = 10, returnData = FALSE, savePredictions = FALSE)
}
\arguments{
\item{object}{An \code{scPred} object with informative PCs obtained using 
the \code{getInformativePCs} function}

\item{model}{Classification model supported via \code{caret} package
Default: support vector machine with polynomial kernel}

\item{resampleMethod}{Resample model used in \code{trainControl} function. Default: K-fold cross validation}

\item{seed}{Seed to apply the resample model}

\item{number}{Number of iterations for resample model. See \code{trainControl} function}

\item{returnData}{If \code{TRUE}, training data is returned}

\item{savePredictions}{If \code{TRUE}, predictions for training data are returned}
}
\value{
A \code{train} object with final results. See \code{train} function for details. An aditional value \code{top} is added to the 
\code{train} object to track the features used to train the model
}
\description{
Trains a prediction model from an \code{scPred} object
}
\author{
José Alquicira Hernández
}
\keyword{model}
\keyword{train,}