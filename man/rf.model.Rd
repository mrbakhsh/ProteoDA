% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rf.model.R
\name{rf.model}
\alias{rf.model}
\title{rf.model}
\usage{
rf.model(
  features,
  gd,
  cv_fold = 5,
  plots = FALSE,
  filename = file.path(tempdir(), "plotsRF.pdf")
)
}
\arguments{
\item{features}{A data frame with bait-prey or pre-prey associations in the
first column, and features to be passed to the classifier in the remaining
columns.}

\item{gd}{A gold reference set including true associations with class labels
indicating if such PPIs are positive or negative.}

\item{cv_fold}{Number of partitions for cross-validation; defaults to 5.}

\item{plots}{Logical value, indicating whether to plot the performance of the
learning algorithm using k-fold cross-validation; defaults to FALSE.If the
argument set to TRUE, plots will be saved in the temp() directory. These
plots are : \itemize{ \item{pr_plot} - Precision-recall PLOT
\item{roc_plot} - ROC plot \item{radar_plot} - Radar plot showing accuracy,
F1-score , positive predictive value (PPV), sensitivity (SE) and MCC. }}

\item{filename}{character string, indicating the output filename as an pdf
object. Defaults to plotsRF.pdf.}
}
\value{
predicted_interactions- Prediction scores for input dataset (i.e.,
features) cm_ouput \itemize{ \item{cm} - A confusion matrix. \item{ACC} -
Accuracy. \item{SE} - Sensitivity. \item{SP} - Specificity. \item{PPV} -
Positive Predictive Value. \item{F1} - F1-score. \item{MCC} - Matthews
correlation coefficient } roc_auc - area under ROC curve pr_auc - area under
PR curve
}
\description{
This function uses Random Forest classifier to assigns each
  bait-prey pair a confidence score, indicating the level of support for that
  pair of proteins to interact. This function also computes and plots
  Precision-Recall (PR), ROC curve, and randar plot including Sensitivity
  (SE), Specificity (SP), Accuracy (ACC), Positive Predictive Value (PPV),
  F1-acore, and Matthews correlation coefficient (MCC) to evaluate the
  performance of the classifier on k-fold cross-validated data.
}
\examples{
data("Testdata_1")
# Score using spoke model
scdata <- score_APMS(Testdata_1, cPASS_m = TRUE, MiST_m = TRUE)
# create random train data
library(dplyr)
t <- sample_n(scdata, 500)
t$label <- # random lables
sample(c("Positive","Negative"), size = nrow(t), replace = TRUE)
t <- t[, c(1,14)] #select pair with the class labels
# perfrorm prediction
pred <- rf.model(scdata, t, cv_fold = 5,  plots = FALSE)
}
\author{
Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
}
