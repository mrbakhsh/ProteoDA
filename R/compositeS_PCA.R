  #'compositeS_PCA
  #' @param Sfile Data frame containing numerical features.
  #' @param scale	A logical value indicating whether the variables should be
  #' scaled to have unit variance before the analysis takes place. See also
  #' (\code{\link{prcomp}}).
  #' @return Data frame containing single composite score (last column)
  #' generated from principal component analysis (PCA).
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @importFrom stats prcomp
  #' @description This function performs a Principal Component Analysis (PCA)
  #' on the numerical feature matrix and select weights that project the
  #' feature values on the first principal component. We recommend using this
  #' option only when insufficient training data is available.
  #' @export
  #' @examples
  #' data("Testdata_1")
  #' # drop the length column
  #' Testdata_1 <- Testdata_1[, -6]
  #' # score file containing spectral counts
  #' scoredfile1 <- cPASS_m(Testdata_1)
  #' scoredfile1 <- scoredfile1[, -1]
  #' # apply PCA
  #' compS <- compositeS_PCA(scoredfile1,  scale = FALSE)


  compositeS_PCA <-
    function(Sfile , scale = TRUE) {


      if(any(is.na(Sfile)) == TRUE) {
        stop("Data includes NA values")
      }


      if(any(is.character(Sfile)) == TRUE) {
        stop("Data includes character vectors")
      }

      pca <-
        prcomp(Sfile[, 2:ncol(Sfile)], scale. = scale)
      # select the first component
      Sfile$CompositeScore <- pca$x[,1]

      # normalize the data
      Sfile$CompositeScore <-
        (Sfile$CompositeScore - min(Sfile$CompositeScore))/
        (max(Sfile$CompositeScore)-min(Sfile$CompositeScore))


      return(Sfile)

    }


