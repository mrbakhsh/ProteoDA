
  .merge.all <- function(x, y) {
    merge(x, y, all=TRUE, by="b_p")
  }


  #' score_APMS
  #' @title Scoring Affinity Purification Mass Spectrometry (AP-MS) Datasets
  #' via spoke and matrix models.
  #' @param datInput Data frame with column names:
  #' Experiment.id, Replicate, Bait, Prey, counts (i.e., summed spectral
  #' counts or log2-transformed summed MS1-intensities), and
  #' prey protein length.
  #' @param cPASS_m If TRUE, scores the associations
  #' using Comparative Proteomic Analysis Software Suite (CompPASS) model.
  #' This algorithm is based on spoke model.
  #' @param MiST_m If TRUE, scores the associations using
  #' Mass spectrometry interaction STatistics (MiST) model.
  #' This algorithm is based on spoke model.
  #' @param HG_m If TRUE, scores the associations
  #' via hypergeometric probability with
  #' incorporation of NSAF. This algorithm is based on matrix model.
  #' @param PCC_m If TRUE, scores the associations
  #' using PCC score. This algorithm is based on matrix model.
  #' @param WM_m If TRUE, scores the associations
  #' using the weighted matrix model.This algorithm is based on matrix model.
  #' @param Dice_m If TRUE, scores the associations
  #' using Dice metric.This algorithm is based on matrix model.
  #' @param Jaccard_m If TRUE, scores the associations
  #' using Jaccard metric.This algorithm is based on matrix model.
  #' @param Overlap_m If TRUE, scores the associations
  #' using Overlap score.This algorithm is based on matrix model.
  #' @param Zscore_m If TRUE, scores the associations
  #' using Z-score.This algorithm is based on matrix model.
  #' @param iter Number of bootstraps to execute to estimate z scores.
  #' Defaults to 3.
  #' @return A data frame containing the calculated features for Bait-Prey
  #' or Prey-Prey associations.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @description This function computes Bait-prey and Prey-Prey associations
  #' via nine metrics: (Comparative Proteomic Analysis Software Suite (CompPASS)
  #' (\code{\link{cPASS_m}}),
  #' Mass spectrometry interaction STatistics (\code{\link{MiST_m}}),
  #' Hypergeometric probability (\code{\link{HG_m}}),
  #' Weighted matrix model (\code{\link{WM_m}}),
  #' Dice (\code{\link{Dice_m}}), Jaccard (\code{\link{Jaccard_m}}),
  #' ovelrap (\code{\link{Overlap_m}}),Z_score
  #' (\code{\link{Zscore_m}}), and Pearson correlation coefficients
  #' (\code{\link{PCC_m}}).By default ComPASS and MiST scoring is enabled,
  #' optionally, one or more of these can be enabled or disabled.
  #' @export
  #' @examples
  #' # score file containing spectral counts
  #' data(Testdata_1)
  #' scoredfile1 <-  score_APMS(Testdata_1)
  #'
  #' # score file containing MS1-intensities
  #' # 1. first we log transform the counts
  #' library(dplyr)
  #' data("Testdata_2")
  #' Testdata_2 <-
  #'   Testdata_2 %>%
  #'   group_by(Bait, Replicate) %>%
  #'   mutate(counts = log2(counts))
  #' # 2. then we apply score_APMS
  #' scoredfile2 <- score_APMS(Testdata_2)


  score_APMS <-
    function(datInput,
             cPASS_m = TRUE,
             MiST_m = TRUE,
             HG_m = FALSE,
             WM_m = FALSE,
             PCC_m = FALSE,
             Dice_m = FALSE,
             Jaccard_m = FALSE,
             Overlap_m = FALSE,
             Zscore_m = FALSE,
             iter = 3){



      clnames <-
        c("Experiment.id","Replicate", "Bait", "Prey", "counts", "Length")
      if(!is.data.frame(datInput)) datInput <- as.data.frame(datInput)

      if(!all(clnames %in% colnames(datInput))) {
        missCols <-
          setdiff(clnames,colnames(datInput))
        stop("columns are missing: ", paste(missCols, collapse = ", "))
      }

      if(is.numeric(datInput$counts) ==  FALSE){
        stop("counts must contain numerical values")
      }

      if(ncol(datInput) != 6){
        stop("Data frame must contain 6 columns including:
          Experiment.id,Replicate,Bait,Prey,counts,Length")
      }

      Length <- NULL
      Replicate <- NULL

      # columns drop
      drops <- c("Replicate","Length", "counts")

      #list of feature matrices
      f_matrices <- list()

      if(cPASS_m) {
        x <- datInput[ ,!(colnames(datInput) == "Length")]
        f_matrices[["cPASS_m"]] <- cPASS_m(x)
      }


      if(MiST_m) {
        x <- datInput
        f_matrices[["MiST_m"]] <- MiST_m(x)
      }

      if(HG_m) {
        x <- datInput
        f_matrices[["HG_m"]] <- HG_m(x)
      }

      if(WM_m) {
        x <- datInput[ ,!(names(datInput) %in% drops)]
        f_matrices[["WM_m"]] <- WM_m(x)
      }

      if(PCC_m) {
        x <- datInput
        f_matrices[["PCC_m"]] <- PCC_m(x)
      }

      if(Dice_m) {
        x <- datInput[ ,!(names(datInput) %in% drops)]
        f_matrices[["Dice_m"]] <- Dice_m(x)
      }
      if(Jaccard_m) {
        x <- datInput[ ,!(names(datInput) %in% drops)]
        f_matrices[["Jaccard_m"]] <- Jaccard_m(x)
      }
      if(Overlap_m) {
        x <- datInput[ ,!(names(datInput) %in% drops)]
        f_matrices[["Overlap_m"]] <- Overlap_m(x)
      }

      if(Zscore_m) {
        x <- datInput[ ,!(names(datInput) %in% drops)]
        f_matrices[["Zscore_m"]] <- Zscore_m(x, iter = iter)
      }

      datPPI <-
        Reduce(.merge.all, f_matrices)

      return(datPPI)

    }
