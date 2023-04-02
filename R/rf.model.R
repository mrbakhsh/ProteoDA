
  .cm_f <- function(preds, real){
    cm <-
      table(preds, real)

    TP <- as.double(cm[1, 1])
    TN <- as.double(cm[2, 2])
    FP <- as.double(cm[1, 2])
    FN <- as.double(cm[2, 1])

    ACC <- (TP + TN) / (TP + TN + FP + FN)
    SE <- TP / (TP + FN)
    SP <- TN / (FP + TN)
    PPV <- TP / (TP + FP)
    F1 <- 2 * TP / (2 * TP + TP + FN)
    MCC <- (TP * TN - FP * FN) /
      sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

    cm_result <- list()

    cm_result[["cm"]] <- cm
    cm_result[["ACC"]] <- ACC
    cm_result[["SE"]] <- SE
    cm_result[["SP"]] <- SP
    cm_result[["PPV"]] <- PPV
    cm_result[["F1"]] <- F1
    cm_result[["MCC"]] <- MCC


    return(cm_result)

  }






  #' rf.model
  #' @param features A data frame with bait-prey or pre-prey associations in the
  #'   first column, and features to be passed to the classifier in the remaining
  #'   columns.
  #' @param gd  A gold reference set including true associations with class labels
  #'   indicating if such PPIs are positive or negative.
  #' @param cv_fold  Number of partitions for cross-validation; defaults to 5.
  #' @param plots Logical value, indicating whether to plot the performance of the
  #'   learning algorithm using k-fold cross-validation; defaults to FALSE.If the
  #'   argument set to TRUE, plots will be saved in the temp() directory. These
  #'   plots are : \itemize{ \item{pr_plot} - Precision-recall PLOT
  #'   \item{roc_plot} - ROC plot \item{radar_plot} - Radar plot showing accuracy,
  #'   F1-score , positive predictive value (PPV), sensitivity (SE) and MCC. }
  #' @param filename character string, indicating the output filename as an pdf
  #'   object. Defaults to plotsRF.pdf.
  #'
  #' @return predicted_interactions- Prediction scores for input dataset (i.e.,
  #' features) cm_ouput \itemize{ \item{cm} - A confusion matrix. \item{ACC} -
  #' Accuracy. \item{SE} - Sensitivity. \item{SP} - Specificity. \item{PPV} -
  #' Positive Predictive Value. \item{F1} - F1-score. \item{MCC} - Matthews
  #' correlation coefficient } roc_auc - area under ROC curve pr_auc - area under
  #' PR curve
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @importFrom pROC roc
  #' @importFrom PRROC pr.curve
  #' @importFrom pROC ggroc
  #' @importFrom dplyr inner_join
  #' @importFrom stats predict
  #' @importFrom ggplot2 aes
  #' @importFrom ggplot2 ggplot
  #' @importFrom ggplot2 annotate
  #' @importFrom ggplot2 coord_equal
  #' @importFrom ggplot2 element_blank
  #' @importFrom ggplot2 element_line
  #' @importFrom ggplot2 element_text
  #' @importFrom ggplot2 geom_abline
  #' @importFrom ggplot2 scale_y_continuous
  #' @importFrom ggplot2 theme
  #' @importFrom ggplot2 theme_bw
  #' @importFrom ggplot2 unit
  #' @importFrom ggplot2 xlab
  #' @importFrom ggplot2 geom_line
  #' @importFrom ggplot2 ylab
  #' @importFrom grDevices dev.off
  #' @importFrom grDevices pdf
  #' @importFrom grDevices recordPlot
  #' @description This function uses Random Forest classifier to assigns each
  #'   bait-prey pair a confidence score, indicating the level of support for that
  #'   pair of proteins to interact. This function also computes and plots
  #'   Precision-Recall (PR), ROC curve, and randar plot including Sensitivity
  #'   (SE), Specificity (SP), Accuracy (ACC), Positive Predictive Value (PPV),
  #'   F1-acore, and Matthews correlation coefficient (MCC) to evaluate the
  #'   performance of the classifier on k-fold cross-validated data.
  #' @export
  #' @examples
  #' data("Testdata_1")
  #' # Score using spoke model
  #' scdata <- score_APMS(Testdata_1, cPASS_m = TRUE, MiST_m = TRUE)
  #' # create random train data
  #' library(dplyr)
  #' t <- sample_n(scdata, 500)
  #' t$label <- # random lables
  #' sample(c("Positive","Negative"), size = nrow(t), replace = TRUE)
  #' t <- t[, c(1,14)] #select pair with the class labels
  #' # perfrorm prediction
  #' pred <- rf.model(scdata, t, cv_fold = 5,  plots = FALSE)




  rf.model <- function (features,
                        gd,
                        cv_fold = 5,
                        plots = FALSE,
                        filename=file.path(tempdir(), "plotsRF.pdf")) {


    if (!is.data.frame(features)) features <- as.data.frame(features)

    if (all(colnames(features) != "b_p") == TRUE) {
      stop("b_p is missing from feature...")
    }

    if (all(colnames(gd) != "label") == TRUE) {
      stop("label attribute is absent from the gold_standard data")
    }

    if (all(colnames(gd) != "b_p") == TRUE) {
      stop("b_p is missing from gold_standard...")
    }
    if (!is.character(gd$label)) {
      stop("Label column must include character vector")
    }
    . <- NULL
    Negative <- NULL
    Positive <- NULL
    V1 <- NULL
    V2 <- NULL
    b_p <- NULL



    t.data <-
      features

    `%nin%` <- Negate(`%in%`)
    p.data <-
      features %>% filter(b_p %nin% gd$b_p)


    # training data preparation
    x_train <-
      inner_join(t.data, gd, by = "b_p") %>% select(-b_p)

    # change charactor to factor and reverse the factor level
    x_train$label <- # chanage charactor to factor
      as.factor(as.character(x_train$label))

    # check the factor level
    x <- unique(x_train$label)
    x <- data.frame(x, as.integer(x))
    x <-
      x %>% filter(x == "Positive") %>% .$as.integer.x.
    if (x == 2) {
      x_train$label <- factor(x_train$label,
                              levels = rev(levels(x_train$label))
      )
    }

    # run the prediction
    clcol <- which(names(x_train) == "label")
    xtr <- x_train[,-clcol]
    ytr <- x_train$label
    mx = dim(xtr)[1]
    rfpred <- x_train$label
    prob = matrix(nrow = length(ytr), ncol = 2)
    index = rep(1:cv_fold, nrow(xtr))
    ind = index[1:nrow(xtr)]


    for (k in 1:cv_fold) {
      cat(".")
      xcal <- xtr[ind != k, ]
      ycal <- ytr[ind != k]
      xtest <- xtr[ind == k, ]


      rfout <- randomForest::randomForest(ycal~.,
                                          data = data.frame(xcal, ycal),
                                          importance = FALSE)
      rfpred[ind == k] = predict(rfout, xtest, type = "response")
      prob[ind == k, ] = predict(rfout, xtest, type = "prob")

    }
    prob <- as.data.frame(prob)
    colnames(prob) <-
      c("Positive", "Negative")



    # predicted interactions
    prob_allD = predict(rfout, features[,-1], type = "prob")

    prob_allD <-
      cbind(features$b_p, prob_allD)
    prob_allD <-
      as.data.frame(prob_allD) %>% rename(b_p=V1) %>%
      select(-Negative)%>%
      mutate(Positive = as.numeric(Positive))


    # Performance evaluation
    cm_out <- .cm_f(rfpred, ytr)
    # roc analysis
    roc_c <-
      roc(ytr, prob$Positive)
    # pr analysis
    pr_c <-
      pr.curve(
        scores.class0 = prob$Positive[ytr == "Positive"],
        scores.class1 =  prob$Positive[ytr == "Negative"],
        curve = TRUE)

    if (plots) {

      pdf(filename)

      # Generate plot for cm result
      df <-
        data.frame(rbind(rep(1,6),
                         rep(0,6), cbind(cm_out$ACC,cm_out$SE,cm_out$SP,
                                         cm_out$PPV,cm_out$F1,cm_out$MCC)))
      colnames(df) <- c("ACC","SE","SP","PPV","F1", "MCC")

      raderplot <-
        fmsb::radarchart(df,cglty = 2, pfcol = c("#99999980"),
                   cglcol = "blue",pcol = 2,plwd = 2, plty = 1)

      p <- recordPlot()



      # roc curve
      roc_plot <-
        ggroc(roc_c, legacy.axes = TRUE) +
        geom_abline(
          slope = 1, intercept = 0,
          linetype = "dashed", alpha = 0.7, color = "darkgreen"
        ) + coord_equal() + theme_bw() +
        theme(
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        ) +
        theme(axis.ticks.length = unit(.5, "cm")) +
        theme(text = element_text(size = 14, color = "black")) +
        theme(axis.text = element_text(size = 12, color = "black")) +
        xlab("False Positive Rate (1-Specificity)") +
        ylab("True Positive Rate (Sensitivity)") +
        annotate("text", x=0.25, y=0.8,
                 label= paste("AUC:",round(roc_c[["auc"]],2)))


      # PR_curve
      PR_Object <- as.data.frame(pr_c$curve)

      pr_plot <-
        ggplot(PR_Object, aes(x = V1, y = V2)) +
        geom_line() +
        theme_bw() + scale_y_continuous(limits = c(0, 1)) +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
        theme(axis.ticks.length = unit(.5, "cm")) +
        theme(text = element_text(size = 14, color = "black")) +
        theme(axis.text = element_text(size = 12, color = "black")) +
        xlab("Recall") + ylab("Percision") +
        annotate("text", x=0.25, y=0.25,
                 label= paste("AUC:",round(pr_c[["auc.integral"]],2)))

      # Make plots.
      plot_list <- list()
      plot_list$roc <- roc_plot
      plot_list$pr <- pr_plot
      plot_list$rader <- p


      for (i in seq_along(plot_list)) {
        print(plot_list[[i]])
      }
      dev.off()

    }


    output <- list()
    output[["predicted_interactions"]] <- prob_allD
    output[["cm_ouput"]] <- cm_out
    output[["roc_auc"]] <- roc_c$auc
    output[["pr_auc"]] <- pr_c$auc.integral


    return(output)


  }






