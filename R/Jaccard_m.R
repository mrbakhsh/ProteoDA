  #' Jaccard_m
  #' @param datInput Data frame with column names:
  #' Experiment.id, Bait, and Prey.
  #' @return Data frame containing bait-prey pairs with Jaccard Index
  #' greater than zero.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @description This function computes the Jaccard Index
  #' for instances (e.g., bait-prey interactions (BPIs)) in
  #' the data.frame. This scoring algorithm is based on matrix model that
  #' treats bait–prey and prey–prey interactions equally.
  #' @importFrom proxy simil
  #' @importFrom reshape2 melt
  #' @export
  #' @examples
  #' data("Testdata_1")
  #' Testdata_1 <- Testdata_1[, c(1,3:4)]
  #' scoredfile1 <- Jaccard_m(Testdata_1)



  Jaccard_m <- function(datInput) {

    if(!is.data.frame(datInput)) datInput <- as.data.frame(datInput)
    clnames <-
      c("Experiment.id","Bait", "Prey")
    if(!is.data.frame(datInput)) datInput <- as.data.frame(datInput)
    if(!all(clnames %in% colnames(datInput))) {
      missCols <-
        setdiff(clnames,colnames(datInput))
      stop("columns are missing: ", paste(missCols, collapse = ", "))
    }
    if(ncol(datInput) != 3){
      stop("Data frame must contain 3 columns including:
          Experiment.id,Bait,Prey")
    }

    . <- NULL
    Jaccard <- NULL
    b_p <- NULL
    Var2 <- NULL
    Var1 <- NULL


    datInput <-
      datInput[!duplicated(datInput),]

    # bait name
    bn <- unique(datInput$Bait)
    # prey name
    pn <- unique(datInput$Prey)

    # create a data input
    c_mat <- .cm(datInput)

    # compute Jaccard distance
    m <-
      as.matrix(proxy::simil(c_mat, by_rows = FALSE, method = "Jaccard"))
    m[lower.tri(m, diag=TRUE)] <- NA
    # convert to data frame
    m_d <-
      melt(m, na.rm=TRUE,value.name="Jaccard") %>%
      as.data.frame(.) %>%
      filter(Jaccard > 0) %>%
      mutate(`b_p` =paste(`Var1`, `Var2`, sep = "~"))

    fdata <-
      .scored_df(m_d, bn, pn)

    return(fdata)
  }
