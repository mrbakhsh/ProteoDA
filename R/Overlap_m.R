  .f <-  function(x, y)
    (sum(x == 1 & y == 1)^2)/(sum(x == 1) * sum(y == 1))

  #' Overlap_m
  #' @param datInput Data frame with column names:
  #' Experiment.id, Replicate, Bait, and Prey.
  #' @return Data frame containing bait-prey pairs with overlap score
  #' greather than zero.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @references Nepusz,T. et al. (2012) Detecting overlapping protein
  #' complexes in protein-protein interaction networks. Nat. Methods.
  #' @description This function computes the overlap similarity scores
  #' for instances (e.g., bait-prey interactions (BPIs)) in
  #' the data.frame. This scoring algorithm is
  #' based on matrix model that treats bait–prey and prey–prey interactions
  #' equally.
  #' @export
  #' @examples
  #' data("Testdata_1")
  #' Testdata_1 <- Testdata_1[, c(1,3:4)]
  #' scoredfile1 <- Overlap_m(Testdata_1)



  Overlap_m <- function(datInput) {


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
    Overlap <- NULL
    b_p <- NULL
    Var2 <- NULL
    Var1 <- NULL


    datInput <-
      datInput[!duplicated(datInput),]

    # bait name
    bn <- unique(datInput$Bait)
    # prey name
    pn <- unique(datInput$Prey)

    # create a matrix data
    c_mat <- .cm(datInput)

    m <-
      as.matrix(proxy::dist(t(c_mat),.f))
    m[lower.tri(m, diag=TRUE)] <- NA
    # convert to data frame
    m_d <-
      melt(m, na.rm=TRUE,value.name="Overlap") %>%
      as.data.frame(.) %>%
      filter(Overlap > 0) %>%
      mutate(`b_p` =paste(`Var1`, `Var2`, sep = "~"))

    fdata <-
      .scored_df(m_d, bn, pn)


    return(fdata)
  }
