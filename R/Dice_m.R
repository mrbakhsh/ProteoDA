  # create a final scoring data
   .scored_df <- function(x, bn, pn) {
   b_p <- NULL
   . <- NULL

    s1 <-
      x %>%
      filter(.[,1] %in% bn) %>%
      select(`b_p`, 3)
    s2 <-
      x %>%
      filter(.[,2] %in% bn) %>%
      select(2,1,3:ncol(.)) %>%
      mutate(`b_p` = paste(.[,1], .[,2], sep = "~")) %>%
      select(`b_p`, 3)
    s3 <-
      x %>%
      filter(.[,1] %in% pn & .[,2] %in% pn) %>%
      select(`b_p`, 3) %>%
      rbind(s1,s2,.)

    return(s3)

   }


  # create an input data
  .cm <- function(x){
    . <- NULL
    Cn <- NULL
    protein <- NULL

    x <-
      x %>%
      gather("key", "protein", 2:3) %>%
      select(-2) %>%
      mutate(Cn = 1)
    x <- x[!duplicated(x), ]
    x <-
      x %>%
      spread(`protein`, `Cn`) %>%
      replace(is.na(.),0) %>%
      select(-1)
      return(x)

  }




  #' Dice_m
  #' @param datInput Data frame with column names:
  #' Experiment.id, Bait, and Prey.
  #' @return Data frame containing bait-prey pairs with Dice coefficient score
  #' greater than zero.
  #' @references Zhang,B. et al. (2008) From pull-down data to protein
  #' interaction networks and complexes with biological
  #' relevance. Bioinformatics.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @description This function applies Dice coefficient to score instances
  #' (e.g., bait-prey interactions (BPIs) in the data.frame.The Dice
  #' coefficient was first applied by Zhang et al., 2008
  #' to score interactions between all identified proteins (baits and preys)
  #' in a given AP-MS experiment. This scoring algorithm is
  #' based on matrix model that treats bait–prey and prey–prey interactions
  #' equally.
  #' @export
  #' @examples
  #' data("Testdata_1")
  #' Testdata_1 <- Testdata_1[, c(1,3:4)]
  #' scoredfile1 <- Dice_m(Testdata_1)

  Dice_m <- function(datInput) {


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

    Dice <- NULL
    Var1 <- NULL
    Var2 <- NULL
    . <- NULL

    datInput <-
      datInput[!duplicated(datInput),]

    # bait name
    bn <- unique(datInput$Bait)
    # prey name
    pn <- unique(datInput$Prey)

    # create a data input
    c_mat <- .cm(datInput)

    # compute Dice
    m <-
      as.matrix(proxy::simil(c_mat, by_rows = FALSE, method = "Dice"))
    m[lower.tri(m, diag=TRUE)] <- NA
    # convert to data frame
    m_d <-
      melt(m, na.rm=TRUE,value.name="Dice") %>%
      as.data.frame(.) %>%
      filter(Dice > 0) %>%
      mutate(`b_p` =paste(`Var1`, `Var2`, sep = "~"))

    fdata <-
      .scored_df(m_d, bn, pn)


    return(fdata)
  }
