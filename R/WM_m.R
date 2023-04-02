  #' WM_m
  #' @param datInput Data frame with column names:
  #' Experiment.id, Bait, and Prey.
  #' @return Data frame containing bait-prey pairs with
  #' k (i.e.,number of co-purifications) greather or equal to 1
  #' & logHG (i.e., -1*log(P-val) of the hypergeometric test).
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @references Drew,K. et al. (2017) Integration of over 9,000 mass
  #' spectrometry experiments builds a global map of human
  #' protein complexes. Mol. Syst. Biol.
  #' @importFrom dplyr mutate
  #' @description This function computes the weighted matrix model for
  #' instances (e.g., bait-prey interactions (BPIs)) in the data.frame.
  #' The output of the weighted matrix model includes the number of experiments
  #' for which the pair of proteins is co-purified (i.e., k) and
  #' -1*log(P-value) of the hypergeometric test (i.e., logHG) given the
  #' experimental overlap value, each protein's total number of observed
  #' experiments, and the total number of experiments (Drew et al., 2017).
  #' This scoring algorithm is based on matrix model that treats
  #' bait–prey and prey–prey interactions equally.
  #' @export
  #' @examples
  #' data("Testdata_1")
  #' Testdata_1 <- Testdata_1[, c(1,3:4)]
  #' scoredfile1 <- WM_m(Testdata_1)


  WM_m <- function(datInput) {


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
    HG <- NULL
    N <- NULL
    ObsExp <- NULL
    UniprotID <- NULL
    Var1 <- NULL
    Var2 <- NULL
    k <- NULL
    m <- NULL

    datInput <-
      datInput[!duplicated(datInput),]

    # bait name
    bn <- unique(datInput$Bait)
    # prey name
    pn <- unique(datInput$Prey)

    # CO-PURIFIED info
    c_mat <- .cm(datInput)
    c_mat <- as.matrix(c_mat)
    c_mat.adj <- t(c_mat) %*% c_mat
    diag(c_mat.adj) <- 1
    c_mat.adj[lower.tri(c_mat.adj, diag=TRUE)] <- NA
    CoPurifed <-
      melt(c_mat.adj, na.rm=TRUE, value.name="k") %>%
      as.data.frame(.) %>% filter(k >= 1) %>%
      mutate(`b_p` =paste(`Var1`, `Var2`, sep = "~"))

    # Sum per protein
    protInt <-
      data.frame(ObsExp = colSums(c_mat[,]),
                 stringsAsFactors = FALSE) %>%
      rownames_to_column("UniprotID")

    m_d <-
      CoPurifed %>%
      rename(UniprotID=`Var1`) %>%
      left_join(., protInt, by = "UniprotID") %>%
      rename(n=ObsExp) %>% rename(`Var1`=UniprotID) %>%
      rename(UniprotID=`Var2`) %>%
      left_join(., protInt,by = "UniprotID") %>%
      dplyr::rename(m=ObsExp) %>%
      rename(`Var2`=UniprotID) %>%
      mutate(N = max(datInput$Experiment.id)) %>%
      mutate(HG = 1- phyper(k-1, n, N-n, m)) %>%
      mutate(logWM = -1 * log(HG)) %>%
      mutate(b_p = paste(`Var1`,`Var2`, sep = "~")) %>% select(1,2,9,4)

    fdata <-
      .scored_df(m_d, bn, pn)

    return(fdata)
    }
