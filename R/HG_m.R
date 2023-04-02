  .min_sum <-
    function(x,y) sum(pmin(x,y))

  #' HG_m
  #' @param datInput Data frame with column names:
  #' Experiment.id, Replicate, Bait, Prey, and counts (i.e., summed spectral
  #' counts or log2-transformed summed MS1-intensities) anf prey protein
  #' length.
  #' @return Data frame containing bait-prey pairs with hypergeometric (HG)
  #' probability scores.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @references Guruharsha,K.G. et al. (2011) A protein complex network of
  #' Drosophila melanogaster. Cell.
  #' @importFrom dplyr bind_rows
  #' @importFrom dplyr summarise_all
  #' @description This function computes the hypergeometric probability
  #' between instances (e.g., bait-prey interactions (BPIs)) with
  #' incorporation of NSAF. This scoring algorithm is based on matrix model
  #' that only considers prey-prey interaction.
  #' @export
  #' @examples
  #' data("Testdata_1")
  #' # score file containing spectral counts
  #' scoredfile1 <- HG_m(Testdata_1)
  #'
  #'
  #' # score file containing MS1-intensities
  #' # 1. first we log transform the counts
  #' library(dplyr)
  #' data("Testdata_2")
  #' Testdata_2 <-
  #'  Testdata_2 %>%
  #'  group_by(Bait, Replicate) %>%
  #'  mutate(counts = log2(counts))
  #' # 2. then we apply HG
  #' scoredfile2 <- HG_m(Testdata_2)


  HG_m <- function(datInput) {

    if(!is.data.frame(datInput)){
      stop("Input data should be data.frame")
    }
    clnames <-
      c("Experiment.id","Replicate","Bait", "Prey", "counts", "Length")
    if(!is.data.frame(datInput)) datInput <- as.data.frame(datInput)
    if(!all(clnames %in% colnames(datInput))) {
      missCols <-
        setdiff(clnames,colnames(datInput))
      stop("columns are missing: ", paste(missCols, collapse = ", "))
    }
    if(ncol(datInput) != 6){
      stop("Data frame must contain 6 columns including:
            Experiment.id,Replicate,Bait,Prey,counts,Length")
    }

    . <- NULL
    Bait <- NULL
    Experiment.id <- NULL
    Var1 <- NULL
    Var2 <- NULL
    Length <- NULL
    NMinTn <- NULL
    NPSM <- NULL
    NS <- NULL
    NSAF <- NULL
    NormalNSAF <- NULL
    PPI <- NULL
    Prey <- NULL
    counts <- NULL
    ppiTN <- NULL
    tnA <- NULL
    tnB <- NULL
    SumNS <- NULL
    Tn <- NULL
    UniprotID <- NULL

    # bait name
    bn <- unique(datInput$Bait)
    # prey name
    pn <- unique(datInput$Prey)


    dN <-
      datInput %>%
      group_by(`Bait`, `Prey`) %>%
      summarise(# Average counts for bait_prey across replicates
        Experiment.id = unique(Experiment.id),
        Length =  unique(Length),
        counts = mean(counts), .groups = 'keep') %>%
      mutate(`NS` = `counts`/`Length`) %>%
      group_by(`Experiment.id`) %>%
      mutate(`SumNS` = sum(`NS`)) %>%
      mutate(`NSAF` = `NS`/`SumNS`) %>%
      group_by(`Experiment.id`) %>%
      mutate(`NormalNSAF` = `NSAF`/min(`NSAF`)) %>%
      mutate(`Tn` = as.integer(sqrt(`NormalNSAF`)))%>%
      select(`Experiment.id`,`Tn`,`Prey`) %>%
      spread(Prey, Tn) %>%
      as.data.frame(.)%>%
      set_rownames(.$Experiment.id)%>%
      select(-`Experiment.id`) %>%
      replace(is.na(.), 0)

    m <-
      as.matrix(proxy::dist(t(dN),.min_sum))
    m[lower.tri(m, diag=TRUE)] <- NA
    # convert to data frame
    m_d <-
      melt(m, na.rm=TRUE,value.name="ppiTN") %>%
      as.data.frame(.) %>%
      filter(ppiTN >= 1) %>%
      mutate(`b_p` =paste(`Var1`, `Var2`, sep = "~"))


    tnIntA <-
      m_d[, c("Var1", "ppiTN")]
    colnames(tnIntA) <-
      c("UniprotID", "ppiTN")
    tnIntB <-
      m_d[, c("Var2", "ppiTN")]
    colnames(tnIntB) <-
      c("UniprotID", "ppiTN")
    tnProtein <-
      bind_rows(tnIntA, tnIntB) %>%
      group_by(`UniprotID`) %>%
      summarise(minTn = sum(`ppiTN`))
    sumIntA <-
      tnProtein
    colnames(sumIntA) <-
      c("Var1", "tnA")
    sumIntB <-
      tnProtein
    colnames(sumIntB) <-
      c("Var2", "tnB")


    dat <-
      m_d %>%
      left_join(., sumIntA, by = "Var1") %>%
      left_join(., sumIntB, by = "Var2") %>%
      mutate(`b_p` = paste(`Var1`, `Var2`, sep = "~")) %>%
      mutate(`NMinTn` = sum(tnProtein$minTn)/2) %>%
      mutate(`HG` = -phyper(`ppiTN`, `tnA`, `NMinTn` - `tnB`,
                            `tnB`, lower.tail = FALSE, log.p = TRUE)) %>%
      select(1,2,8,4)

    fdata <-
      .scored_df(dat, bn, pn)

    return(fdata)

  }




