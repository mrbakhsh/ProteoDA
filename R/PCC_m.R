  #' PCC_m
  #' @param datInput Data frame with column names:
  #' Experiment.id, Replicate, Bait, Prey, counts (i.e., summed spectral
  #' counts or log2-transformed summed MS1-intensities) and  prey protein
  #' length.
  #' @return Data frame containing bait-prey pairs with Pearson correlation
  #' scores greather than zero.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @importFrom stats cor
  #' @description This function computes the correlation
  #' between instances (e.g., bait-prey interactions (BPIs)). This scoring
  #' algorithm is based on matrix model that only considers prey-prey
  #' interaction.
  #' @export
  #' @examples
  #' data("Testdata_1")
  #' # score file containing spectral counts
  #' scoredfile1 <- PCC_m(Testdata_1)
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
  #' # 2. then we apply PCC
  #' scoredfile2 <- PCC_m(Testdata_2)

  PCC_m <- function(datInput) {

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
    Abundance <- NULL
    Bait <- NULL
    Experiment.id <- NULL
    Length <- NULL
    M3D_n <- NULL
    Prey <- NULL
    Replicate <- NULL
    Var1 <- NULL
    Var2 <- NULL
    counts <- NULL
    key <- NULL
    protein <- NULL
    x1 <- NULL
    PCC <- NULL


    # bait name
    bn <- unique(datInput$Bait)
    # prey name
    pn <- unique(datInput$Prey)


    # Normalize by protein length (input for PCC)
    dN <-
      datInput %>%
      mutate(x1 = counts/Length) %>%
      group_by(Bait,Replicate) %>%
      mutate(x1 = scale(x1, center = FALSE,
                        scale = sum(counts))) %>%
      mutate(M3D_n = scale(x1, center = FALSE,
                           scale = sum(x1))) %>%
      ungroup(.) %>%
      group_by(Bait, Prey) %>%
      summarise(# Average counts for bait_prey across replicates
        Experiment.id = unique(Experiment.id),
        Prey =  unique(Prey),
        Abundance = mean(M3D_n), .groups = 'drop') %>%
      select(`Experiment.id`,`Abundance`, `Prey`) %>%
      spread(Prey, Abundance) %>%
      as.data.frame(.) %>%
      set_rownames(.$Experiment.id)%>%
      select(-`Experiment.id`)%>%
      replace(is.na(.), 0)


    cor_m <-
      cor(dN, method = "pearson",
          use = "pairwise.complete.obs")

    cor_m[upper.tri(cor_m, diag=TRUE)] <- NA
    Scored_file <-
      melt(cor_m, na.rm=TRUE, value.name="PCC")
    Scored_file <- Scored_file[, c(2,1,3)]
    colnames(Scored_file)[1:2] <- c("Var1", "Var2")

    dat <-
      Scored_file %>%
      filter(PCC > 0) %>%# remove negative PCC
      mutate(b_p = paste(Var1,Var2, sep = "~"))
    fdata <-
      .scored_df(dat, bn, pn)

    return(fdata)
  }


