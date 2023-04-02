
    .n_f <- function(x){
      x <- apply(x, 2, function(z) z/sum(z, na.rm = T))
      return(as.data.frame(x))
    }


  #' MiST_m
  #' @param datInput Data frame with column names:
  #' Experiment.id, Replicate, Bait, Prey, counts (i.e., summed spectral
  #' counts or log2-transformed summed MS1-intensities), and
  #' prey protein length.
  #' @return Data frame containing bait-prey pairs with  MiST score (abundance,
  #' reproducibility, specificity)
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @references Verschueren,E. et al. (2015) Scoring large‐scale affinity
  #' purification mass spectrometry datasets with MiST.
  #' Curr. Protoc. Bioinforma.
  #' @importFrom dplyr full_join
  #' @importFrom stats na.omit
  #' @importFrom dplyr n_distinct
  #' @description This function applies Mass spectrometry interaction
  #' STatistics (MiST) model to score instances (e.g., bait-prey interactions
  #' (BPIs) in the data.frame. The MiST score is a linear combination of prey
  #' quantity (abundance), abundance invariability across repeated experiments
  #' (reproducibility), and prey uniqueness relative to other baits
  #' (specificity).This scoring algorithm is based on spoke model.
  #' @export
  #' @examples
  #' data("Testdata_1")
  #' # score file containing spectral counts
  #' scoredfile1 <- MiST_m(Testdata_1)
  #'
  #'
  #' # score file containing MS1-intensities
  #' # 1. first we log transform the counts
  #' library(dplyr)
  #' data("Testdata_2")
  #' Testdata_2 <-
  #'     Testdata_2 %>%
  #'     group_by(Bait, Replicate) %>%
  #'     mutate(counts = log2(counts))
  #' # 2. then we apply MIST
  #' scoredfile2 <- MiST_m(Testdata_2)

    MiST_m <- function(datInput) {

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


    . <- NULL
    Abundance <- NULL
    Bait <- NULL
    Length <- NULL
    M3D_n <- NULL
    Prey <- NULL
    R1 <- NULL
    Replicate <- NULL
    b_p <- NULL
    counts <- NULL
    nreplicate <- NULL
    x1 <- NULL
    Specificity <- NULL
    Experiment.id <- NULL
    NPSM <- NULL



    # Normalize by protein length (input for )
    datInput_n <-
      datInput %>%# Drop the Ex.id as we don't need it for MiST scoring
      select(-`Experiment.id`) %>%
      mutate(x1 = counts/Length) %>%
      group_by(Bait,Replicate) %>%
      mutate(x1 = scale(x1, center = FALSE,
                        scale = sum(counts))) %>%
      mutate(M3D_n = scale(x1, center = FALSE,
                           scale = sum(x1))) %>%
      select(`Replicate`, `Bait`, `Prey`,`M3D_n`)



    # Get abundance (same as average in cPASS)
    getAbundance <- # same as average
      datInput_n %>%
      group_by(Bait, Prey) %>%
      summarise(Abundance = mean(M3D_n)) %>%
      mutate(b_p = paste(Bait, Prey, sep = "~")) %>%
      ungroup(.) %>%
      select(`b_p`,`Abundance`)


    # Get reproducibility
    # Number of replicates per bait-prey
    R <-
      datInput_n %>%
      group_by(Bait, Prey) %>%
      dplyr::summarise(nreplicate = n_distinct(Replicate)) %>%
      mutate(b_p = paste(Bait, Prey, sep = "~")) %>%
      ungroup(.) %>%
      select(`b_p`,`nreplicate`)


    getReproducibility <-
      datInput_n %>%
      mutate(b_p = paste(Bait, Prey, sep = "~")) %>%
      left_join(., R) %>%
      group_by(Bait, Prey) %>%
      mutate(NPSM = ifelse(nreplicate >1,(M3D_n/sum(M3D_n)), M3D_n)) %>%
      replace(is.na(.), 0) %>%
      ungroup(.) %>%
      select(`b_p`,`nreplicate`,`NPSM`) %>%
      mutate(R1 = with(.,ifelse(NPSM != 1 & NPSM>0 & nreplicate >1,
                                NPSM * log2(NPSM),
                                ifelse(NPSM == 1 & nreplicate >1,
                                       ((1-1e-10)*log2(1-1e-10)),
                                       ifelse(nreplicate == 1, 0, NPSM))))) %>%
      group_by(b_p) %>%
      summarise(Reproducibility = sum(R1)/log2(1/mean(nreplicate))) %>%
      replace(is.na(.), 0) %>%
      ungroup(.)



   # Get specificity
    getSpecificity <-
      datInput_n %>%
      group_by(Bait,Prey) %>%
      summarise(Abundance =  mean(M3D_n)) %>%
      tidyr::spread(., Prey, Abundance) %>%
      ungroup(.) %>%
      as.data.frame(.) %>%
      magrittr::set_rownames(.$Bait) %>%
      dplyr::select(-`Bait`)%>%
      .n_f(.) %>%
      tibble::rownames_to_column("Bait") %>%
      gather("Prey", "Specificity", 2:ncol(.)) %>%
      na.omit(.)%>%
      mutate(b_p = paste(Bait, Prey, sep = "~")) %>%
      select(`b_p`,`Specificity`)


    # Create the final data.frame

    Scored_file <-
      full_join(getAbundance,getReproducibility, by = "b_p") %>%
      full_join(., getSpecificity, by = "b_p")

    return(Scored_file)
  }
