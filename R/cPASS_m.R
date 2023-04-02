    #function for entropy

    .entropy <- function(xs) {
      p <- (xs + 1/length(xs)) / (sum(xs) + 1)
      ent <- sum(sapply(p, function(x) { -1*x*log(x, 2) }))
      return(ent)
    }


  #'cPASS_m
  #' @param datInput Data frame with column names:
  #' Experiment.id, Replicate, Bait, Prey, and counts (i.e., summed spectral
  #' counts or log2-transformed summed MS1-intensities).
  #' @return Data frame containing bait-prey pairs with average peptide
  #' spectrum match (PSMs), total PSMs, ratio total PSMs, Z-score, S-score,
  #' D-score and WD-score.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @references Huttlin,E.L. et al. (2015) The BioPlex Network: A Systematic
  #' Exploration of the Human Interactome. Cell.
  #' @importFrom dplyr mutate
  #' @importFrom dplyr n
  #' @importFrom magrittr %>%
  #' @importFrom dplyr filter
  #' @importFrom dplyr row_number
  #' @importFrom dplyr rename
  #' @importFrom dplyr group_by
  #' @importFrom dplyr left_join
  #' @importFrom dplyr summarise_each
  #' @importFrom dplyr funs
  #' @importFrom dplyr select
  #' @importFrom utils combn
  #' @importFrom tibble rownames_to_column
  #' @importFrom stats phyper
  #' @importFrom tidyr separate
  #' @importFrom tidyr gather
  #' @importFrom tidyr spread
  #' @importFrom magrittr set_rownames
  #' @importFrom dplyr arrange
  #' @importFrom dplyr distinct
  #' @importFrom dplyr summarise
  #' @importFrom dplyr ungroup
  #' @importFrom plyr desc
  #' @importFrom stats quantile
  #' @importFrom dplyr summarize
  #' @description This function applies Comparative Proteomic Analysis Software
  #' Suite (CompPASS) model to score instances (e.g., bait-prey interactions
  #' (BPIs) in the data.frame.This scoring algorithm is based on spoke model.
  #' @export
  #' @examples
  #' data("Testdata_1")
  #' # drop the length column
  #' Testdata_1 <- Testdata_1[, -6]
  #' # score file containing spectral counts
  #' scoredfile1 <- cPASS_m(Testdata_1)
  #'
  #'
  #' # score file containing MS1-intensities
  #' # 1. first we log transform the counts
  #' library(dplyr)
  #' data("Testdata_2")
  #' Testdata_2 <- Testdata_2[, -6]
  #' Testdata_2 <-
  #'   Testdata_2 %>%
  #'   group_by(Bait, Replicate) %>%
  #'   mutate(counts = log2(counts))
  #' # 2. then we apply cPASS
  #' scoredfile2 <- cPASS_m(Testdata_2)



  cPASS_m <- function(datInput) {

    clnames <-
      c("Experiment.id","Replicate", "Bait", "Prey", "counts")



    if(!is.data.frame(datInput)) datInput <- as.data.frame(datInput)


    if(!all(clnames %in% colnames(datInput))) {
      missCols <-
        setdiff(clnames,colnames(datInput))
      stop("columns are missing: ", paste(missCols, collapse = ", "))
    }


    if(is.numeric(datInput$counts) ==  FALSE){
      stop("counts must contain numerical values")
    }

    if(ncol(datInput) != 5){
      stop("Data frame must contain 5 columns including:
          Experiment.id,Replicate,Bait,Prey,counts")
    }

    ObsExp <- NULL
    Bait <- NULL
    Prey <- NULL
    Replicate <- NULL
    counts <- NULL
    MaxSpec <- NULL
    TotalPSM <- NULL
    intBait <- NULL
    . <- NULL
    meanPint <- NULL
    meandiff <- NULL
    freq_b <- NULL
    SD <- NULL
    WD_inner <- NULL
    WD_raw <- NULL



    # Total number of Bait
    N <-
      length(unique(datInput$Bait))

    # Drop the experiment.id
    datInput <-
      datInput %>%
      select(-"Experiment.id")


    # Number prey captured across runs
    p_sum <-
      unique(datInput[, c("Bait", "Prey")])
    p_sum <-
      p_sum %>%
      group_by(Prey) %>%
      summarise(p_sum = n())


    # Bait-prey pair captured across runs
    l <-
      unique(datInput[, c("Bait", "Prey")])

    l <-
      l %>%
      group_by(Bait, Prey) %>%
      summarise(l = n()) %>%
      mutate(b_p = paste(Bait, Prey, sep = "~")) %>%
      select(4,1,2,3)


    # Compute Entropy using maximum count
    getentropy <-
      datInput %>%
      group_by(Bait, Prey, Replicate) %>%
      summarise(Bait = unique(Bait),
                MaxSpec = max(counts)) %>%
      ungroup() %>%
      group_by(Bait,Prey) %>%
      summarise(Bait = unique(Bait),
                Entropy = .entropy(MaxSpec)) %>%
      ungroup() %>%
      as.data.frame() %>%
      mutate(b_p = paste(Bait, Prey, sep = "~")) %>%
      dplyr::select(4,1,2,3)


    # Average counts for bait_prey across replicates
    AvePSM <-
      datInput %>%
      group_by(Bait, Prey) %>%
      summarise(AvePSM = mean(counts)) %>%
      mutate(b_p = paste(Bait, Prey, sep = "~")) %>%
      dplyr::select(4,1,2,3)

    other_Stats <-
      datInput %>%
      group_by(Bait, Prey) %>%
      summarise(AvePSM = mean(counts)) %>%
      ungroup() %>%
      group_by(Prey) %>%
      mutate(TotalPSM = sum(AvePSM)) %>%
      mutate(Ratio_totalPSM =AvePSM/TotalPSM) %>%
      dplyr::mutate(Ratio = n()/N) %>%
      mutate(b_p = paste(Bait, Prey, sep = "~"))%>%
      dplyr::select(7,1,2,4:6)


    # Number of prey that interact with one bait only
    P_int_a_Bait <-
      AvePSM %>%
      group_by(Prey) %>%
      dplyr::summarise(intBait = n()) %>%
      filter(intBait == 1)


    datStat <-
      AvePSM %>%
      left_join(., p_sum, by = "Prey") %>%
      left_join(., l[, c("b_p", "l")], by = "b_p") %>%
      left_join(., getentropy[, c("b_p", "Entropy")], by = "b_p") %>%
      left_join(.,
                other_Stats[, c("b_p", "TotalPSM","Ratio_totalPSM",
                                "Ratio")], by = "b_p")%>%
      mutate(N = N) %>%
      mutate(freq_b =
               ifelse((Prey %in% P_int_a_Bait$Prey),N, N-p_sum))

    datStat1 <-
      AvePSM %>%
      dplyr::select(3,2,4) %>%
      spread(Bait, AvePSM) %>%
      as.data.frame(.) %>%
      set_rownames(.$Prey) %>%
      dplyr::select(-1)

    aveIntPrey <-
      rowSums(as.matrix(datStat1), na.rm = TRUE) / N



    Prey_Calc <-
      data.frame(
        Prey = names(aveIntPrey),
        meanPint = aveIntPrey,
        meandiff = rowSums((as.matrix(datStat1) - aveIntPrey)^2,
                           na.rm = TRUE)) %>%
      left_join(.,datStat[, c("Prey", "freq_b")], by = "Prey") %>%
      distinct(Prey, meanPint,meandiff,freq_b, .keep_all = TRUE) %>%
      mutate(SD = sqrt((meandiff + ((meanPint^2) * (freq_b)))/(N-1)))

    dff <-
      left_join(datStat, Prey_Calc[, c("Prey", "meanPint", "SD")], by = "Prey")

    Scored_file <-
      dff %>%
      mutate(S_score = sqrt((AvePSM) * (N) / (p_sum))) %>%
      mutate(Z_score = (AvePSM - meanPint) / (SD)) %>%
      mutate(D_score = sqrt((AvePSM) * (((N) / (p_sum))^l))) %>%
      mutate(WD_inner = (N / p_sum) * (SD / meanPint)) %>%
      mutate(WD_raw = sqrt(AvePSM * (WD_inner^l)))


    # Weighted WD score
    s <-
      Scored_file %>%
      arrange(desc(WD_raw))
    WD_raw.factor <-
      as.numeric(quantile(s$WD_raw, 0.98)[1])

    Scored_file$WD_score <-
      Scored_file$WD_raw /WD_raw.factor


    Scored_file <-
      as.data.frame(Scored_file[, c("b_p", "AvePSM",
                                    "TotalPSM","Ratio_totalPSM","Ratio",
                                    "S_score", "Z_score", "D_score",
                                    "WD_score", "Entropy")])
    return(Scored_file)
  }

