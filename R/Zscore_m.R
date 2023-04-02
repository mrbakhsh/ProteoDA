  .Random_f <- function(x, iter = iter, obs_data){

    Prey<- NULL
    .<- NULL
    Var1<- NULL
    Var2<- NULL
    b_p <- NULL

    datalist = list()
    for (i in 1:iter){
      # create suffled data
      rand_data <-
        x %>% mutate(Prey=sample(Prey))
      # create adj matrix data
      cm_i <- .cm(rand_data)
      cm_i <- as.matrix(cm_i)
      cm_i.adj <- t(cm_i) %*% cm_i
      diag(cm_i.adj) <- 1
      cm_i.adj[lower.tri(cm_i.adj, diag=TRUE)] <- NA
      rand_df <-
        melt(cm_i.adj, na.rm=TRUE, value.name="rand_int") %>%
        as.data.frame(.) %>%
        mutate(`b_p` = paste(`Var1`, `Var2`, sep = "~")) %>%
        filter(`b_p` %in% obs_data$b_p) %>%
        set_rownames(.$`b_p`) %>% select(3)

      datalist[[i]] <- rand_df # add it to your list
    }
    r_d <-
      bind_cols(datalist)

    average <- # get the mean of intersection
      data.frame(apply(r_d,1,mean), stringsAsFactors = FALSE)
    colnames(average) <- "mean"

    sd_c <-  # get sd intersection
      data.frame(apply(r_d,1,mean), stringsAsFactors = FALSE)
    colnames(sd_c) <- "std"
    datR_PPI <-
      data.frame(cbind(average, sd_c),
                 stringsAsFactors = FALSE) %>%
      rownames_to_column("b_p")

    return(datR_PPI)
  }




  #' Zscore_m
  #' @param datInput Data frame with column names:
  #' Experiment.id, Bait, and Prey.
  #' @param iter  Number of bootstraps to execute to estimate z scores.
  #' Defaults to 3.
  #' @return Data frame containing bait-prey pairs with normalized Z-score.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @description This function applies Z-score to score instances
  #' (e.g., bait-prey interactions (BPIs) in the data.frame. This scoring
  #' algorithm is based on matrix model that treats bait–prey and
  #' prey–prey interactions equally.
  #' @importFrom dplyr bind_cols
  #' @export
  #' @examples
  #' data("Testdata_1")
  #' Testdata_1 <- Testdata_1[, c(1,3:4)]
  #' scoredfile1 <- Zscore_m(Testdata_1, iter = 2)


  Zscore_m <- function(datInput, iter = 3) {

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
    obs_int<- NULL
    Var1<- NULL
    Var2<- NULL
    std<- NULL
    Zscore<- NULL
    b_p<- NULL

    datInput <-
      datInput[!duplicated(datInput),]

    # bait name
    bn <- unique(datInput$Bait)
    # prey name
    pn <- unique(datInput$Prey)

    # create a matrix data
    c_mat <- .cm(datInput)
    c_mat <- as.matrix(c_mat)
    c_mat.adj <- t(c_mat) %*% c_mat
    diag(c_mat.adj) <- 1
    c_mat.adj[lower.tri(c_mat.adj, diag=TRUE)] <- NA
    obs_data <-
      melt(c_mat.adj, na.rm=TRUE, value.name="obs_int") %>%
      as.data.frame(.) %>% filter(obs_int > 0) %>%
      mutate(`b_p` =paste(`Var1`, `Var2`, sep = "~")) %>% select(4,3)



    # create random set and get mean and sd
     R_df <- .Random_f(datInput, iter = iter, obs_data)

     sfile <-
      left_join(obs_data, R_df, by = "b_p") %>%
      mutate(Zscore = (obs_int - mean)/std) %>%
      na.omit(.) %>%
      filter(Zscore != "Inf" & Zscore != "-Inf") %>%
      mutate(N_Zscore = (Zscore - min(Zscore))/
                           (max(Zscore)-min(Zscore))) %>%
       separate(`b_p`, c("Var1","Var2"), sep = "~", remove = FALSE) %>%
       select(2,3,8,1)


     fdata <-
      .scored_df(sfile, bn, pn)


    return(fdata)
  }
