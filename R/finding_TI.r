#' =========================================================================
#' finding_TI   
#' -------------------------------------------------------------------------
#' finding_TI flags potential candidates for transcription interference
#' 
#' finding_TI uses 'score_fun_ave' to make groups by the mean of "probe_TI".
#' "TI" is added to the "flag" column. TI is characterized by relative
#' intensities at time points later than "0". The rowRanges need to contain at
#' least "ID", "probe_TI" and "position_segment"!
#'
#' @param inp SummarizedExperiment: the input.
#' @param cores integer: the number of assigned cores for the task
#' @param pen numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer fragments. Advised to be kept at 10.
#' Default is 10.
#' @param thrsh numeric: an internal parameter that allows fragments with a
#' certain amount of IDs with higher relative intensities at time points later
#' than "0" to be flagged as "_TI_". Higher values result in fewer candidates.
#' -0.5 is 25 %, 0 is 50%, 0.5 is 75%. Advised to be kept at 0.5.
#' Default is 0.5.
#' @param add integer: range of nucleotides before and after a potential TI
#' event wherein IDs are fitted with the TI fit.
#' 
#' @return The SummarizedExperiment object: with "_TI_" added to the flag
#' column.
#'       
#' @examples
#' data(preprocess_minimal)
#' finding_TI(inp = preprocess_minimal, cores = 2, pen = 10, thrsh = 0.5,
#' add = 1000)
#' 
#' @export


finding_TI <-
  function(inp, cores, pen = 10, thrsh = 0.5, add = 1000) {
    num_args <- c(pen, thrsh, add)
    names(num_args) <- c("pen", "thrsh", "add")
    assert(all(is.numeric(num_args)),
           paste0("one of the following arguments is not numeric: ",
                  paste0(names(num_args),collapse = ", ")))
    registerDoMC(cores) # cores for DoMC
    #order the input
    inp <- inp_order(inp)
    #make the tmp_df
    tmp_df <- inp_df(inp, "ID", "probe_TI", "position_segment")
    # makes a vector of all position segments (S_1,S_2,...)
    unique_seg <- unique(tmp_df$position_segment)

    frags <-
      foreach(k = seq_along(unique_seg)) %dopar% {
        # the loop iterates over each unique segment

        # only the part of the tmp_df that responds to the respective
        # segment is picked
        corr_IDs <- tmp_df[tmp_df$position_segment == unique_seg[k], "ID"]
        section <- tmp_df[match(corr_IDs, tmp_df$ID), ]
        # best_frags collects all scores that the dp is referring to
        best_frags <- c()
        # best names collects the names, and its last element is returned
        # as result
        best_names <- c()

        for (i in seq_len(nrow(section))) {
          # the loop iterates over each value in the segment
          # this part always goes from position 1 to the referred position
          # 1:1,1:2...
          tmp_score <-
            score_fun_ave(section[seq_len(i), "probe_TI"],
                          section[seq_len(i), "ID"], 0, 0) #outlier penalty and
          #allowed outliers are set to 0!
          tmp_name <- names(tmp_score)
          # in this loop all smaller parts are scored eg (i = 4) 4:4,3:4 and 2:4
          # they are then combined with the former score eg 1,2,3|4, 1,2|3,4...
          if (i > 1) {
            for (j in i:2) {
              tmp_probe_TI <- section[j:i, "probe_TI"]
              tmp_ID <- section[j:i, "ID"]
              tmp <-
                score_fun_ave(tmp_probe_TI, tmp_ID, 0, 0) + pen + best_frags[j - 1]
              #penalty for a new fragment and former scores are added
              tmp_score <- c(tmp_score, tmp) # the score is cached
              tmp_n <-
                paste0(best_names[j - 1], "|", names(tmp)) # the new fragment
              #is pasted to its corresponding former fragment
              tmp_name <- c(tmp_name, tmp_n) # the names is cached
            }
          }
          # from the first score eg 1:4 and the smaller scores from the loop
          # 1,2,3|4; 1,2|3,4; 1,2,3|4 the smallest is chosen and passed to
          # best_frags and best_names for the next iteration
          pos <-
            which(tmp_score == min(tmp_score))[1] # lowest score is collected
          tmp_score <- tmp_score[pos]
          tmp_name <- tmp_name[pos]
          best_frags <- c(best_frags, tmp_score)
          best_names <- c(best_names, tmp_name)
        }
        # the final result put into a list called frags
        best_names[length(best_names)]
      }
    rowRanges(inp)$flag <- gsub("_TI", "", rowRanges(inp)$flag)

    for (k in seq_along(frags)) {
      # the final result best_names[length(best_names)] is splitted into a
      # list by "|"
      na <- strsplit(frags[[k]], "\\|")[[1]]

      # the loop gives out scores by the mean
      for (i in seq_along(na)) {
        tmp_trgt <- strsplit(na[i], "_")[[1]][1] # gives IDs
        trgt <- strsplit(tmp_trgt, ",")[[1]]
        rows <-
          match(trgt, rowRanges(inp)$ID) # matches the row in the inp df
        pos_seg <- unique(rowRanges(inp)$position_segment[rows])
        first_pos <- rowRanges(inp)$position[rows[1]]
        if (all(decode(strand(inp))[rowRanges(inp)$position_segment == pos_seg] == "-", na.rm = TRUE)) {
          ID_before <-
            rowRanges(inp)$ID[which(rowRanges(inp)$position_segment == pos_seg)][
              na.omit(match(c((first_pos + add):(first_pos + 1)),
            rowRanges(inp)$position[which(rowRanges(inp)$position_segment == pos_seg)]))]
        } else {
          ID_before <-
            rowRanges(inp)$ID[which(rowRanges(inp)$position_segment == pos_seg)][na.omit(match(c((
              first_pos - add
            ):(
              first_pos - 1
            )), rowRanges(inp)$position[which(rowRanges(inp)$position_segment == pos_seg)]))]
        }
        rows_before <- match(ID_before, rowRanges(inp)$ID)
        rows <- c(rows_before, rows)
        score <- as.numeric(strsplit(na[i], "_")[[1]][2]) # gives the mean
        if (score >= thrsh) {
          # this threshold can be chosen in the function and is set
          # to 0.5 by default, meaning 75% of probes need to be TI
          # the range is -1(0%) (to 0(50%)) to 1(100%)
          rowRanges(inp)$flag[rows][grep("TI", rowRanges(inp)$flag[rows], invert = TRUE)] <-
            paste0(rowRanges(inp)$flag[rows][grep("TI", rowRanges(inp)$flag[rows],
                                            invert = TRUE)], "TI_") # all TI
          #candidates are flagged with "TI"
        }
      }
    }
    # the inp based df is returned as result
    inp
  }
