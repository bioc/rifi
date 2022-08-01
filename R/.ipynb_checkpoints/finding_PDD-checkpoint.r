# =========================================================================
# finding_PDD   Flags potential candidates for post transcription decay
# -------------------------------------------------------------------------
#'
#' 
#' 'finding_PDD' uses 'score_fun_linear_PDD' to make groups by the difference
#' to the slope. The slope is further checked for steepness to decide for PDD.
#' '_PDD_' is added to the 'flag' column.

#' Post transcription decay is characterized by a strong decrease of intensity
#' by position.

#' The rowRanges need to contain at least 'ID', 'intensity', 'position' and
#' 'position_segment'!
#'
#' @param inp SummarizedExperiment: the input.
#' @param cores integer: the number of assigned cores for the task
#' @param pen numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer fragments. Advised to be kept at 2.
#' Default is 2.
#' @param pen_out numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer possible outliers. Advised to be kept at 1.
#' Default is 1.
#' @param thrsh numeric: an internal parameter that allows fragments with slopes
#' steeper than the thrsh to be flagged with '_PDD_'. Higher values result in
#' fewer candidates. Advised to be kept at 0.001. Default is 0.001.
#' 
#' @return the SummarizedExperiment object: with "_PDD_" added to the flag
#' column.
#'       
#' @examples
#' data(preprocess_minimal)
#' finding_PDD(inp = preprocess_minimal, cores = 2, pen = 2,
#' pen_out = 1, thrsh = 0.001)
#' 
#' @export


finding_PDD <- function(inp, cores = 1, pen = 2, pen_out = 1, thrsh = 0.001) {
    stranded <- 1
    num_args <- c(pen, pen_out, thrsh)
    names(num_args) <- c("pen", "pen_out", "thrsh")
    assert(all(is.numeric(num_args)),
           paste0("one of the following arguments is not numeric: ",
                  paste0(names(num_args),collapse = ", ")))
    registerDoMC(cores)  #cores for DoMC
    #order the input
    inp <- inp_order(inp)
    #make the tmp_df
    tmp_df <- inp_df(inp, "ID", "position", "intensity", "position_segment")
    #revert the order in plus
    tmp_df <- tmp_df_rev(tmp_df, "+")
    #penalty is cached
    tmp_pen <- pen  
    tmp_pen_out <- pen_out
    tmp_df$intensity <- tmp_df$intensity / mean(tmp_df$intensity)
    # makes a vector of all position segments (S_1,S_2,...)
    unique_seg <- unlist(unique(tmp_df$position_segment))
    # the foreach loop iterates over each unique segment
    frags <- foreach(k = seq_along(unique_seg)) %dopar% {
        # only the part of the tmp_df that responds to the respective segment
        # is picked
        corr_IDs <- tmp_df[tmp_df$position_segment == unique_seg[k], "ID"]
        section <- tmp_df[match(corr_IDs, tmp_df$ID), ]
        # the penalties are dynamically adjusted to the mean of the intensity
        pen <- tmp_pen * mean(section$intensity)
        pen_out <- tmp_pen_out * mean(section$intensity)
        # best_frags collects all scores that the dp is referring to
        best_frags <- c()
        # best names collects the names, and its last element is returned as
        # result
        best_names <- c()
        if (nrow(section) > 2) {
            # only segments with more than two values are grouped...*

            for (i in 3:nrow(section)) {
                # the loop iterates over each value in the segment this part
                # always goes from position 1 to the referred position
                # 1:3,1:4...
                tmp_score <- score_fun_linear_PDD(section[seq_len(i),
                  "intensity"], section[seq_len(i),
                  "position"], section[seq_len(i), "ID"], pen_out, stranded)
                tmp_name <- names(tmp_score)
                # in this loop all smaller parts are scored eg (i = 6)
                # 6:6,5:6,4:6... they are then combined with the former score
                # eg 1,2,3,4,5|6, 1,2,3,4|5,6... only parts bigger than 6 are
                # accepted as three is the smallest possible fragment size
                if (i > 5) {
                  for (j in (i - 2):4) {
                    tmp_intensity <- section[j:i, "intensity"]
                    tmp_position <- section[j:i, "position"]
                    tmp_ID <- section[j:i, "ID"]
                    # penalty for a new fragment and former scores are added
                    tmp <- score_fun_linear_PDD(tmp_intensity, tmp_position,
                            tmp_ID, pen_out, stranded) + pen + best_frags[j - 3]
                    tmp_score <- c(tmp_score, tmp)  #the score is cached
                    # the new fragment is pasted to its corresponding former
                    # fragment
                    tmp_n <- paste0(best_names[j - 3], "|", names(tmp))
                    tmp_name <- c(tmp_name, tmp_n)  #the names is cached
                  }
                }
                # from the first score eg 1:6 and the smaller scores from the
                # loop 1,2,3,4,5|6, 1,2,3,4|5,6... the smallest is chosen and
                # passed to best_frags and best_names for the next iteration
                pos <- which(tmp_score == min(tmp_score))[1]  #lowest score is
                #collected
                tmp_score <- tmp_score[pos]
                tmp_name <- tmp_name[pos]
                best_frags <- c(best_frags, tmp_score)
                best_names <- c(best_names, tmp_name)
            }
        } else {
            # *...all segments with less than three values are grouped
            # automatically
            tmp_score <- score_fun_linear_PDD(section[, "intensity"],
                                              section[, "position"],
                section[, "ID"], pen_out, stranded)
            tmp_name <- names(tmp_score)
            best_names <- c(best_names, tmp_name)
        }
        # the final result put into a list called frags
        best_names[length(best_names)]
    }
    rowRanges(inp)$flag <- gsub("_PDD", "", rowRanges(inp)$flag)
    # this loop iterates over the segments
    for (k in seq_along(frags)) {
        # the single fragments are split by |
        na <- strsplit(frags[[k]], "\\|")[[1]]
        # the loop gives out scores by the slope
        for (i in seq_along(na)) {
            tmp_trgt <- strsplit(na[i], "_")[[1]][1]  #gives IDs
            trgt <- strsplit(tmp_trgt, ",")[[1]]  #wants to be numeric
            rows <- match(trgt, rowRanges(inp)$ID)  #matches the row in the inp
            score <- as.numeric(strsplit(na[i], "_")[[1]][2])  #gives the slope
            if (score > thrsh) {
                # this threshold can be chosen in the function and is set to
                # 0.001 by default, meaning the slope needs to be lower than
                # -0.001
                rowRanges(inp)$flag[rows] <- 
                  paste0(rowRanges(inp)$flag[rows], "PDD_")
                #all PDD candidates are flagged with 'PDD'
            }
        }
    }
    # the inp based df is returned as result
    inp
}
