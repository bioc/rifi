#' finding_PDD: flags potential candidates for post transcription decay.
#' 'finding_PDD' uses 'score_fun_linear_PDD' to make groups by the difference
#' to the slope.
#' Then the slope is checked for steepness to decide for PDD.
#' '_PDD_' is added to the 'flag' column.
#' Post transcription decay is characterized by a strong decrease of intensity
#' by position.
#' The data frame needs to contain at least 'ID', 'intensity', 'position' and
#' 'position_segment'!
#' @param probe data frame: the probe based data frame.
#' @param cores interger: the number of assigned cores for the task
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
#' @return the probe based data frame with the modified flag:
#'  \describe{
#'     \item{ID:}{The bin/probe specific ID}
#'     \item{position:}{The bin/probe specific position}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{intensity:}{The relative intensity at time point 0}
#'     \item{probe_TI:}{An internal value to determine which fitting model is
#'     applied}
#'     \item{flag:}{Information on which fitting model is applied}
#'     \item{postion_segment:}{The position based segment}
#'       }
#'       
#' @examples
#' data(preprocess_minimal)
#' finding_PDD(probe = preprocess_minimal$probe_df,  cores = 2, pen = 2,
#' pen_out = 1, thrsh = 0.001)
#' 
#' @export


# finding_PDD uses score_fun_linear_PDD to make groups by the difference to the
# slope. Then the slope is checked for steepness to decide for PDD. We are
# looking for steep decreasing slopes. The only input is the probe based df.
# pen is the penalty for the dp, pen_out is the outlier penalty, thrsh is the
# threshold for the slope. stranded should be TRUE if strand is not (all) NA.
finding_PDD <- function(probe, cores = 1, pen = 2, pen_out = 1, thrsh = 0.001) {
    stranded <- 1
    num_args <- list(cores, pen, pen_out, thrsh)
    names(num_args) <- c("cores", "pen", "pen_out", "thrsh")
    assert(all(unlist(lapply(num_args, FUN = function(x)
      (is.numeric(x) & length(x) ==
        1)))), paste0("'", names(which(unlist(lapply(num_args,
                                                     FUN = function(x)
                                                       (is.numeric(x) &
        length(x) == 1))) == FALSE))[1], "' must be numeric of length one"))
    assert(cores > 0, "'cores' must be a positive integer")
    req_cols_probe <- c("ID", "intensity", "position",
                        "position_segment", "strand")
    assert(all(req_cols_probe %in% colnames(probe)),
           paste0("'", req_cols_probe[which(!req_cols_probe %in%
           colnames(probe))], "' must be a column in 'probe'!"))
    registerDoMC(cores)  #cores for DoMC
    tmp_pen <- pen  #penalty is cached
    tmp_pen_out <- pen_out  #pen out is cached
    # a temporary df with ID, value (intensity), position and strand
    tmp_df <- data.frame(ID = probe$ID, val = probe$intensity,
                         position = probe$position,
                         seg = probe$position_segment)
    if (stranded == TRUE) {
        tmp_df$strand <- probe$strand
        # the positions are inverted for the '+' strand so that the slope is
        # seen as increasing by the scoring function.
        tmp_df[tmp_df$strand == "+", "position"] <-
          ((tmp_df[tmp_df$strand == "+",
            "position"]) - (tmp_df[tmp_df$strand == "+", ]
                            [nrow(tmp_df[tmp_df$strand == "+", ]),
                              "position"])) * -1
    }
    tmp_df <- na.omit(tmp_df)
    tmp_df$val <- tmp_df$val / mean(tmp_df$val)
    # makes a vector of all position segments (S_1,S_2,...)
    unique_seg <- unlist(unique(tmp_df$seg))
    # the foreach loop iterates over each unique segment
    frags <- foreach(k = seq_along(unique_seg)) %dopar% {
        # only the part of the tmp_df that responds to the respective segment
        # is picked
        corr_IDs <- tmp_df[tmp_df$seg == unique_seg[k], "ID"]
        section <- tmp_df[match(corr_IDs, tmp_df$ID), ]
        # the penalties are dynamically adjusted to the mean of the intensity
        pen <- tmp_pen * mean(section$val)
        pen_out <- tmp_pen_out * mean(section$val)
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
                tmp_score <- score_fun_linear_PDD(section[seq_len(i), "val"],
                                                  section[seq_len(i),
                  "position"], section[seq_len(i), "ID"], pen_out, stranded)
                tmp_name <- names(tmp_score)
                # in this loop all smaller parts are scored eg (i = 6)
                # 6:6,5:6,4:6... they are then combined with the former score
                # eg 1,2,3,4,5|6, 1,2,3,4|5,6... only parts bigger than 6 are
                # accepted as three is the smallest possible fragment size
                if (i > 5) {
                  for (j in (i - 2):4) {
                    tmp_val <- section[j:i, "val"]
                    tmp_position <- section[j:i, "position"]
                    tmp_ID <- section[j:i, "ID"]
                    # penalty for a new fragment and former scores are added
                    tmp <- score_fun_linear_PDD(tmp_val, tmp_position, tmp_ID,
                            pen_out, stranded) + pen + best_frags[j - 3]
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
            tmp_score <- score_fun_linear_PDD(section[, "val"],
                                              section[, "position"],
                section[, "ID"], pen_out, stranded)
            tmp_name <- names(tmp_score)
            best_names <- c(best_names, tmp_name)
        }
        # the final result put into a list called frags
        best_names[length(best_names)]
    }
    probe[, "flag"] <- gsub("_PDD", "", probe[, "flag"])
    # this loop iterates over the segments
    for (k in seq_along(frags)) {
        # the single fragments are split by |
        na <- strsplit(frags[[k]], "\\|")[[1]]
        # the loop gives out scores by the slope
        for (i in seq_along(na)) {
            tmp_trgt <- strsplit(na[i], "_")[[1]][1]  #gives IDs
            trgt <- strsplit(tmp_trgt, ",")[[1]]  #wants to be numeric
            rows <- match(trgt, probe[, "ID"])  #matches the row in the probe df
            score <- as.numeric(strsplit(na[i], "_")[[1]][2])  #gives the slope
            if (score > thrsh) {
                # this threshold can be chosen in the function and is set to
                # 0.001 by default, meaning the slope needs to be lower than
                # -0.001
                probe[rows, "flag"] <- paste0(probe[rows, "flag"], "PDD_")
                #all PDD candidates are flagged with 'PDD'
            }
        }
    }
    # the probe based df is returned as result
    probe
}
