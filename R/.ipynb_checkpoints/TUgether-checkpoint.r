# =========================================================================
# TUgether               Combines delay fragments into TUs
# -------------------------------------------------------------------------
#'
#'
#' TUgether combines delay fragments into TUs. The column "TU" is added.
#'
#' TUgether combines delay fragments into TUs. It uses score fun_increasing
#' on the start and end points of delay_fragments.

#' The function used is:

#' .score_fun_increasing

#' The input is the SummarizedExperiment object.

#' pen is the penalty for new fragments in the dynamic programming. Since high
#' scores are aimed, pen is negative.
#'
#' @param inp SummarizedExperiment: the input data frame with correct format.
#' @param cores cores: integer: the number of assigned cores for the task.
#' @param pen numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer fragments. Default -0.75.
#'
#' @return the SummarizedExperiment with the columns regarding the TU:
#' \describe{
#'   \item{ID:}{The bin/inp specific ID}
#'   \item{position:}{The bin/inp specific position}
#'   \item{position_segment:}{The position based segment}
#'   \item{delay_fragment:}{The delay fragment the bin belongs to}
#'   \item{intercept:}{The vintercept of fit through the respective delay
#'   fragment}
#'   \item{slope:}{The slope of the fit through the respective delay fragment}
#'   \item{TU:}{The overarching transcription unit}
#' }
#'
#' @examples
#' data(fragmentation_minimal)
#' TUgether(inp = fragmentation_minimal, cores = 2, pen = -0.75)
#' 
#' @export

TUgether <- function(inp,
                     cores = 1,
                     pen = -0.75) {
  # I.Preparations: the dataframe is configured and some other variables are
  # assigned
  registerDoMC(cores) # cores for DoMC
  
  rowRanges(inp)$TU <- NA
  
  # the dataframe is sorted by strand and position.
  inp <- inp_order(inp)
  #make the tmp_df
  tmp_df <- inp_df(inp, "ID", "position", "slope", "intercept",
                   "position_segment", "delay_fragment")
  
  #revert the order in plus
  tmp_df <- tmp_df_rev(tmp_df, "-")

  tmp_df <- tmp_df[!grepl("_NA|_T", tmp_df$delay_fragment), ]

  # _NA and _T fragments are included here, but are later discarded due to the
  # lack of velocity and intercept
  tmp_df[, "delay_fragment"] <- gsub("_O", "", tmp_df$delay_fragment)
  
  tmp_df <- na.omit(tmp_df)

  unique_seg <- unlist(unique(tmp_df$position_segment))

  count <- 1

  # II. Dynamic Programming: the scoring function is interpreted

  frags <- foreach(k = seq_along(unique_seg)) %dopar% {
    section <- tmp_df[which(tmp_df$position_segment == unique_seg[k]), ]

    # here we select the different delay fragments in addition to the...
    # ...position segment in the step before
    unique_frg <- unlist(unique(section$delay_fragment))

    # we check if we have at leat two delay fragments...*
    if (length(unique_frg) > 1) {
      # here we select all the slopes
      m <- as.numeric(section$slope[match(unique_frg, section$delay_fragment)])

      # and all the intercepts
      t <- as.numeric(section$intercept[match(unique_frg, section$delay_fragment)])

      # vectors for the start and end points are initiated
      start_v <- rep(NA, length(unique_frg))
      end_v <- rep(NA, length(unique_frg))

      # in this loop we calculate the start and end points
      for (i in seq_along(unique_frg)) {
        # start position (first position in the fragment)
        s <- section$position[section$delay_fragment == unique_frg[i]][1]
        # end position (last position in the fragment)
        e <- section$position[
          section$delay_fragment == unique_frg[i]][length(section$position[
            section$delay_fragment == unique_frg[i]])]
        # start_v is the calculated delay at the first point of each fragment
        start_v[i] <- m[i] * s + t[i]
        # end_v is the calculated delay at the last point of each fragment
        end_v[i] <- m[i] * e + t[i]
      }

      # put is a matrix, with all end points but the last one...
      # ...and all start points but the first one. We cant compare the first...
      # ...start and the last end to anything.
      put <- rbind(end_v[-length(end_v)], start_v[-1])

      # best frags (the score) is initiated with the penalty for a new fragment
      # because the scoring function now compares two fragments and does not
      # look at single points. We initiate with a single first fragment.
      best_frags <- c(pen)
      # best frags is initiated with the name of the single first
      # delay_fragment
      namer <- unique_frg[1]
      best_names <- c(namer)

      for (i in seq_len(ncol(put))) {
        # in the very first loop, already two fragments are scores, that's why
        # we initiated with a single fragment
        tmp_score <-
          score_fun_increasing(put[, seq_len(i)], unique_frg[seq_len(i + 1)])
        tmp_name <- names(tmp_score)
        # this checks if we have at least three fragments, despite being i>1!
        if (i > 1) {
          for (j in i:2) {
            # here we look for all possible combinations, but we cannot look
            # at the case where all fragments are on their own.
            tmp <-
              score_fun_increasing(put[, j:i], unique_frg[j:(i + 1)]) + pen +
              best_frags[j - 1]
            tmp_score <- c(tmp_score, tmp)
            tmp_n <- paste0(best_names[j - 1], "|", names(tmp))
            tmp_name <- c(tmp_name, tmp_n)
          }
        }
        # here we now look at the score if the last fragment that we are looking
        # at right now is alone.
        tmp_score <- c(tmp_score, pen + best_frags[i])
        namer <- unique_frg[i + 1]
        tmp_name <- c(tmp_name, paste0(best_names[i], "|", namer))
        # here we then decide for the highest score
        pos <- which(tmp_score == max(tmp_score))[1]
        tmp_score <- tmp_score[pos]
        tmp_name <- tmp_name[pos]
        best_frags <- c(best_frags, tmp_score)
        best_names <- c(best_names, tmp_name)
      }
      #* ...otherwise they are one TU by default
    } else {
      best_names <- unique_frg
    }
    
    best_names[length(best_names)]
  }

  # III. Fill the dataframe

  # the filling of the df here soes not care about outliers, as they don't exist
  # for TUs

  for (k in seq_along(frags)) {
    na <- strsplit(frags[[k]], "\\|")[[1]]
    for (i in seq_along(na)) {
      trgt <- strsplit(na[i], ",")[[1]]
      for (j in seq_along(trgt)) {
        # rows are selected by full delay fragments here
        rows <-
          which(trgt[j] == gsub("_O", "", rowRanges(inp)$delay_fragment))
        nam <- paste0("TU_", count)
        rowRanges(inp)$TU[rows] <- nam
      }
      count <- count + 1
    }
  }

  # NAs are treated the same way as in e.g. fragment_delay, so for example
  # terminal outlier fragments within one TU are combined with it

  if (sum(cumprod(is.na(rowRanges(inp)$TU))) > 0) {
    rowRanges(inp)$TU[seq_len(sum(cumprod(is.na(rowRanges(inp)$TU))))] <- "TU_0.5_NA"
  }
  row_NA <- which(is.na(rowRanges(inp)$TU))
  if (length(row_NA) > 0) {
    group <- c(row_NA[1])
    for (i in seq_along(row_NA)) {
      if (is.na(rowRanges(inp)$TU[row_NA[i] + 1])) {
        group <- c(group, row_NA[i] + 1)
      } else if (rowRanges(inp)$TU[group[1] - 1] == rowRanges(inp)$TU[group[length(group)] +
        1]) {
        rowRanges(inp)$TU[group] <- paste0(rowRanges(inp)$TU[group[1] - 1], "_NA")
        group <- row_NA[i + 1]
      } else if (rowRanges(inp)$TU[group[1] - 1] != rowRanges(inp)$TU[group[length(group)] +
        1]) {
        rowRanges(inp)$TU[group] <-
          paste0("TU_", as.numeric(gsub("TU_", "", rowRanges(inp)$TU[group[1] - 1])) +
                   0.5, "_NA")
        group <- row_NA[i + 1]
      }
    }
    rowRanges(inp)$TU[is.na(rowRanges(inp)$TU)] <-
      paste0("TU_", as.numeric(gsub(
        "TU_", "", rowRanges(inp)$TU[!is.na(rowRanges(inp)$TU)]
        [length(rowRanges(inp)$TU[!is.na(rowRanges(inp)$TU)])])) + 0.5, "_NA")
  }
  inp
}
