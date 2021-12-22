#' TUgether: combines delay fragments into TUs.
#'
#' TUgether combines delay fragments into TUs.
#' The column "TU" is added.
#'
#' TUgether combines delay fragments into TUs. It uses score fun_increasing
#' on the start and end points of delay_fragments.
#' The function used is:
#' .score_fun_increasing
#' The input is the probe, a dataframe with ID, velocity_fragment, intercept,
#' position, strand, position_segment, and delay_fragment.
#' pen is the penalty for new fragments in the dynamic programming. Since high
#' scores are aimed, Pen is negative.
#'
#' @param probe data frame: the probe based data frame.
#' @param cores cores: integer: the number of assigned cores for the task.
#' @param pen numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer fragments. Default -0.75.
#'
#' @return the probe data frame with the columns regarding the TU:
#' \describe{
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{position:}{The bin/probe specific position}
#'   \item{strand:}{The bin/probe specific strand}
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
#' TUgether(probe = fragmentation_minimal, cores = 2, pen = -0.75)
#' 
#' @export

TUgether <- function(probe,
                     cores = 1,
                     pen = -0.75) {
  num_args <- list(pen)
  names(num_args) <- c("pen")
  assert(
    all(unlist(lapply(
      num_args,
      FUN = function(x) {
        (is.numeric(x) &
          length(x) == 1)
      }
    ))),
    paste0(
      "'",
      names(which(unlist(
        lapply(
          num_args,
          FUN = function(x) {
            (is.numeric(x) &
              length(x) == 1)
          }
        )
      ) == FALSE))[1],
      "' must be numeric of length one or given by the logs"
    )
  )
  assert(cores > 0, "'cores' must be a positive integer")
  req_cols_probe <-
    c(
      "ID",
      "slope",
      "intercept",
      "position",
      "strand",
      "position_segment",
      "delay_fragment"
    )
  assert(
    all(req_cols_probe %in% colnames(probe)),
    paste0("'", req_cols_probe[which(!req_cols_probe %in% colnames(probe))],
           "' must be a column in 'probe'!")
  )

  # I.Preparations: the dataframe is configured and some other variables are
  # assigned

  registerDoMC(cores)

  probe <- probe[with(probe, order(-xtfrm(probe$strand), probe$position)), ]
  probe[probe$strand == "-", ] <- probe[probe$strand == "-", ][
    order(probe[probe$strand == "-", ]$position, decreasing = TRUE), ]

  probe[, "TU"] <- NA

  tmp_df <-
    data.frame(
      ID = probe$ID,
      m = probe$slope,
      t = probe$intercept,
      position = probe$position,
      strand = probe$strand,
      seg = probe$position_segment,
      frg = probe$delay_fragment
    )

  tmp_df <- tmp_df[!grepl("_NA|_T", tmp_df$frg), ]

  # _NA and _T fragments are included here, but are later discarded due to the
  # lack of velocity and intercept
  tmp_df[, "frg"] <- gsub("_O", "", tmp_df$frg)

  # like in fragment delay
  tmp_df[tmp_df$strand == "-", "position"] <-
    ((tmp_df[tmp_df$strand == "-", "position"]) -
      (tmp_df[tmp_df$strand == "-", ][1, "position"])) * -1

  tmp_df <- na.omit(tmp_df)

  unique_seg <- unlist(unique(tmp_df$seg))

  count <- 1

  # II. Dynamic Programming: the scoring function is interpreted

  frags <- foreach(k = seq_along(unique_seg)) %dopar% {
    section <- tmp_df[which(tmp_df$seg == unique_seg[k]), ]

    # here we select the different delay fragments in addition to the...
    # ...position segment in the step before
    unique_frg <- unlist(unique(section$frg))

    # we check if we have at leat two delay fragments...*
    if (length(unique_frg) > 1) {
      # here we select all the slopes
      m <- as.numeric(section$m[match(unique_frg, section$frg)])

      # and all the intercepts
      t <- as.numeric(section$t[match(unique_frg, section$frg)])

      # vectors for the start and end points are initiated
      start_v <- rep(NA, length(unique_frg))
      end_v <- rep(NA, length(unique_frg))

      # in this loop we calculate the start and end points
      for (i in seq_along(unique_frg)) {
        # start position (first position in the fragment)
        s <- section$position[section$frg == unique_frg[i]][1]
        # end position (last position in the fragment)
        e <- section$position[
          section$frg == unique_frg[i]][length(section$position[
            section$frg == unique_frg[i]])]
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
          which(trgt[j] == gsub("_O", "", probe[, "delay_fragment"]))
        nam <- paste0("TU_", count)
        probe[rows, "TU"] <- nam
      }
      count <- count + 1
    }
  }

  # NAs are treated the same way as in e.g. fragment_delay, so for example
  # terminal outlier fragments within one TU are combined with it

  if (sum(cumprod(is.na(probe$TU))) > 0) {
    probe$TU[seq_len(sum(cumprod(is.na(probe$TU))))] <- "TU_0.5_NA"
  }
  row_NA <- which(is.na(probe$TU))
  if (length(row_NA) > 0) {
    group <- c(row_NA[1])
    for (i in seq_along(row_NA)) {
      if (is.na(probe$TU[row_NA[i] + 1])) {
        group <- c(group, row_NA[i] + 1)
      } else if (probe$TU[group[1] - 1] == probe$TU[group[length(group)] +
        1]) {
        probe$TU[group] <- paste0(probe$TU[group[1] - 1], "_NA")
        group <- row_NA[i + 1]
      } else if (probe$TU[group[1] - 1] != probe$TU[group[length(group)] +
        1]) {
        probe$TU[group] <-
          paste0("TU_", as.numeric(gsub("TU_", "", probe$TU[group[1] - 1])) +
                   0.5, "_NA")
        group <- row_NA[i + 1]
      }
    }
    probe$TU[is.na(probe$TU)] <-
      paste0("TU_", as.numeric(gsub(
        "TU_", "", probe$TU[!is.na(probe$TU)]
        [length(probe$TU[!is.na(probe$TU)])])) + 0.5, "_NA")
  }
  probe
}
