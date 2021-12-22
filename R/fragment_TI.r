#' fragment_TI: performs the TI fragmentation.
#' fragment_TI makes TI_fragments based on TUs and assigns all gathered
#' information to the probe based data frame.
#' The columns "TI_termination_fragment" and the TI_mean_termination_factor
#' are added.
#' The function used is:
#'  .score_fun_ave.
#' The input is the probe, a dataframe with ID, TI_termination_factor and TU.
#' pen is the penalty for new fragments in the dynamic programming, pen_out is
#' the outlier penalty.
#'
#' @param probe data frame: the probe based data frame.
#' @param cores cores: integer: the number of assigned cores for the task.
#' @param pen numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer fragments. Default is the auto generated value.
#' @param pen_out numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer allowed outliers. Default is the auto generated
#' value.
#' 
#' @return the probe data frame with the columns regarding the TI:
#' TI_termination_fragment and TI_mean_termination_fragment:
#' \describe{
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{position:}{The bin/probe specific position}
#'   \item{strand:}{The bin/probe specific strand}
#'   \item{flag:}{Information on which fitting model is applied}
#'   \item{TI_termination_factor:}{The termination factor of the bin/probe}
#'   \item{TU:}{The overarching transcription unit}
#'   \item{TI_termination_fragment:}{The TI fragment the bin belongs to}
#'   \item{TI_mean_termination_factor:}{The mean termination factor of the
#'   respective TI fragment}
#' }
#' 
#' @examples
#' data(fragmentation_minimal)
#' data(penalties_minimal)
#' fragment_TI(
#'   probe = fragmentation_minimal, cores = 2,
#'   pen = penalties_minimal["TI_penalty"],
#'   pen_out = penalties_minimal["TI_outlier_penalty"]
#' )
#' 
#' @export

fragment_TI <- function(probe, cores = 1, pen, pen_out) {
  num_args <- list(pen, pen_out)
  names(num_args) <- c("pen", "pen_out")
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
  req_cols_probe <- c("ID", "TU", "TI_termination_factor", "flag")
  assert(
    all(req_cols_probe %in% colnames(probe)),
    paste0("'", req_cols_probe[which(!req_cols_probe %in% colnames(probe))],
           "' must be a column in 'probe'!")
  )

  # I.Preperations: the dataframe is configured and some other variables are
  # assigned

  registerDoMC(cores)

  probe <- probe[with(probe, order(-xtfrm(probe$strand), probe$position)), ]
  probe[probe$strand == "-", ] <- probe[probe$strand == "-", ][
      order(probe[probe$strand == "-", ]$position, decreasing = TRUE), ]

  probe[, "TI_termination_fragment"] <- NA
  probe[, "TI_mean_termination_factor"] <- NA

  # We additionally need the flag information here
  tmp_df <-
    data.frame(
      ID = probe$ID,
      val = probe$TI_termination_factor,
      seg = probe$TU,
      flag = probe$flag
    )

  # this is the same way of selecting the IDs as the selection for the TI_fit
  corr_IDs <- tmp_df$ID[grep("TI", tmp_df$flag)]

  tmp_df <- tmp_df[tmp_df$ID %in% corr_IDs, ]

  tmp_df <- na.omit(tmp_df)

  # only the TUs that are flagged with TI at any position are considered
  unique_seg <- unlist(unique(tmp_df$seg))

  # only true TUs are taken (no TU_NA)
  unique_seg <- unique_seg[grep("_NA", unique_seg, invert = TRUE)]

  count <- 1

  if (length(unique_seg) == 0) {
    return(probe)
  }
  # II. Dynamic Programming: the scoring function is interpreted

  frags <- foreach(k = seq_along(unique_seg)) %dopar% {
    section <- tmp_df[which(tmp_df$seg == unique_seg[k]), ]

    best_frags <- c()
    best_names <- c()
    if (nrow(section) > 1) {
      # only segments with more than one value are grouped...*

    for (i in 2:nrow(section)) {
      tmp_score <-
        score_fun_ave(section[seq_len(i), "val"],
                      section[seq_len(i), "ID"], pen_out)
      tmp_name <- names(tmp_score)
      if (i > 3) {
        for (j in (i - 1):3) {
          tmp_val <- section[j:i, "val"]
          tmp_ID <- section[j:i, "ID"]
          tmp <-
            score_fun_ave(tmp_val, tmp_ID, pen_out) + pen + best_frags[j - 2]
          tmp_score <- c(tmp_score, tmp)
          tmp_n <- paste0(best_names[j - 2], "|", names(tmp))
          tmp_name <- c(tmp_name, tmp_n)
        }
      }
      pos <- which(tmp_score == min(tmp_score))[1]
      tmp_score <- tmp_score[pos]
      tmp_name <- tmp_name[pos]
      best_frags <- c(best_frags, tmp_score)
      best_names <- c(best_names, tmp_name)
    }
} else {
      #* ...all segments with less than two values are grouped automatically
      tmp_score <-
        score_fun_ave(section[, "val"], section[, "ID"], pen_out)
      tmp_name <- names(tmp_score)
      best_names <- c(best_names, tmp_name)
    }
    best_names[length(best_names)]
  }

  # III. Fill the dataframe

  for (k in seq_along(frags)) {
    na <- strsplit(frags[[k]], "\\|")[[1]]

    for (i in seq_along(na)) {
      tmp_trgt <- strsplit(na[i], "_")[[1]][1]
      trgt <- strsplit(tmp_trgt, ",")[[1]]
      if (length(strsplit(na[i], "_")[[1]]) == 3) {
        tmp_outl <- strsplit(na[i], "_")[[1]][3]
        outl <- strsplit(tmp_outl, ",")[[1]]
        trgt <- trgt[-which(trgt %in% outl)]
        rows <- match(outl, probe[, "ID"])
        nam <- paste0("TI_", count, "_O")
        probe[rows, "TI_termination_fragment"] <- nam
        probe[rows, "TI_mean_termination_factor"] <-
          as.numeric(strsplit(na[i], "_")[[1]][2])
      }
      rows <- match(trgt, probe[, "ID"])
      nam <- paste0("TI_", count)
      probe[rows, "TI_termination_fragment"] <- nam
      probe[rows, "TI_mean_termination_factor"] <-
        as.numeric(strsplit(na[i], "_")[[1]][2])
      count <- count + 1
    }
  }

  if (sum(cumprod(is.na(probe$TI_termination_fragment))) > 0) {
    # this is a trick to fill the first position
    probe$TI_termination_fragment[
      seq_len(sum(cumprod(is.na(probe$TI_termination_fragment))))] <- "bla"
  }

  row_NA <- which(is.na(probe$TI_termination_fragment))
  if (length(row_NA) > 0) {
    group <- c(row_NA[1])
    for (i in seq_along(row_NA)) {
      if (is.na(probe$TI_termination_fragment[row_NA[i] + 1])) {
        group <- c(group, row_NA[i] + 1)
      }
      # this time only the ones within a fragment are considered
      else if (gsub("_O", "", probe$TI_termination_fragment[group[1] - 1]) ==
        gsub("_O", "",
             probe$TI_termination_fragment[group[length(group)] + 1])) {
        probe$TI_termination_fragment[group] <-
          paste0(gsub("_O", "",
                      probe$TI_termination_fragment[group[1] - 1]), "_NA")
        probe$TI_mean_termination_factor[group] <-
          probe$TI_mean_termination_factor[group[1] - 1]
        group <- row_NA[i + 1]
      }
      # the ones between are still NA
      else if (probe$TI_termination_fragment[group[1] - 1] !=
               probe$TI_termination_fragment[group[length(group)] + 1]) {
        group <- row_NA[i + 1]
      }
    }
  }
  # here the first position is reverted to NA
  probe$TI_termination_fragment[which(probe$TI_termination_fragment ==
                                        "bla")] <- NA
  probe
}
