#' fragment_HL: performs the half_life fragmentation
#' fragment_HL makes HL_fragments based on delay_fragments and assigns all
#' gathered information to the probe based data frame.
#' The columns "HL_fragment" and "HL_mean_fragment" are added.
#' fragment_HL makes half-life_fragments and assigns the mean of each fragment.
#' The function used is:
#' .score_fun_ave.
#' The input is the probe, a dataframe with ID, half_life and delay_fragment.
#' pen is the penalty for new fragments in the dynamic programming, pen_out is
#' the outlier penalty.
#'
#' @param probe data frame: the probe based data frame.
#' @param cores integer: the number of assigned cores for the task.
#' @param pen numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer fragments. Default is the auto generated value.
#' @param pen_out numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer allowed outliers. Default is the auto
#' generated value.
#' 
#' @return the probe data frame with the columns regarding the half_life:
#' HL_fragment and HL_mean_fragment:
#' \describe{
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{position:}{The bin/probe specific position}
#'   \item{strand:}{The bin/probe specific strand}
#'   \item{half_life:}{The half-life of the bin/probe}
#'   \item{delay_fragment:}{The delay fragment the bin belongs to}
#'   \item{HL_fragment:}{The half-life fragment the bin belongs to}
#'   \item{HL_mean_fragment:}{The mean half-life value of the respective
#'   half-life fragment}
#' }
#'
#' @examples
#' data(fragmentation_minimal)
#' data(penalties_minimal)
#' fragment_HL(
#'   probe = fragmentation_minimal, cores = 2,
#'   pen = penalties_minimal["half_life_penalty"],
#'   pen_out = penalties_minimal["half_life_outlier_penalty"]
#' )
#' 
#' @export

fragment_HL <- function(probe, cores = 1, pen, pen_out) {
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
  req_cols_probe <- c("ID", "delay_fragment", "half_life")
  assert(
    all(req_cols_probe %in% colnames(probe)),
    paste0("'", req_cols_probe[which(!req_cols_probe %in% colnames(probe))],
           "' must be a column in 'probe'!")
  )

  # I.Preperations: the dataframe is configured and some other variables are
  # assigned

  registerDoMC(cores)

  probe <- probe[with(probe, order(-xtfrm(probe$strand), probe$position)), ]

  probe[probe$strand == "-", ] <-
    probe[probe$strand == "-", ][order(probe[probe$strand == "-", ]$position,
                                       decreasing = TRUE), ]

  probe[, "HL_fragment"] <- NA
  probe[, "HL_mean_fragment"] <- NA

  # the position is not relevant in this case
  tmp_df <-
    data.frame(
      ID = probe$ID,
      val = probe$half_life,
      seg = probe$delay_fragment
    )

  # the fragmentation is performed on the delay_fragments, independend on
  # if they are (terminal) outlires or NAs.
  tmp_df[, "seg"] <- gsub("_O|_NA", "", tmp_df$seg)

  # although _NA will most probably will be dismissed here, because there
  # should be no half-life, if there is no delay
  tmp_df <- na.omit(tmp_df)

  unique_seg <- unlist(unique(tmp_df$seg))

  count <- 1

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
          # fragments of the size at lest 4 are allowed for half-life
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
        nam <- paste0("Dc_", count, "_O")
        probe[rows, "HL_fragment"] <- nam
        # the value assigned here is the mean (compare velocity and intercept
        # for delay)
        probe[rows, "HL_mean_fragment"] <-
          as.numeric(strsplit(na[i], "_")[[1]][2])
      }
      rows <- match(trgt, probe[, "ID"])
      nam <- paste0("Dc_", count)
      probe[rows, "HL_fragment"] <- nam
      probe[rows, "HL_mean_fragment"] <-
        as.numeric(strsplit(na[i], "_")[[1]][2])
      count <- count + 1
    }
  }
  row_NA <- which(is.na(probe$HL_fragment))
  if (length(row_NA) > 0) {
    group <- c(row_NA[1])
    for (i in seq_along(row_NA)) {
      if (is.na(probe$HL_fragment[row_NA[i] + 1])) {
        group <- c(group, row_NA[i] + 1)
      } else if (gsub("_O", "", probe$HL_fragment[group[1] - 1]) ==
                 gsub("_O", "", probe$HL_fragment[group[length(group)] + 1])) {
        probe$HL_fragment[group] <-
          paste0(gsub("_O", "", probe$HL_fragment[group[1] - 1]), "_NA")
        probe$HL_mean_fragment[group] <-
          probe$HL_mean_fragment[group[1] - 1]
        group <- row_NA[i + 1]
      } else if (probe$HL_fragment[group[1] - 1] != probe$HL_fragment[
        group[length(group)] + 1]) {
        probe$HL_fragment[group] <-
          paste0("Dc_", as.numeric(
            gsub("Dc_|_O", "", probe$HL_fragment[group[1] - 1])) + 0.5, "_NA")
        group <- row_NA[i + 1]
      }
    }
  }
  probe$HL_fragment[is.na(probe$HL_fragment)] <- paste0("Dc_", as.numeric(
      gsub("Dc_|_O", "", probe$HL_fragment[!is.na(probe$HL_fragment)]
           [length(probe$HL_fragment[!is.na(probe$HL_fragment)])])) + .5, "_NA")
  probe
}
