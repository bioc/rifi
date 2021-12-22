#' fragment_inty: performs the intensity fragmentation
#' fragment_inty makes intensity_fragments based on HL_fragments and assigns
#' all gathered information to the probe based data frame.
#' The columns "intensity_fragment" and "intensity_mean_fragment" are added.
#' fragment_inty makes intensity_fragments and assigns the mean of each
#' fragment.
#' The function used is:
#' .score_fun_ave.
#' The input is the probe, a dataframe with ID, intensity and HL_fragment.
#' pen is the penalty for new fragments in the dynamic programming, pen_out
#' is the outlier penalty.

#' @param probe data frame: the probe based data frame.
#' @param cores cores: integer: the number of assigned cores for the task.
#' @param pen numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer fragments. Default is the auto generated value.
#' @param pen_out numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer allowed outliers. Default is the auto generated
#' value.
#'
#' @return the probe data frame with the columns regarding the intensity:
#' intensity_fragment and intensity_mean_fragment:
#' \describe{
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{position:}{The bin/probe specific position}
#'   \item{strand:}{The bin/probe specific strand}
#'   \item{intensity:}{The relative intensity at time point 0}
#'   \item{HL_fragment:}{The half-life fragment the bin belongs to}
#'   \item{intensity_fragment:}{The intensity fragment the bin belongs to}
#'   \item{intensity_mean_fragment:}{The mean intensity value of the respective
#'   intensity fragment}
#' }
#'
#' @examples
#' data(fragmentation_minimal)
#' data(penalties_minimal)
#' fragment_inty(
#'   probe = fragmentation_minimal, cores = 2,
#'   pen = penalties_minimal["intensity_penalty"],
#'   pen_out = penalties_minimal["intensity_outlier_penalty"]
#' )
#' 
#' @export

fragment_inty <- function(probe, cores = 1, pen, pen_out) {
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
  req_cols_probe <- c("ID", "HL_fragment", "intensity")
  assert(
    all(req_cols_probe %in% colnames(probe)),
    paste0("'", req_cols_probe[which(
      !req_cols_probe %in% colnames(probe))], "' must be a column in 'probe'!")
  )

  # I.Preparations: the dataframe is configured and some other variables
  # are assigned

  registerDoMC(cores)

  probe <- probe[with(probe, order(-xtfrm(probe$strand), probe$position)), ]

  probe[probe$strand == "-", ] <-
    probe[probe$strand == "-", ][order(probe[probe$strand == "-",
                                             ]$position, decreasing = TRUE), ]

  probe[, "intensity_fragment"] <- NA
  probe[, "intensity_mean_fragment"] <- NA

  # here the assigned penalty is cached to be changed dynamically in the loop
  tmp_df <-
    data.frame(
      ID = probe$ID,
      val = probe$intensity,
      seg = probe$HL_fragment
    )

  # here it is important that (terminal) outliers and NAs are considered, as all
  # bins have an intensity that needs to be grouped.
  tmp_df[, "seg"] <- gsub("_O|_NA", "", tmp_df$seg)

  # this is just for safety, nothing should be omitted here
  tmp_df <- na.omit(tmp_df)

  tmp_df$val <- log2(tmp_df$val)

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
        nam <- paste0("I_", count, "_O")
        probe[rows, "intensity_fragment"] <- nam
        # the mean needs to be reverted by the correction factor here
        probe[rows, "intensity_mean_fragment"] <- 2**
          (as.numeric(strsplit(na[i], "_")[[1]][2]))
      }
      rows <- match(trgt, probe[, "ID"])
      nam <- paste0("I_", count)
      probe[rows, "intensity_fragment"] <- nam
      # the mean needs to be reverted by the correction factor here
      probe[rows, "intensity_mean_fragment"] <- 2**
        (as.numeric(strsplit(na[i], "_")[[1]][2]))
      count <- count + 1
    }
  }

  if (sum(cumprod(is.na(probe$intensity_fragment))) > 0) {
    probe$intensity_fragment[seq_len(
      sum(cumprod(is.na(probe$intensity_fragment))))] <- "I_0.5_NA"
  }
  row_NA <- which(is.na(probe$intensity_fragment))
  if (length(row_NA) > 0) {
    group <- c(row_NA[1])
    for (i in seq_along(row_NA)) {
      if (is.na(probe$intensity_fragment[row_NA[i] + 1])) {
        group <- c(group, row_NA[i] + 1)
      } else if (gsub("_O", "", probe$intensity_fragment[group[1] - 1]) ==
        gsub("_O", "", probe$intensity_fragment[group[length(group)] + 1])) {
        probe$intensity_fragment[group] <-
          paste0(gsub("_O", "", probe$intensity_fragment[group[1] - 1]), "_NA")
        probe$intensity_mean_fragment[group] <-
          probe$intensity_mean_fragment[group[1] - 1]
        group <- row_NA[i + 1]
      } else if (probe$intensity_fragment[group[1] - 1] !=
                 probe$intensity_fragment[group[length(group)] + 1]) {
        probe$intensity_fragment[group] <-
          paste0("I_", as.numeric(gsub(
            "I_|_O", "", probe$intensity_fragment[group[1] - 1]
          )) + 0.5, "_NA")
        group <- row_NA[i + 1]
      }
    }
    probe$intensity_fragment[is.na(probe$intensity_fragment)] <-
      paste0("I_", as.numeric(gsub(
        "I_|_O", "", probe$intensity_fragment[!is.na(probe$intensity_fragment)]
        [length(probe$intensity_fragment[!is.na(probe$intensity_fragment)])]
      )) + 0.5, "_NA")
  }
  probe
}
