# =========================================================================
# fragment_inty            Performs the intensity fragmentation
# -------------------------------------------------------------------------
#'
#' fragment_inty makes intensity_fragments based on HL_fragments and assigns
#' all gathered information to the SummarizedExperiment object.

#' The columns "intensity_fragment" and "intensity_mean_fragment" are added.

#' fragment_inty makes intensity_fragments and assigns the mean of each
#' fragment.

#' The function used is:

#' .score_fun_ave.

#' The input is the the SummarizedExperiment object.

#' pen is the penalty for new fragments in the dynamic programming, pen_out
#' is the outlier penalty.
#'
#' @param inp SummarizedExperiment: the input data frame with correct format.
#' @param cores cores: integer: the number of assigned cores for the task.
#' @param pen numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer fragments. Default is the auto generated value.
#' @param pen_out numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer allowed outliers. Default is the auto generated
#' value.
#'
#' @return the SummarizedExperiment object: with intensity_fragment and
#' intensity_mean_fragment added to the rowRanges.
#'
#' @examples
#' data(fragmentation_minimal)
#' fragment_inty(inp = fragmentation_minimal, cores = 2, pen = 2, pen_out = 1)
#' 
#' @export

fragment_inty <- function(inp, cores = 1, pen, pen_out) {
  
  # I.Preparations: the dataframe is configured and some other variables are
  # assigned
  registerDoMC(cores) # cores for DoMC
  
  rowRanges(inp)$intensity_fragment <- NA
  rowRanges(inp)$intensity_mean_fragment <- NA
  
  # the dataframe is sorted by strand and position.
  inp <- inp_order(inp)
  #make the tmp_df
  tmp_df <- inp_df(inp, "ID", "intensity", "HL_fragment")
  #revert the order in plus
  tmp_df <- tmp_df_rev(tmp_df, "-")

  # here it is important that (terminal) outliers and NAs are considered, as all
  # bins have an intensity that needs to be grouped.
  tmp_df[, "HL_fragment"] <- gsub("_O|_NA", "", tmp_df$HL_fragment)

  # this is just for safety, nothing should be omitted here
  tmp_df <- na.omit(tmp_df)

  tmp_df$intensity <- log2(tmp_df$intensity)

  unique_seg <- unlist(unique(tmp_df$HL_fragment))

  count <- 1

  # II. Dynamic Programming: the scoring function is interpreted

  frags <- foreach(k = seq_along(unique_seg)) %dopar% {
    section <- tmp_df[which(tmp_df$HL_fragment == unique_seg[k]), ]

    best_frags <- c()
    best_names <- c()

    if (nrow(section) > 1) {
      # only segments with more than one value are grouped...*

      for (i in 2:nrow(section)) {
        tmp_score <-
          score_fun_ave(section[seq_len(i), "intensity"],
                        section[seq_len(i), "ID"], pen_out)
        tmp_name <- names(tmp_score)
        if (i > 3) {
          for (j in (i - 1):3) {
            tmp_val <- section[j:i, "intensity"]
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
        score_fun_ave(section[, "intensity"], section[, "ID"], pen_out)
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
        rows <- match(outl, rowRanges(inp)$ID)
        nam <- paste0("I_", count, "_O")
        rowRanges(inp)$intensity_fragment[rows] <- nam
        # the mean needs to be reverted by the correction factor here
        rowRanges(inp)$intensity_mean_fragment[rows] <- 2**
          (as.numeric(strsplit(na[i], "_")[[1]][2]))
      }
      rows <- match(trgt, rowRanges(inp)$ID)
      nam <- paste0("I_", count)
      rowRanges(inp)$intensity_fragment[rows] <- nam
      # the mean needs to be reverted by the correction factor here
      rowRanges(inp)$intensity_mean_fragment[rows] <- 2**
        (as.numeric(strsplit(na[i], "_")[[1]][2]))
      count <- count + 1
    }
  }

  if (sum(cumprod(is.na(rowRanges(inp)$intensity_fragment))) > 0) {
    rowRanges(inp)$intensity_fragment[seq_len(
      sum(cumprod(is.na(rowRanges(inp)$intensity_fragment))))] <- "I_0.5_NA"
  }
  row_NA <- which(is.na(rowRanges(inp)$intensity_fragment))
  if (length(row_NA) > 0) {
    group <- c(row_NA[1])
    for (i in seq_along(row_NA)) {
      if (is.na(rowRanges(inp)$intensity_fragment[row_NA[i] + 1])) {
        group <- c(group, row_NA[i] + 1)
      } else if (gsub("_O", "", rowRanges(inp)$intensity_fragment[group[1] - 1]) ==
        gsub("_O", "", rowRanges(inp)$intensity_fragment[group[length(group)] + 1])) {
        rowRanges(inp)$intensity_fragment[group] <-
          paste0(gsub("_O", "", rowRanges(inp)$intensity_fragment[group[1] - 1]), "_NA")
        rowRanges(inp)$intensity_mean_fragment[group] <-
          rowRanges(inp)$intensity_mean_fragment[group[1] - 1]
        group <- row_NA[i + 1]
      } else if (rowRanges(inp)$intensity_fragment[group[1] - 1] !=
                 rowRanges(inp)$intensity_fragment[group[length(group)] + 1]) {
        rowRanges(inp)$intensity_fragment[group] <-
          paste0("I_", as.numeric(gsub(
            "I_|_O", "", rowRanges(inp)$intensity_fragment[group[1] - 1]
          )) + 0.5, "_NA")
        group <- row_NA[i + 1]
      }
    }
    rowRanges(inp)$intensity_fragment[is.na(rowRanges(inp)$intensity_fragment)] <-
      paste0("I_", as.numeric(gsub(
        "I_|_O", "", rowRanges(inp)$intensity_fragment[!is.na(rowRanges(inp)$intensity_fragment)]
        [length(rowRanges(inp)$intensity_fragment[!is.na(rowRanges(inp)$intensity_fragment)])]
      )) + 0.5, "_NA")
  }
  inp
}
