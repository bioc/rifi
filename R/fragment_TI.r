#' fragment_TI: performs the TI fragmentation.
#' fragment_TI makes TI_fragments based on TUs and assigns all gathered
#' information to the SummarizedExperiment object.
#' The columns "TI_termination_fragment" and the TI_mean_termination_factor
#' are added.
#' The function used is:
#'  .score_fun_ave.
#' The input is the SummarizedExperiment object.
#' pen is the penalty for new fragments in the dynamic programming, pen_out is
#' the outlier penalty.
#'
#' @param inp SummarizedExperiment: the input data frame with correct format.
#' @param cores cores: integer: the number of assigned cores for the task.
#' @param pen numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer fragments. Default is the auto generated value.
#' @param pen_out numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer allowed outliers. Default is the auto generated
#' value.
#' 
#' @return the SummarizedExperiment object: with TI_termination_fragment and
#' TI_termination_mean_fragment added to the rowRanges.
#' 
#' @examples
#' data(fragmentation_minimal)
#' fragment_TI(inp = fragmentation_minimal, cores = 2, pen = 2, pen_out = 1)
#' 
#' @export

fragment_TI <- function(inp, cores = 1, pen, pen_out) {
  
  # I.Preparations: the dataframe is configured and some other variables are
  # assigned
  registerDoMC(cores) # cores for DoMC
  
  rowRanges(inp)$TI_termination_fragment <- NA
  rowRanges(inp)$TI_mean_termination_factor <- NA
  
  # the dataframe is sorted by strand and position.
  inp <- inp_order(inp)
  #make the tmp_df
  tmp_df <- inp_df(inp, "ID", "TI_termination_factor", "TU", "flag")
  #revert the order in plus
  tmp_df <- tmp_df_rev(tmp_df, "-")
  
  # this is the same way of selecting the IDs as the selection for the TI_fit
  corr_IDs <- tmp_df$ID[grep("TI", tmp_df$flag)]

  tmp_df <- tmp_df[tmp_df$ID %in% corr_IDs, ]

  tmp_df <- na.omit(tmp_df)

  # only the TUs that are flagged with TI at any position are considered
  unique_seg <- unlist(unique(tmp_df$TU))

  # only true TUs are taken (no TU_NA)
  unique_seg <- unique_seg[grep("_NA", unique_seg, invert = TRUE)]

  count <- 1

  if (length(unique_seg) == 0) {
    return(inp)
  }
  # II. Dynamic Programming: the scoring function is interpreted

  frags <- foreach(k = seq_along(unique_seg)) %dopar% {
    section <- tmp_df[which(tmp_df$TU == unique_seg[k]), ]

    best_frags <- c()
    best_names <- c()
    if (nrow(section) > 1) {
      # only segments with more than one value are grouped...*

    for (i in 2:nrow(section)) {
      tmp_score <-
        score_fun_ave(section[seq_len(i), "TI_termination_factor"],
                      section[seq_len(i), "ID"], pen_out)
      tmp_name <- names(tmp_score)
      if (i > 3) {
        for (j in (i - 1):3) {
          tmp_val <- section[j:i, "TI_termination_factor"]
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
        score_fun_ave(section[, "TI_termination_factor"], section[, "ID"], pen_out)
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
        nam <- paste0("TI_", count, "_O")
        rowRanges(inp)$TI_termination_fragment[rows] <- nam
        rowRanges(inp)$TI_mean_termination_factor[rows] <-
          as.numeric(strsplit(na[i], "_")[[1]][2])
      }
      rows <- match(trgt, rowRanges(inp)$ID)
      nam <- paste0("TI_", count)
      rowRanges(inp)$TI_termination_fragment[rows] <- nam
      rowRanges(inp)$TI_mean_termination_factor[rows] <-
        as.numeric(strsplit(na[i], "_")[[1]][2])
      count <- count + 1
    }
  }

  if (sum(cumprod(is.na(rowRanges(inp)$TI_termination_fragment))) > 0) {
    # this is a trick to fill the first position
    rowRanges(inp)$TI_termination_fragment[
      seq_len(sum(cumprod(is.na(rowRanges(inp)$TI_termination_fragment))))] <- "bla"
  }

  row_NA <- which(is.na(rowRanges(inp)$TI_termination_fragment))
  if (length(row_NA) > 0) {
    group <- c(row_NA[1])
    for (i in seq_along(row_NA)) {
      if (is.na(rowRanges(inp)$TI_termination_fragment[row_NA[i] + 1])) {
        group <- c(group, row_NA[i] + 1)
      }
      # this time only the ones within a fragment are considered
      else if (gsub("_O", "", rowRanges(inp)$TI_termination_fragment[group[1] - 1]) ==
        gsub("_O", "",
             rowRanges(inp)$TI_termination_fragment[group[length(group)] + 1])) {
        rowRanges(inp)$TI_termination_fragment[group] <-
          paste0(gsub("_O", "",
                      rowRanges(inp)$TI_termination_fragment[group[1] - 1]), "_NA")
        rowRanges(inp)$TI_mean_termination_factor[group] <-
          rowRanges(inp)$TI_mean_termination_factor[group[1] - 1]
        group <- row_NA[i + 1]
      }
      # the ones between are still NA
      else if (rowRanges(inp)$TI_termination_fragment[group[1] - 1] !=
               rowRanges(inp)$TI_termination_fragment[group[length(group)] + 1]) {
        group <- row_NA[i + 1]
      }
    }
  }
  # here the first position is reverted to NA
  rowRanges(inp)$TI_termination_fragment[which(rowRanges(inp)$TI_termination_fragment ==
                                        "bla")] <- NA
  inp
}
