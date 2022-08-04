#' =========================================================================
#' fragment_HL             
#' -------------------------------------------------------------------------
#' fragment_HL performs the half_life fragmentation
#' 
#' fragment_HL makes HL_fragments based on delay_fragments and assigns all
#' gathered information to the SummarizedExperiment object.
#' 
#' The columns "HL_fragment" and "HL_mean_fragment" are added.
#' 
#' fragment_HL makes half-life_fragments and assigns the mean of each fragment.
#' 
#' The function used is:
#' 
#' .score_fun_ave.
#' 
#' The input the SummarizedExperiment object.
#' 
#' pen is the penalty for new fragments in the dynamic programming, pen_out is
#' the outlier penalty.
#'
#' @param inp SummarizedExperiment: the input data frame with correct format.
#' @param cores integer: the number of assigned cores for the task.
#' @param pen numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer fragments. Default is the auto generated value.
#' @param pen_out numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer allowed outliers. Default is the auto
#' generated value.
#' 
#' @return The SummarizedExperiment object: 
#'  \describe{
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{position:}{The bin/probe specific position}
#'   \item{intensity:}{The relative intensity at time point 0}
#'   \item{probe_TI:}{An internal value to determine which fitting model is
#'    applied}
#'   \item{flag:}{Information on which fitting model is applied}
#'   \item{position_segment:}{The position based segment}
#'   \item{delay:}{The delay value of the bin/probe}
#'   \item{half_life:}{The half-life of the bin/probe}
#'   \item{TI_termination_factor:}{String, the factor of TI fragment}
#'   \item{delay_fragment:}{The delay fragment the bin belongs to}
#'   \item{velocity_fragment:}{The velocity value of the respective delay
#'   fragment}
#'   \item{intercept:}{The vintercept of fit through the respective delay
#'   fragment}
#'   \item{slope:}{The slope of the fit through the respective delay fragment}
#'   \item{HL_fragment:}{The half-life fragment the bin belongs to}
#'   \item{HL_mean_fragment:}{The mean half-life value of the respective
#'   half-life fragment}
#'   }
#'
#' @examples
#' data(fragmentation_minimal)
#' fragment_HL(inp = fragmentation_minimal, cores = 2, pen = 2, pen_out = 1)
#' 
#' @export

fragment_HL <- function(inp, cores = 1, pen, pen_out) {
  
  # I.Preparations: the dataframe is configured and some other variables are
  # assigned
  registerDoMC(cores) # cores for DoMC
  
  rowRanges(inp)$HL_fragment <- NA
  rowRanges(inp)$HL_mean_fragment <- NA

  # the dataframe is sorted by strand and position.
  inp <- inp_order(inp)
  #make the tmp_df
  tmp_df <- inp_df(inp, "ID", "half_life", "delay_fragment")
  #revert the order in plus
  tmp_df <- tmp_df_rev(tmp_df, "-")
  
  # the fragmentation is performed on the delay_fragments, independend on
  # if they are (terminal) outlires or NAs.
  tmp_df[, "delay_fragment"] <- gsub("_O|_NA", "", tmp_df$delay_fragment)

  # although _NA will most probably will be dismissed here, because there
  # should be no half-life, if there is no delay
  tmp_df <- na.omit(tmp_df)

  unique_seg <- unlist(unique(tmp_df$delay_fragment))

  count <- 1

  # II. Dynamic Programming: the scoring function is interpreted

  frags <- foreach(k = seq_along(unique_seg)) %dopar% {
    section <- tmp_df[which(tmp_df$delay_fragment == unique_seg[k]), ]

    best_frags <- c()
    best_names <- c()

    if (nrow(section) > 1) {
      # only segments with more than one value are grouped...*

      for (i in 2:nrow(section)) {
        tmp_score <-
          score_fun_ave(section[seq_len(i), "half_life"],
                        section[seq_len(i), "ID"], pen_out)
        tmp_name <- names(tmp_score)
        if (i > 3) {
          # fragments of the size at lest 4 are allowed for half-life
          for (j in (i - 1):3) {
            tmp_val <- section[j:i, "half_life"]
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
        score_fun_ave(section[, "half_life"], section[, "ID"], pen_out)
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
        nam <- paste0("Dc_", count, "_O")
        rowRanges(inp)$HL_fragment[rows] <- nam
        # the value assigned here is the mean (compare velocity and intercept
        # for delay)
        rowRanges(inp)$HL_mean_fragment[rows] <-
          as.numeric(strsplit(na[i], "_")[[1]][2])
      }
      rows <- match(trgt, rowRanges(inp)$ID)
      nam <- paste0("Dc_", count)
      rowRanges(inp)$HL_fragment[rows] <- nam
      rowRanges(inp)$HL_mean_fragment[rows] <-
        as.numeric(strsplit(na[i], "_")[[1]][2])
      count <- count + 1
    }
  }
  if (sum(cumprod(is.na(rowRanges(inp)$HL_fragment))) > 0) {
    rowRanges(inp)$HL_fragment[seq_len(
      sum(cumprod(is.na(rowRanges(inp)$HL_fragment))))] <- "Dc_0.5_NA"
  }
  row_NA <- which(is.na(rowRanges(inp)$HL_fragment))
  if (length(row_NA) > 0) {
    group <- c(row_NA[1])
    for (i in seq_along(row_NA)) {
      if (is.na(rowRanges(inp)$HL_fragment[row_NA[i] + 1])) {
        group <- c(group, row_NA[i] + 1)
      } else if (gsub("_O", "", rowRanges(inp)$HL_fragment[group[1] - 1]) ==
                 gsub("_O", "", rowRanges(inp)$HL_fragment[group[length(group)] + 1])) {
        rowRanges(inp)$HL_fragment[group] <-
          paste0(gsub("_O", "", rowRanges(inp)$HL_fragment[group[1] - 1]), "_NA")
        rowRanges(inp)$HL_mean_fragment[group] <-
          rowRanges(inp)$HL_mean_fragment[group[1] - 1]
        group <- row_NA[i + 1]
      } else if (rowRanges(inp)$HL_fragment[group[1] - 1] != rowRanges(inp)$HL_fragment[
        group[length(group)] + 1]) {
        rowRanges(inp)$HL_fragment[group] <-
          paste0("Dc_", as.numeric(
            gsub("Dc_|_O", "", rowRanges(inp)$HL_fragment[group[1] - 1])) + 0.5, "_NA")
        group <- row_NA[i + 1]
      }
    }
  }
  rowRanges(inp)$HL_fragment[is.na(rowRanges(inp)$HL_fragment)] <- paste0("Dc_", as.numeric(
      gsub("Dc_|_O", "", rowRanges(inp)$HL_fragment[!is.na(rowRanges(inp)$HL_fragment)]
           [length(rowRanges(inp)$HL_fragment[!is.na(rowRanges(inp)$HL_fragment)])])) + .5, "_NA")
  inp
}
