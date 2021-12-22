#' make_df: converts the input data frame to the template probe based dataframe.
#' 'make_df' creates the probe based data frame with the columns: "ID",
#' "position","strand","intensity","probe_TI" and "flag" from the input
#' dataframe.
#' The replicates are collapsed into their respective means.
#' "ID", "position" and "strand" are taken over as is.
#' "intensity" is the mean intensity from time point 0.
#' "probe_TI" is a value needed for the distribution for the different fitting
#' models.
#' "flag" contains information or the distribution for the different fitting
#' models.
#' Here, probes that don't reach the background level expression are flagged
#' as "_ABG_" ("above background").
#' This is only needed for microarray data and is controlled by the bg
#' parameter.
#' The default for bg = 0, resulting in all probes to be above background
#' (0 is advised for RNAseq data).
#' Probes where all replicates were filtered in the optional filtration step
#' can be fully removed by rm_FLT = TRUE!
#'
#' @param inp data frame: the (checked) input data frame.
#' @param cores integer: the number of assigned cores for the task.
#' @param bg numeric: threshold over which the last timepoint has to be fitted
#' with the above background mode.
#' @param rm_FLT logical: remove IDs where all replicates are marked as
#' filtered. Default is FALSE.
#'
#' @return the probe based data frame:
#'  \describe{
#'     \item{ID:}{The bin/probe specific ID}
#'     \item{position:}{The bin/probe specific position}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{intensity:}{The relative intensity at time point 0}
#'     \item{probe_TI:}{An internal value to determine which fitting model is
#'     applied}
#'     \item{flag:}{Information on which fitting model is applied}
#'       }
#'       
#' @examples
#' data(preprocess_minimal)
#' make_df(inp = preprocess_minimal$input_df, cores = 2, bg = 0,
#' rm_FLT = FALSE)
#' 
#' @export

make_df <- function(inp,
                    cores = 1,
                    bg = 0,
                    rm_FLT = FALSE) {
  num_args <- list(cores, bg)
  names(num_args) <- c("cores", "bg")
  assert(all(unlist(lapply(
    num_args,
    FUN = function(x) {
      (is.numeric(x) &
         length(x) == 1)
    }
  ))),
  paste0("'", names(which(
    unlist(lapply(
      num_args,
      FUN = function(x) {
        (is.numeric(x) &
           length(x) == 1)
      }
    )) == FALSE
  ))[1], "' must be numeric of length one"))
  assert(cores > 0, "'cores' must be a positive integer")
  assert(is.logical(rm_FLT), "'rm_FLT' must be a logical")
  req_cols_inp <- c("0", "ID", "position", "strand", "filtration")
  assert(
    all(req_cols_inp %in% colnames(inp)),
    paste0("'", req_cols_inp[which(!req_cols_inp %in% colnames(inp))],
           "' must be a column in 'inp'!")
  )
  registerDoMC(cores)
  time <-
    as.numeric(names(inp[seq_len(which(colnames(inp) %in% "ID") - 1)]))
  names(time) <- seq_along(time)
  ord <- as.numeric(names(sort(time)))
  ord <- c(ord, seq((max(ord) + 1), ncol(inp), 1))
  inp <- inp[, c(ord)]
  # unique IDs
  unique_ID <- unique(inp$ID)
  # assign a data frame
  # the mean of each replicate is calculated
  tmp_df <-
    foreach(i = seq_along(unique_ID), .combine = rbind) %dopar% {
      tmp <- inp[which(inp$ID %in% unique_ID[i]), ]
      tmp[grep("FLT", tmp$filtration), seq_len(which(colnames(inp) %in% "ID") -
                                                 1)] <- NA
      tmp_mean <- t(colMeans(tmp[seq_len(which(colnames(inp) %in% "ID") - 1)],
                   na.rm = TRUE))
      tmp <- cbind(tmp_mean, tmp[1, (which(colnames(inp) %in% "ID")):
                              (which(colnames(inp) %in% "strand"))])
      tmp
    }
  if (rm_FLT == TRUE) {
    bad_IDs <- tmp_df$ID[which(is.na(tmp_df[1]))]
    tmp_df <- tmp_df[!is.na(tmp_df[1]), ]
    if (length(bad_IDs) > 0) {
      message(paste0(
        "The following IDs have been removed because of rm_FLT == TRUE: ",
        paste0(bad_IDs, collapse = ", ")
      ))
    }
    if (nrow(tmp_df) == 0) {
      stop("All IDs have been removed because of rm_FLT == TRUE")
    }
  }
  val4df <- rep(NA, nrow(tmp_df))
  for (i in seq_len(nrow(tmp_df))) {
    # each row of the input is iterated
    if (is.na(tmp_df[i, 1])) {
      val4df[i] <- NA
    } else {
      val4df[i] <- char_TI(as.numeric(
        tmp_df[i, seq_len(which(colnames(inp) %in% "ID") - 1)]))
    }
  }
  flag <- finding_above_backg(tmp_df, bg_lastTime = bg, time = time)
  out_df <-
    as.data.frame(cbind(
      tmp_df$ID,
      tmp_df$position,
      tmp_df$strand,
      tmp_df[1],
      val4df,
      flag
    ))
  names(out_df) <-
    (c("ID", "position", "strand", "intensity", "probe_TI", "flag"))
  return(out_df)
}
