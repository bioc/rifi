# =========================================================================
# make_df      Adds important columns to the SummarizedExperiment object
# -------------------------------------------------------------------------
#'
#' 
#' 'make_df' adds to the SummarizedExperiment object with the columns:
#' "intensity", "probe_TI" and "flag".

#' The replicates are collapsed into their respective means.

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
#' can be fully removed by rm_FLT = TRUE! If you wish to keep all information
#' in the assay set to FALSE!
#'
#' @param inp SummarizedExperiment: the (checked) input.
#' @param cores integer: the number of assigned cores for the task.
#' @param bg numeric: threshold over which the last timepoint has to be fitted
#' with the above background mode.
#' @param rm_FLT logical: remove IDs where all replicates are marked as
#' filtered. Default is FALSE.
#'
#' @return the SummarizedExperiment object: with intensity, probe_TI and
#' flag added to the rowRanges.
#'       
#' @examples
#' data(preprocess_minimal)
#' make_df(inp = preprocess_minimal, cores = 2, bg = 0, rm_FLT = TRUE)
#' 
#' @export

make_df <- function(inp,
                    cores = 1,
                    bg = 0,
                    rm_FLT = TRUE) {
  assert(is.numeric(bg) & length(bg) == 1,
         "bg must be numeric of length one")
  assert(is.logical(rm_FLT),
         "'rm_FLT' must be a logical")
  registerDoMC(cores)
  #order the input
  inp <- inp_order(inp)
  # unique time points
  time <- unique(metadata(inp)$timepoints)
  # the mean of each replicate is calculated
  tmp_df <-
    foreach(i = seq_along(time), .combine = cbind) %dopar% {
      tmp_inp <-inp[,colData(inp)$timepoint == time[i]]
      assay(tmp_inp)[decode_FLT(tmp_inp)] <- NA
      tmp <- rowMeans(assay(tmp_inp), na.rm = TRUE)
      tmp
    }
  if (rm_FLT == TRUE) {
    remv <- is.na(tmp_df)[,1]
    tmp_df <- tmp_df[!remv, ]
    inp <- inp[!remv, ]
    if (nrow(tmp_df) == 0) {
      stop("All probes have been removed because of rm_FLT == TRUE")
    }
  }
  val4df <- rep(NA, nrow(tmp_df))
  for (i in seq_len(nrow(tmp_df))) {
    # each row of the input is iterated
    if (is.na(tmp_df[i, 1])) {
      val4df[i] <- NA
    } else {
      val4df[i] <- char_TI(as.numeric(
        tmp_df[i, ]))
    }
  }
  flag <- finding_above_backg(df = tmp_df, bg = bg)
  rowRanges(inp)$intensity <- tmp_df[,1]
  rowRanges(inp)$probe_TI <- val4df
  rowRanges(inp)$flag <- flag
  inp
}
