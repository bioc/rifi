# =========================================================================
# apply_t_test   Uses the statistical t_test to check if the fold-change 
#'                  of half-life (HL) fragments and the fold-change
#'                  intensity fragments respectively are significant.
# -------------------------------------------------------------------------
#'
#' apply_t_test compares the mean of two neighboring fragments within the
#' same TU to check if the fold-change is significant.

#' Fragments with distance above threshold are not subjected to t-test.
#' Dataframes with less than 3 rows are excluded.
#'
#' The functions used are:
#'
#' 1. fragment_function: checks number of fragments inside TU, less
#' than 2 are excluded otherwise they are gathered for analysis.

#' 2. t_test_function: excludes dataframes with less than 3 rows,
#' makes fold-change and apply t-test, assign fragments names
#' and ratio, add columns with the corresponding p_values.
#'
#' @param inp SummarizedExperiment: the input data frame with correct format.
#' @param threshold integer: threshold.
#' 
#' @return the SummarizedExperiment with the columns regarding statistics:
#' \describe{
#'   \item{FC_HL:}{Integer, the fold change value of 2 HL fragments}
#'   \item{FC_fragment_HL:}{String, the fragments corresponding to HL fold change}
#'   \item{p_value_HL:}{Integer, the p_value added to the input of 2 HL fragments}
#'   \item{FC_intensity:}{Integer, the fold change value of 2 intensity fragments}
#'   \item{FC_fragment_intensity:}{String, the fragments corresponding to intensity fold change}
#'   \item{p_value_intensity:}{Integer, the p_value added to the input of 2 intensity fragments}
#' }
#'
#' @examples
#' data(stats_minimal)
#' apply_t_test(inp = stats_minimal, threshold = 300)
#' 
#' @export

apply_t_test <- function(inp, threshold = 300) {
  rowRanges(inp)$FC_fragment_HL <- NA
  rowRanges(inp)$FC_HL <- NA 
  rowRanges(inp)$p_value_HL <- NA
  rowRanges(inp)$FC_fragment_intensity <- NA
  rowRanges(inp)$FC_intensity <- NA 
  rowRanges(inp)$p_value_intensity <- NA
  uniqueTU <- unique(rowRanges(inp)$TU)
  uniqueTU <- uniqueTU[grep("_NA|_T", uniqueTU, invert = TRUE)]
  for (i in seq_along(uniqueTU)) {
    # select ID, position, HL, HL fragments, intensity and intensity
    # fragments for the corresponding TU
    tu <-
      as.data.frame(rowRanges(inp)[which(rowRanges(inp)$TU %in% uniqueTU[i]), c(
        "ID",
        "position",
        "half_life",
        "TU",
        "HL_fragment",
        "intensity",
        "intensity_fragment",
        "HL_mean_fragment",
        "intensity_mean_fragment"
      )])
    # HL and intensity segments in the TU
    hl_segs <-
      tu[grep(paste0("\\Dc_\\d+", "$"), 
                                     tu$HL_fragment), "HL_fragment"]
    int_segs <- tu[grep(paste0("\\I_\\d+", "$"),
                        tu$intensity_fragment), "intensity_fragment"]
    hl_segs <- fragment_function(hl_segs)
    int_segs <- fragment_function(int_segs)
    # loop into all HL segments and apply t_test between consecutive segments
   tryCatch({
      inp <-
        t_test_function(
        data = inp,
        seg = hl_segs,
        param = "half_life",
        o = "HL",
        tu = tu,
        threshold = threshold
      )
        }, error = function(e) {
      })
    if(length(int_segs) > 1){
    inp <-
      t_test_function(
        data = inp,
        seg = int_segs,
        param = "intensity",
        o = "intensity",
        tu = tu,
        threshold = threshold
      )
    }
  }
  return(inp)
}
