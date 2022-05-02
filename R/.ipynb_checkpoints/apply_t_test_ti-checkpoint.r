# =========================================================================
# apply_t_test_ti   Compares the mean of two neighboring TI fragments
#'                  within the same TU.
# -------------------------------------------------------------------------
#'
#'
#' apply_t_test_ti uses the statistical t_test to check if two neighboring 
#' TI fragments are significant.
#'
#' @param inp SummarizedExperiment: the input data frame with correct format.
#'
#' @return the SummarizedExperiment with the columns regarding statistics:
#' \describe{
#'   \item{p_value_TI:}{Integer, the p_value added to the input}
#'   \item{TI_fragments_p_value:}{String, the fragments subjected to statistical test}
#' }
#'
#' @examples
#' data(stats_minimal)
#' apply_t_test_ti(inp = stats_minimal)
#'
#' @export
#'
apply_t_test_ti <- function(inp) {
  #new columns are added
  rowRanges(inp)$p_value_TI <- NA
  rowRanges(inp)$TI_fragments_p_value <- NA
  #grep TI fragments excluding outliers
  data_1 <-
    rowRanges(inp)[!grepl("_T|_O|_NA", rowRanges(inp)$TI_termination_fragment),]
  #select unique TUs
  unique_TU <- unique(rowRanges(inp)$TU)
  #exclude TU terminales, outliers and NAs
  unique_TU <- na.omit(unique_TU[!grepl("_T|_O|_NA", unique_TU)])
  for (i in seq_along(unique_TU)) {
    #grep only TI fragments
    tu <- data_1[which(unique_TU[i] == data_1$TU),
                 c(
                   "ID",
                   "position",
                   "flag",
                   "TI_termination_fragment",
                   "TI_termination_factor",
                   "intensity"
                 )]
    tu <- tu[!is.na(tu$TI_termination_fragment),]
    ti <- tu[grep("_TI_", tu$flag),]
    #adjust the fragments for t-test
    if (length(na.omit(ti$TI_termination_fragment)) == 0) {
      next ()
    } else {
      ti_frag <- unique(tu$TI_termination_fragment)
      if (length(ti_frag) > 1) {
        for (k in seq_len(length(ti_frag) - 1)) {
          seg1 <- tu[which(ti_frag[k] == tu$TI_termination_fragment),
                     "TI_termination_factor"]
          seg2 <-
            tu[which(ti_frag[k + 1] == tu$TI_termination_fragment),
               "TI_termination_factor"]
          if (length(seg1) < 2 | length(seg2) < 2) {
            next ()
          }
          tryCatch({
            #t-test
            ti_test <- t.test(seg1$TI_termination_factor,
                              seg2$TI_termination_factor,
                              alternative = "two.sided",
                              var.equal = FALSE)
            #add 2 columns, fragments column and p_value from t-test
            p_value_tiTest <- ti_test[[3]]
            rowRanges(inp)$TI_fragments_p_value[
              which(ti_frag[k] == rowRanges(inp)$TI_termination_fragment)] <-
              paste0(ti_frag[k], ":", ti_frag[k + 1])
            rowRanges(inp)$p_value_TI[
              which(ti_frag[k] == rowRanges(inp)$TI_termination_fragment)] <-
              p_value_tiTest
          }, error = function(e) {
          })
        }
      }
    }
  }
  return(inp)
}
