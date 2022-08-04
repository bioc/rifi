#' =========================================================================
#' apply_t_test_ti   
#' -------------------------------------------------------------------------
#' apply_t_test_ti compares the mean of two neighboring TI fragments within the
#' same TU.
#'
#' apply_t_test_ti uses the statistical t_test to check if two neighboring 
#' TI fragments are significant.
#'
#' @param inp SummarizedExperiment: the input data frame with correct format.
#'
#' @return the SummarizedExperiment with the columns regarding statistics:
#' \describe{
#'   \item{ID:}{The bin/probe specific ID.}
#'   \item{position:}{The bin/probe specific position.}
#'   \item{strand:}{The bin/probe specific strand.}
#'   \item{intensity:}{The relative intensity at time point 0.}
#'   \item{probe_TI:}{An internal value to determine which fitting model
#'   is applied.}
#'   \item{flag:}{Information on which fitting model is applied.}
#'   \item{position_segment:}{The position based segment.}
#'   \item{delay:}{The delay value of the bin/probe.}
#'   \item{half_life:}{The half-life of the bin/probe.}
#'   \item{TI_termination_factor:}{String, the factor of TI fragment.}
#'   \item{delay_fragment:}{The delay fragment the bin belongs to.}
#'   \item{velocity_fragment:}{The velocity value of the respective delay
#'   fragment.}
#'   \item{intercept:}{The vintercept of fit through the respective delay
#'   fragment.}
#'   \item{slope:}{The slope of the fit through the respective delay fragment.}
#'   \item{HL_fragment:}{The half-life fragment the bin belongs to.}
#'   \item{HL_mean_fragment:}{The mean half-life value of the respective
#'   half-life fragment.}
#'   \item{intensity_fragment:}{The intensity fragment the bin belongs to.}
#'   \item{intensity_mean_fragment:}{The mean intensity value of the respective
#'   intensity fragment.}
#'   \item{TU:}{The overarching transcription unit.}
#'   \item{TI_termination_fragment:}{The TI fragment the bin belongs to.}
#'   \item{TI_mean_termination_factor:}{The mean termination factor of the
#'   respective TI fragment.}
#'   \item{seg_ID:}{The combined ID of the fragment.}
#'   \item{pausing_site:}{presence of pausing site indicated by +/-.}
#'   \item{iTSS_I:}{presence of iTSS_I indicated by +/-.}
#'   \item{ps_ts_fragment:}{The fragments involved in pausing site or iTSS_I.}
#'   \item{event_duration:}{Integer, the duration between two delay fragments.}
#'   \item{event_ps_itss_p_value_Ttest:}{p_value of pausing site or iTSS_I.}
#'   \item{p_value_slope:}{Integer, the p_value added to the inp.}
#'   \item{delay_frg_slope:}{Integer, the slope value of the fit through the 
#'   respective delay fragment.}
#'   \item{velocity_ratio:}{Integer, the ratio value of velocity from 2 delay 
#'   fragments.}
#'   \item{event_position:}{Integer, position of the event added to the input.}
#'   \item{FC_HL:}{Integer, the fold change value of 2 HL fragments.}
#'   \item{FC_fragment_HL:}{String, the fragments corresponding to HL fold 
#'   change.}
#'   \item{p_value_HL:}{Integer, the p_value added to the input of 2 HL 
#'   fragments.}
#'   \item{FC_intensity:}{Integer, the fold change value of 2 intensity 
#'   fragments.}
#'   \item{FC_fragment_intensity:}{String, the fragments corresponding to 
#'   intensity fold change.}
#'   \item{p_value_intensity:}{Integer, the p_value added to the input of 2 
#'   intensity fragments.}
#'   \item{synthesis_ratio:}{Integer, the value correspomding to synthesis rate.}
#'   \item{synthesis_ratio_event:}{String, the event assigned by synthesis 
#'   rate either Termination or iTSS.}
#'   \item{FC_HL_intensity:}{Integer, the value corresponding to HL and
#'   intensity fold change.}
#'   \item{FC_HL_intensity_fragment:}{String, the fragments corresponding to 
#'   intensity and HL fold change.}
#'   \item{FC_HL_adapted:}{Integer, the fold change of half-life/ fold change 
#'   of intensity,position of the half-life fragment is adapted to intensity 
#'   fragment.}
#'   \item{p_value_Manova:}{Integer, the p_value added to the input.}
#'   \item{p_value_TI:}{Integer, the p_value added to the input.}
#'   \item{TI_fragments_p_value:}{String, the fragments subjected to statistical
#'   test.}
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
