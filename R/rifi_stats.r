#' rifi_stats: conveniently wraps all statistical prediction steps.
#' Wraps the functions: predict_ps_itss, apply_Ttest_delay, apply_ancova,
#' apply_event_position, apply_t_test, fold_change, apply_manova, 
#' apply_t_test_ti and gff3_preprocess.
#' 
#' @param inp SummarizedExperiment: the input data frame with correct format.
#' @param dista integer: the maximal distance allowed between two successive
#' fragments. Default is the auto generated value.
#' @param path path: to the directory containing the gff3 file.
#' 
#' @return the probe data frame with the columns regarding statistics:
#' \describe{
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{position:}{The bin/probe specific position}
#'   \item{strand:}{The bin/probe specific strand}
#'   \item{intensity:}{The relative intensity at time point 0}
#'   \item{probe_TI:}{An internal value to determine which fitting model is
#'   applied}
#'   \item{flag:}{Information on which fitting model is applied}
#'   \item{position_segment:}{The position based segment}
#'   \item{delay:}{The delay value of the bin/probe}
#'   \item{half_life:}{The half-life of the bin/probe}
#'   \item{TI_termination_factor:}{The termination factor of the bin/probe}
#'   \item{delay_fragment:}{The delay fragment the bin belongs to}
#'   \item{velocity_fragment:}{The velocity value of the respective delay
#'   fragment}
#'   \item{intercept:}{The vintercept of fit through the respective delay
#'   fragment}
#'   \item{slope:}{The slope of the fit through the respective delay fragment}
#'   \item{HL_fragment:}{The half-life fragment the bin belongs to}
#'   \item{HL_mean_fragment:}{The mean half-life value of the respective
#'   half-life fragment}
#'   \item{intensity_fragment:}{The intensity fragment the bin belongs to}
#'   \item{intensity_mean_fragment:}{The mean intensity value of the respective
#'   intensity fragment}
#'   \item{TU:}{The overarching transcription unit}
#'   \item{TI_termination_fragment:}{The TI fragment the bin belongs to}
#'   \item{TI_mean_termination_factor:}{The mean termination factor of the
#'   respective TI fragment}
#'   \item{seg_ID:}{The combined ID of the fragment}
#'   \item{pausing_site:}{}
#'   \item{iTSS_I:}{}
#'   \item{ps_ts_fragment:}{}
#'   \item{event_ps_itss_p_value_Ttest:}{}
#'   \item{p_value_slope:}{}
#'   \item{delay_frg_slope:}{}
#'   \item{velocity_ratio:}{}
#'   \item{event_duration:}{}
#'   \item{event_position:}{}
#'   \item{FC_HL:}{}
#'   \item{FC_fragment_HL:}{}
#'   \item{p_value_HL:}{}
#'   \item{FC_intensity:}{}
#'   \item{FC_fragment_intensity:}{}
#'   \item{p_value_intensity:}{}
#'   \item{FC_HL_intensity:}{}
#'   \item{FC_HL_intensity_fragment:}{}
#'   \item{FC_HL_adapted:}{}
#'   \item{synthesis_ratio:}{}
#'   \item{synthesis_ratio_event:}{}
#'   \item{p_value_Manova:}{}
#'   \item{p_value_TI:}{}
#'   \item{TI_fragments_p_value:}{}
#' }
#' 
#' @seealso `predict_ps_itss`
#' @seealso `apply_Ttest_delay`
#' @seealso `apply_ancova`
#' @seealso `apply_event_position`
#' @seealso `apply_t_test`
#' @seealso `fold_change`
#' @seealso `apply_manova`
#' @seealso `apply_t_test_ti`
#' @seealso `gff3_preprocess`
#' 
#' @examples
#' data(fragmentation_minimal)
#' rifi_stats(inp = fragmentation_minimal, dista = 300, path)
#' 
#' @export

rifi_stats <- function(inp, dista = 300, gff) {
  message("running predict_ps_itss...")
  probe <- predict_ps_itss(inp = inp, maxDis = dista)
  message("running apply_Ttest_delay...")
  probe <- apply_Ttest_delay(inp = probe)
  message("running apply_ancova...")
  probe <- apply_ancova(inp = probe)
  message("running apply_event_position...")
  probe <- apply_event_position(inp = probe)
  message("running apply_t_test...")
  probe <- apply_t_test(inp = probe)
  message("running fold_change...")
  probe <- fold_change(inp = probe)
  message("running apply_manova...")
  probe <- apply_manova(inp = probe)
  message("running apply_t_test_ti...")
  probe <- apply_t_test_ti(inp = probe)
  metadata(probe)$annot <- gff3_preprocess(path)
  probe
}
