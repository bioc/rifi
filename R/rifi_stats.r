# =========================================================================
# rifi_stats        
# -------------------------------------------------------------------------
#' rifi_stats wraps all statistical prediction steps conveniently
#' 
#' rifi_stats wraps the functions: 

#' 1. predict_ps_itss 
#' 2. apply_Ttest_delay
#' 3. apply_ancova
#' 4. apply_event_position
#' 5. apply_t_test
#' 6. fold_change
#' 7. apply_manova
#' 8. apply_t_test_ti 
#' 9. gff3_preprocess
#' 
#' @param inp SummarizedExperiment: the input data frame with correct format.
#' @param dista integer: the maximal distance allowed between two successive
#' fragments. Default is the auto generated value.
#' @param path path: to the directory containing the gff3 file.
#' 
#' @return The SummarizedExperiment object: ID with position, strand, intensity, 
#'  probe_TI, flag, position_segment, delay, half_life, TI_termination_factor, 
#'  delay_fragment, velocity_fragment, intercept, slope, HL_fragment,
#'  HL_mean_fragment, intensity_fragment, intensity_mean_fragment, TU,
#'  TI_termination_fragment, TI_mean_termination_factor, seg_ID, pausing_site,
#'  iTSS_I, ps_ts_fragment, event_ps_itss_p_value_Ttest, p_value_slope,
#'  delay_frg_slope, velocity_ratio, event_duration, event_position, FC_HL,
#'  FC_fragment_HL, p_value_HL, FC_intensity, FC_fragment_intensity,
#'  p_value_intensity, FC_HL_intensity, FC_HL_intensity_fragment, FC_HL_adapted,
#'  synthesis_ratio, synthesis_ratio_event, p_value_Manova, p_value_TI,
#'  TI_fragments_p_value
#' 
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
#' rifi_stats(inp = fragmentation_minimal, dista = 300, 
#' path = gzfile(system.file("extdata", "gff_e_coli.gff3.gz",
#' package = "rifi")))
#' 
#' @export

rifi_stats <- function(inp, dista = 300, path) {
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
