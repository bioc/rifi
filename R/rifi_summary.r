#' rifi_summary: conveniently wraps and summarize all rifi outputs.
#' Wraps the functions: event_dataframe, dataframe_summary,
#' dataframe_summary_events, dataframe_summary_events_HL_int,
#' dataframe_summary_events_ps_itss, dataframe_summary_events_velocity and
#' dataframe_summary_TI.
#' 
#' @return WIP
#' 
#' @param probe SummarizedExperiment: the input data frame with correct format.
#' @param data_annotation dataframe: the coordinates are extracted from gff3.
#'
#' @seealso `event_dataframe`
#' @seealso `dataframe_summary`
#' @seealso `dataframe_summary_events`
#' @seealso `dataframe_summary_events_HL_int`
#' @seealso `dataframe_summary_events_ps_itss`
#' @seealso `dataframe_summary_events_velocity`
#' @seealso `dataframe_summary_TI`
#' 
#' @examples
#' data(stats_minimal)
#' data(annot_g_minimal)
#' rifi_summary(probe = stats_minimal, data_annotation =
#' annot_g_minimal[[1]])
#' 
#' @export

rifi_summary <- function(probe, data_annotation) {
  res1 <- event_dataframe(data = probe, data_annotation = data_annotation)
  res2 <- dataframe_summary(data = probe, input = res1)
  metadata(probe)$res2 <- res2
  res3 <- dataframe_summary_events(data = probe, data_annotation =
                                     data_annotation)
  metadata(probe)$res3 <- res3
  res4 <- dataframe_summary_events_HL_int(data = probe, data_annotation =
                                            data_annotation)
  metadata(probe)$res4 <- res4
  res5 <- dataframe_summary_events_ps_itss(data = probe, data_annotation =
                                             data_annotation)
  metadata(probe)$res5 <- res5
  res6 <- dataframe_summary_events_velocity(data = probe, data_annotation =
                                              data_annotation)
  metadata(probe)$res6 <- res6
  res7 <- dataframe_summary_TI(data = probe, input = res1)
  metadata(probe)$res7 <- res7
  probe
}
