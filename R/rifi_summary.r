#' rifi_summary: conveniently wraps and summarize all rifi outputs.
#' Wraps the functions: event_dataframe, dataframe_summary,
#' dataframe_summary_events, dataframe_summary_events_HL_int,
#' dataframe_summary_events_ps_itss, dataframe_summary_events_velocity and
#' dataframe_summary_TI.
#' 
#' @return WIP
#' 
#' @param probe data frame: the probe based data frame.
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
  res3 <- dataframe_summary_events(data = probe, data_annotation =
                                     data_annotation)
  res4 <- dataframe_summary_events_HL_int(data = probe, data_annotation =
                                            data_annotation)
  res5 <- dataframe_summary_events_ps_itss(data = probe, data_annotation =
                                             data_annotation)
  res6 <- dataframe_summary_events_velocity(data = probe, data_annotation =
                                              data_annotation)
  res7 <- dataframe_summary_TI(data = probe, input = res1)
  res <- list(res2[[1]], res2[[2]], res3, res4, res5, res6, res7)
  res
}
