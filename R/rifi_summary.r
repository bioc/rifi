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
#' rifi_summary(inp = stats_minimal, data_annotation = 
#' metadata(stats_minimal)$annot[[1]])
#' 
#' @export

rifi_summary <- function(inp, data_annotation) {
  res1 <- event_dataframe(data = inp, data_annotation = data_annotation)
  res2 <- dataframe_summary(data = inp, input = res1)
  metadata(inp)$res2 <- res2
  res3 <- dataframe_summary_events(data = inp, data_annotation =
                                     data_annotation)
  metadata(inp)$res3 <- res3
  res4 <- dataframe_summary_events_HL_int(data = inp, data_annotation =
                                            data_annotation)
  metadata(inp)$res4 <- res4
  res5 <- dataframe_summary_events_ps_itss(data = inp, data_annotation =
                                             data_annotation)
  metadata(inp)$res5 <- res5
  res6 <- dataframe_summary_events_velocity(data = inp, data_annotation =
                                              data_annotation)
  metadata(inp)$res6 <- res6
  res7 <- dataframe_summary_TI(data = inp, input = res1)
  metadata(inp)$res7 <- res7
  inp
}
