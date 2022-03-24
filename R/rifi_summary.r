#' rifi_summary: conveniently wraps and summarize all rifi outputs.
#' Wraps the functions: event_dataframe, dataframe_summary,
#' dataframe_summary_events, dataframe_summary_events_HL_int,
#' dataframe_summary_events_ps_itss, dataframe_summary_events_velocity and
#' dataframe_summary_TI.
#' 
#' @param inp SummarizedExperiment: the input data frame with correct format.
#' @param data_annotation dataframe: gff3 dataframe after processing.
#'
#' @return WIP
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
#' if(!require(SummarizedExperiment)){
#' suppressPackageStartupMessages(library(SummarizedExperiment))
#' }
#' rifi_summary(inp = stats_minimal, data_annotation = 
#' metadata(stats_minimal)$annot[[1]])
#' 
#' @export

rifi_summary <- function(inp, data_annotation = metadata(inp)$annot[[1]]) {
  res <- event_dataframe(data = inp, data_annotation = data_annotation)
  res1 <- dataframe_summary(data = inp, input = res)
  metadata(inp)$dataframe_summary_1 <- res1[[1]]
  metadata(inp)$dataframe_summary_2 <- res1[[2]]
  res3 <- dataframe_summary_events(data = inp, data_annotation =
                                     data_annotation)
  metadata(inp)$dataframe_summary_events <- res3
  res4 <- dataframe_summary_events_HL_int(data = inp, data_annotation =
                                            data_annotation)
  metadata(inp)$dataframe_summary_events_HL_int <- res4
  res5 <- dataframe_summary_events_ps_itss(data = inp, data_annotation =
                                             data_annotation)
  metadata(inp)$dataframe_summary_events_ps_itss <- res5
  res6 <- dataframe_summary_events_velocity(data = inp, data_annotation =
                                              data_annotation)
  metadata(inp)$dataframe_summary_events_velocity <- res6
  res7 <- dataframe_summary_TI(data = inp, input = res)
  metadata(inp)$dataframe_summary_TI <- res7
  inp
}
