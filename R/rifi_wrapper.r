#' rifi_wrapper: conveniently wraps all functions included on rifi workflow.
#' Wraps the functions: check_input, rifi_preprocess, rifi_fit, rifi_penalties,
#' rifi_fragmentation, rifi_stats, rifi_summary and rifi_visualization.
#'
#' @param inp data frame: the input data frame with correct format.
#' @param cores integer: the number of assigned cores for the task.
#' @param gff path: path to an annotation file in gff format
#' @param bg numeric: threshold over which the last time point has to be to be
#' fitted with the above background mode.
#' @param restr numeric: a parameter that restricts the freedom of the fit to
#' avoid wrong TI-term_factors, ranges from 0 to 0.2
#'
#' @return All intermediate objects
#'
#' @seealso `check_input`
#' @seealso `rifi_preprocess`
#'
#' @export


rifi_wrapper <- function(inp, cores, gff, bg, restr) {
  prepro <- rifi_preprocess(inp = inp, cores = cores, bg = bg)
  probe <-
    rifi_fit(prepro[[2]],
             prepro[[1]],
             cores = cores,
             viz = FALSE,
             restr = restr)
  pen <- rifi_penalties(
    probe,
    details = FALSE,
    viz = FALSE,
    top_i = 0,
    cores = cores
  )
  probe_fra <-
    rifi_fragmentation(probe, cores = cores, logbook = pen)
  annot <- gff3_preprocess(gff)
  probe_sta <- rifi_stats(probe_fra)
  probe_summary <-
    rifi_summary(probe_sta, data_annotation = annot[[1]])
  rifi_visualization(data = probe_sta,
                     genomeLength = annot[[2]],
                     annot = annot[[1]])
  res <- list(prepro, probe, pen, probe_fra, probe_sta, probe_summary)
  return(res)
}
