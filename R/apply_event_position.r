#' =========================================================================
#' apply_event_position 
#' -------------------------------------------------------------------------
#' apply_event_position extracts event time duration for pausing site or
#' iTSS
#'
#' apply_event_position is a short version of apply_Ttest_delay function to
#' extract event time duration for pausing site or iTSS. Its adds a new column 
#' with the duration.
#'
#' @param inp SummarizedExperiment: the input data frame with correct 
#' format.
#' 
#' @return The SummarizedExperiment with the columns regarding statistics:
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
#' }
#'
#' @examples
#' data(stats_minimal)
#' apply_event_position(inp = stats_minimal)
#' 
#' @export
#'

apply_event_position <- function(inp) {
  event_1 <- which(rowRanges(inp)$pausing_site == "+")
  event_2 <- which(rowRanges(inp)$iTSS_I == "+") 
  event <- c(event_1, event_2)
  rowRanges(inp)$event_position <- NA
  
  for (i in seq_along(event)) {
    ps <- rowRanges(inp)$ps_ts_fragment[event[i]]
    ps_1 <- unlist(str_split(ps, ":"))
    ps_2 <- ps_1[2]
    ps_1 <- ps_1[1]
    seg_1_d <-
      rowRanges(inp)$delay[which(rowRanges(inp)$delay_fragment %in% ps_1)]
    seg_2_d <-
      rowRanges(inp)$delay[which(rowRanges(inp)$delay_fragment %in% ps_2)]
    seg_1_p <-
      rowRanges(inp)$position[which(rowRanges(inp)$delay_fragment %in% ps_1)]
    seg_2_p <-
      rowRanges(inp)$position[which(rowRanges(inp)$delay_fragment %in% ps_2)]
    # in case of negative strand, the positions are shifted
    if (unique(strand(inp)[which(rowRanges(inp)$delay_fragment %in%
                          ps_1)]) == "-") {
      seg_1_d <- seg_1_d[rev(seq_len(length(seg_1_d)))]
      seg_2_d <- seg_2_d[rev(seq_len(length(seg_2_d)))]
    }
    # if segment had a length of one, its not considered for further
    # analysis
    if (length(seg_1_d) == 1 | length(seg_2_d) == 1) {
      (next)()
    } else {
      df_1 <- cbind.data.frame(seg_1_d, seg_1_p)
      df_2 <- cbind.data.frame(seg_2_d, seg_2_p)
      colnames(df_1) <- c("delay", "position")
      colnames(df_2) <- c("delay", "position")
      # select the last point from the first fragment and the first...
      # ...point from the second fragment
      del_p1 <- last(df_1$delay)
      del_p2 <- df_2$delay[1]
      rowRanges(inp)$event_position[
        which(rowRanges(inp)$ps_ts_fragment %in% ps)] <-
        (rowRanges(inp)$position[last(which(rowRanges(inp)$delay_fragment ==
                           ps_1))] +
           rowRanges(inp)$position[
             which(rowRanges(inp)$delay_fragment == ps_2)][1]) / 2
    }
  }
  return(inp)
}
