# =========================================================================
# apply_event_position Extract event time duration for pausing site or iTSS
# -------------------------------------------------------------------------
#'
#'
#' apply_event_position is a short version of apply_Ttest_delay function
#' to extract event time duration for pausing site or iTSS.
#'
#' apply_event_position adds a new column with the duration.
#'
#' @param inp SummarizedExperiment: the input data frame with correct format.
#' 
#' @return the SummarizedExperiment with the columns regarding statistics:
#' \describe{
#'   \item{position:}{Integer, the bin/probe specific position}
#'   \item{delay:}{Integer, the delay value of the bin/probe}
#'   \item{pausing_site:}{Boolean, pausing site event if happend, either + or -}
#'   \item{iTSS_I:}{Boolean, iTSS_I event if happend, either + or -}
#'   \item{ps_ts_fragment:}{String, the fragment where the event happend}
#'   \item{event_position:}{Integer, position of the event added to the input}
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
