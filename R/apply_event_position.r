#' apply_event_position: is a short version of apply_Ttest_delay function
#' to extract event time duration as pausing site or iTSS happens.
#'
#' apply_event_position adds a new column with the duration.
#'
#' @param data dataframe: the probe based dataframe.
#' 
#' @return the probe data frame with the columns regarding statistics:
#' \describe{
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{position:}{The bin/probe specific position}
#'   \item{strand:}{The bin/probe specific strand}
#'   \item{delay:}{The delay value of the bin/probe}
#'   \item{pausing_site:}{}
#'   \item{iTSS_I:}{}
#'   \item{ps_ts_fragment:}{}
#'   \item{event_ps_itss_p_value_Ttest:}{}
#'   \item{p_value_slope:}{}
#'   \item{delay_frg_slope:}{}
#'   \item{velocity_ratio:}{}
#'   \item{event_position:}{}
#' }
#'
#' @examples
#' data(stats_minimal)
#' apply_event_position(data = stats_minimal)
#' 
#' @export
#'

apply_event_position <- function(data) {
  event_1 <- which(data[, "pausing_site"] == "+")
  event_2 <- which(data[, "iTSS_I"] == "+")
  event <- c(event_1, event_2)

  data[, "event_position"] <- NA
  for (i in seq_along(event)) {
    ps <- data[event[i], "ps_ts_fragment"]
    ps_1 <- unlist(str_split(ps, ":"))
    ps_2 <- ps_1[2]
    ps_1 <- ps_1[1]
    seg_1_d <-
      data[which(data$delay_fragment %in% ps_1), "delay"]
    seg_2_d <-
      data[which(data$delay_fragment %in% ps_2), "delay"]
    seg_1_p <-
      data[which(data$delay_fragment %in% ps_1), "position"]
    seg_2_p <-
      data[which(data$delay_fragment %in% ps_2), "position"]
    # in case of negative strand, the positions are shifted
    if (unique(data[which(data$delay_fragment %in%
                          ps_1), "strand"]) == "-") {
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
      data[which(data$ps_ts_fragment %in% ps),
           "event_position"] <-
        (data[last(which(data$delay_fragment ==
                           ps_1)), "position"] +
           data[which(data$delay_fragment == ps_2), "position"][1]) / 2
    }
  }
  return(data)
}
