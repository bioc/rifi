#' apply_Ttest_delay: is a statistical test to check the significance
#' of the point between 2 segments showing pausing site (ps) and
#' internal starting site (ITSS) independently.
#' apply_Ttest_delay uses t-test. The last point from the first segment
#' and the first point from the second segment are selected and added
#' to the residuals of each model. The sum is subjected to t-test.
#' @param data dataframe: the probe based dataframe.
#' 
#' @return the probe data frame with the columns regarding statistics:
#' \describe{
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{position:}{The bin/probe specific position}
#'   \item{strand:}{The bin/probe specific strand}
#'   \item{delay:}{The delay value of the bin/probe}
#'   \item{delay_fragment:}{The delay fragment the bin belongs to}
#'   \item{pausing_site:}{}
#'   \item{iTSS_I:}{}
#'   \item{ps_ts_fragment:}{}
#'   \item{event_ps_itss_p_value_Ttest:}{}
#' }
#'
#' @examples
#' data(stats_minimal)
#' apply_Ttest_delay(inp = stats_minimal)
#' 
#' @export
#' 
apply_Ttest_delay <- function(inp) {
  event_1 <- which(rowRanges(inp)$pausing_site == "+")
  event_2 <- which(rowRanges(inp)$iTSS_I == "+")
  event <- c(event_1, event_2)
  rowRanges(inp)$event_ps_itss_p_value_Ttest <- NA
  for (i in seq_len(length(event) - 1)) {
    ps <- unlist(str_split(rowRanges(inp)$ps_ts_fragment[event[i]], ":"))
    seg_1_d <-
      rowRanges(inp)$delay[which(rowRanges(inp)$delay_fragment %in% ps[1])]
    seg_2_d <-
      rowRanges(inp)$delay[which(rowRanges(inp)$delay_fragment %in% ps[2])]
    seg_1_p <-
      rowRanges(inp)$position[which(rowRanges(inp)$delay_fragment %in% ps[1])]
    seg_2_p <-
      rowRanges(inp)$position[which(rowRanges(inp)$delay_fragment %in% ps[2])]
    # in case of negative strand, positions are shifted
    if (unique(strand(inp)[which(rowRanges(inp)$delay_fragment %in% ps[1])]) == "-") {
      seg_1_d <- seg_1_d[rev(seq_len(length(seg_1_d)))]
      seg_2_d <- seg_2_d[rev(seq_len(length(seg_2_d)))]
    }
    # in case of segment has a length of one, its not considered for...
    # ...further analysis
    if (length(seg_1_d) == 1 | length(seg_2_d) == 1) {
      (next)()
    } else {
      df_1 <- cbind.data.frame(seg_1_d, seg_1_p)
      df_2 <- cbind.data.frame(seg_2_d, seg_2_p)
      colnames(df_1) <- c("delay", "position")
      colnames(df_2) <- c("delay", "position")
      # linear model for both segments separately
      model1 <- lm(delay ~ position, data = df_1)
      model2 <- lm(delay ~ position, data = df_2)
      # select the last point from the first fragment and the first...
      # ...point from the second fragment
      del_p1 <- last(df_1$delay)
      del_p2 <- df_2$delay[1]
      # sum up absolute residuals to the point indicated above
      res_model1 <- abs(residuals(model1)) + del_p1
      res_model2 <- abs(residuals(model2)) + del_p2
      tryCatch({
        # run t-test between both sum of residuals and the point...
        # ...selected for statistics
        t_h <-
          t.test(res_model1,
                 res_model2,
                 alternative = "two.sided",
                 var.equal = FALSE)
        # extract the p_value from t-test
        p_value_Ttest <- t_h[[3]]
        rowRanges(inp)$event_ps_itss_p_value_Ttest[event[i]] <-
          p_value_Ttest
      }, error = function(e) {
      })
    }
  }
  return(inp)
}
