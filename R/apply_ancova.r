#' apply_ancova: is a statistical test to check variances between 2 segments
#' showing pausing site (ps) or internal starting site (ITSS) independently.
#' apply_ancova: is a statistical test to check if fragments showing ps and
#' ITSS events have significant slope using Ancova test.
#' The function uses ancova test. Ancova is applied when the data contains
#' independent variables, dependent variables and covariant variables.
#' In this case, segments are independent variables, position is the dependent
#' variable and the delay is the covariant.
#' The dataframe is prepared as depicted below. The lm fit is applied and
#' p_value is extracted
#'   delay  position  segment
#'                         S1
#'                         S1
#'                         S1
#'                         S1
#'                         S2
#'                         S2
#'                         S2
#'                         S2
#' @param data dataframe: the probe based data frame.
#'
#' @return the probe data frame with the columns regarding statistics:
#' \describe{
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{position:}{The bin/probe specific position}
#'   \item{strand:}{The bin/probe specific strand}
#'   \item{delay:}{The delay value of the bin/probe}
#'   \item{intercept:}{The vintercept of fit through the respective delay
#'   fragment}
#'   \item{slope:}{The slope of the fit through the respective delay
#'   fragment}
#'   \item{pausing_site:}{}
#'   \item{iTSS_I:}{}
#'   \item{ps_ts_fragment:}{}
#'   \item{event_ps_itss_p_value_Ttest:}{}
#'   \item{p_value_slope:}{}
#'   \item{delay_frg_slope:}{}
#'   \item{velocity_ratio:}{}
#' }
#' 
#' @examples
#' data(stats_minimal)
#' apply_ancova(inp = stats_minimal)
#' 
#' @export

apply_ancova <- function(inp) {
    # exclude TUs with 'T', 'NA' or 'O'
    data <- inp[grep("_T|_O|_NA", rowRanges(inp)$TU, invert = TRUE), ]
    data <-
    data[grep("_T|_O|_NA", rowRanges(data)$delay_fragment, invert = TRUE), ]
    # select unique TU
    uniqueTU <- unique(rowRanges(data)$TU)
    rowRanges(inp)$p_value_slope <- NA
    rowRanges(inp)$delay_frg_slope <- NA
    rowRanges(inp)$velocity_ratio <- NA
    
    for (i in seq_along(uniqueTU)) {
      df <- data.frame()
      del_1 <- rowRanges(data)[
        which(rowRanges(data)$TU %in% uniqueTU[i]),"delay_fragment"]
      del_1 <- unique(del_1$delay_fragment)
      if (length(del_1) < 2) {
         next ()
      } else {
        for (j in seq_len(length(del_1) - 1)) {
          seg_1_d <- rowRanges(data)[
            which(rowRanges(data)$delay_fragment %in% del_1[j]), "delay"]
          seg_2_d <- rowRanges(data)[
            which(rowRanges(data)$delay_fragment %in% del_1[j + 1]), "delay"]
          seg_1_p <- rowRanges(data)[
            which(rowRanges(data)$delay_fragment %in% del_1[j]), "position"]
          seg_2_p <- rowRanges(data)[
            which(rowRanges(data)$delay_fragment %in% del_1[j + 1]), "position"]
          
          #if negative strand, data is reversed
          if (unique(strand(data)[which(rowRanges(data)$delay_fragment %in%
                                del_1[1])]) == "-") {
            seg_1_d <- seg_1_d[rev(seq_len(length(seg_1_d)))]
            seg_2_d <- seg_2_d[rev(seq_len(length(seg_2_d)))]
          }
          
          if (length(seg_1_d) == 1 | length(seg_2_d) == 1) {
              next()
          } else {
            df_1 <- cbind.data.frame(seg_1_d, seg_1_p$position)
            df_2 <- cbind.data.frame(seg_2_d, seg_2_p$position)
            colnames(df_1)[7] <- "position"
            colnames(df_2)[7] <- "position"
            # linear model for both segments separately
            model1 <- lm(delay ~ position, data = df_1)
            model2 <- lm(delay ~ position, data = df_2)
            # apply the coefficients of both models to the last point...
            # ...of segment 1 and subtract the distance separating both
            # segments to bring them later... ...to 0
            y1 <-
              model1$coefficients[1] + model1$coefficients[2] *
              last(df_1$position)
            y2 <-
              model2$coefficients[1] + model2$coefficients[2] *
              last(df_1$position)
            dif <- abs(abs(y1) - abs(y2))
            # set the dataframe with delay fragment and delay values
            # subtract the positions from both segments from last
            # position... ...of segment_1
            df_2$position <-
              abs(last(df_1$position) - df_2$position)
            df_1$position <-
              abs(df_1$position[1] - df_1$position)
            # rerun the model with adjusted positions
            model1 <- lm(delay ~ position, data = df_1)
            model2 <- lm(delay ~ position, data = df_2)
            df_1$delay <- df_1$delay - model1$coefficients[1]
            df_2$delay <-
              df_2$delay - (model2$coefficients[1] - dif)
            df <- rbind(df_1, df_2)
            df <-
              cbind(df, c(
                rep("seg.1", times = length(seg_1_d)),
                rep("seg.2", times = length(seg_2_d))
              ))
            colnames(df)[8] <- "seg"
            model1 <-
              lm(delay ~ position + seg + position:seg, data = df)
            tryCatch({
              p_value_slope <- Anova(model1, type = "II")$"Pr(>F)"[3]
              rowRanges(inp)$delay_frg_slope[
                which(rowRanges(inp)$delay_fragment %in% del_1[j])] <-
                paste0(del_1[j], ":", del_1[j + 1])
              rowRanges(inp)$velocity_ratio[
                which(rowRanges(inp)$delay_fragment %in% del_1[j])] <-
              rowRanges(inp)$velocity_fragment[which(rowRanges(inp)$delay_fragment %in% 
                                                       del_1[j + 1])[1]] /
              rowRanges(inp)$velocity_fragment[which(rowRanges(inp)$delay_fragment %in%
                             del_1[j])[1]]
              rowRanges(inp)$p_value_slope[which(rowRanges(inp)$delay_fragment %in%
                           del_1[j])] <-
                p_value_slope
            }, error = function(e) {
            })
          }
        }
      }
    }
    return(inp)
  }
