#' =========================================================================
#' predict_ps_itss 
#' -------------------------------------------------------------------------
#'predict_ps_itss predicts pausing sites (ps) and internal starting sites
#' (ITSS) between delay fragments.
#' 
#' predict_ps_itss predicts ps and ITSS within the same TU. Neighboring delay
#' segments are compared to each other by positioning the intercept of the 
#' second segment into the first segment using slope and intercept coefficients.
#' 
#' predict_ps_itss uses 3 steps to identify ps and ITSS:
#' 
#' 1. select unique TU.
#' 
#' 2. select from the input dataframe the columns: ID, position, strand, delay.
#' delay fragment, TU and slope coordinates, velocity_fragment and intercept.
#' 
#' 3. select delay segments in the TU.
#' 
#' 4. loop into all delay segments and estimate the coordinates of the last
#' point of the first segment using the coefficients of the second segment 
#' and vice versa. We get two predicted positions, the difference between
#' them is compared to the threshold.
#' 
#' In case the strand is "-", additional steps are added:
#' 
#' The positions of both segments are ordered from the last position to the
#' first one.
#' 
#' All positions are merged in one column and subtracted from the maximum
#' position. the column is split in 2. The first and second correspond to
#' the positions of the first and second segments respectively.
#' 
#' Both segments are subjected to lm fit and the positions predicted are used
#' on the same way as the opposite strand.
#' 
#' If the difference between the positions predicted is lower than negative
#' threshold, ps is assigned otherwise, and if the difference is higher than
#' the positive threshold, ITSS is assigned.
#'
#' @param inp SummarizedExperiment: the input data frame with correct format.
#' @param maxDis integer: the maximal distance allowed between two successive
#' fragments.
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
#'   }
#' }
#' 
#' @examples
#' data(fragmentation_minimal)
#' predict_ps_itss(inp = fragmentation_minimal, maxDis = 300) 
#' @export

predict_ps_itss <- function(inp, maxDis = 300) {
  rowRanges(inp)$pausing_site <- "-"
  rowRanges(inp)$iTSS_I <- "-"
  rowRanges(inp)$ps_ts_fragment <- NA
  rowRanges(inp)$event_duration <- NA

  # select unique TUs
  uniqueTU <- unique(rowRanges(inp)$TU)
  uniqueTU <- uniqueTU[grep("_NA|_T", uniqueTU, invert = TRUE)]
  uniqueTU <- na.omit(uniqueTU)

  for (i in seq_along(uniqueTU)) {
    # select ID, position, delay, delay fragments, coordinates of the slope
    # and TU
    tu <- rowRanges(inp)[
      which(rowRanges(inp)$TU %in% uniqueTU[i]),
      c(
        "ID",
        "position",
        "delay",
        "delay_fragment",
        "slope",
        "intercept",
        "TU"
      )
    ]
    
    if (unique(strand(tu)) == "-") {
      rows_tu <- order(tu[strand(tu) == "-",]$position,
                    decreasing = FALSE)
      tu[strand(tu) == "-",] <- tu[strand(tu) == "-",][rows_tu,]
    }
    # select delay segments from TU
    del_segs <-
      unique(tu[grep(paste0("\\D_\\d+", "$"), tu$delay_fragment)]$delay_fragment)

    if (length(del_segs) > 1) {
      for (j in seq_len(length(del_segs) - 1)) {
        if (unique(strand(tu)) == "+") {
          del.1 <- tu[which(tu$delay_fragment == del_segs[j]), ]
          y1 <- last(del.1$position) * unique(del.1$slope) +
            unique(del.1$intercept)
          del.2 <- tu[which(tu$delay_fragment == del_segs[j + 1]), ]
          y2 <- last(del.1$position) * unique(del.2$slope) +
            unique(del.2$intercept)
        } else {
          del.1 <- tu[which(tu$delay_fragment == del_segs[j + 1]), ]
          rows <- length(del.1):1
          del.1 <- del.1[rows,]
          del.2 <- tu[which(tu$delay_fragment == del_segs[j]), ]
          rows <- length(del.2):1
          del.2 <- del.2[rows,]
          del <- c(del.1, del.2)
          del$position <- abs(del$position - max(del$position)) + 1
          del.1 <- del[seq_along(del.1), ]
          del.2 <- del[length(del.1) + seq_along(del.2), ]
          coef.del.1 <- coef(lm(del.1$delay ~ del.1$position))
          y1 <- last(del.1$position) * coef.del.1[2] + coef.del.1[1]
          coef.del.2 <- coef(lm(del.2$delay ~ del.2$position))
          y2 <- last(del.1$position) * coef.del.2[2] + coef.del.2[1]
        }
        dis <- abs(last(del.1$position) - del.2$position[1])
        if (dis > maxDis) {
          next()
        } else {
          y.dif <- y2 - y1
          if (y.dif >= 0) {
              rows <- match(last(del.1$ID), rowRanges(inp)$ID)
              rowRanges(inp)$pausing_site[rows] <- "+"
              rowRanges(inp)$event_duration[rows] <- abs(y.dif)
              rowRanges(inp)$ps_ts_fragment[rows] <-
                paste0(
                  del.1$delay_fragment[1],
                  ":",
                  del.2$delay_fragment[2]
                )
          } else if (y.dif < 0) {
              rows <- match(last(del.1$ID), rowRanges(inp)$ID)
              rowRanges(inp)$iTSS_I[rows] <- "+"
              rowRanges(inp)$event_duration[rows] <- y.dif
              rowRanges(inp)$ps_ts_fragment[rows] <-
                paste0(
                  del.1$delay_fragment[1],
                  ":",
                  del.2$delay_fragment[2]
                )
            }
          }
      }
    }
  }
  return(inp)
}
