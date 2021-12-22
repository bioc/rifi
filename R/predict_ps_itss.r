#' predict_ps_itss: predicts pausing sites (ps) and internal starting sites
#' (ITSS) between delay fragments.
#' predict_ps_itss predicts ps and ITSS within the same TU. Neighboring delay
#' segments are compared to
#' each other by positioning the intercept of the second segment into the first
#' segment using slope and intercept coefficients.#'
#' predict_ps_itss uses 3 steps to identify ps and ITSS:
#' 1. select unique TU.
#' 2. select from the input dataframe the columns: ID, position, strand, delay,
#' delay fragment, TU and slope coordinates, velocity_fragment and intercept.
#' 3. select delay segments in the TU.
#' 4. loop into all delay segments and estimate the coordinates of the last
#' point of the first segment using the coefficients
#' of the second segment and vice versa. We get two predicted positions, the
#' difference between them is compared to the threshold.
#' In case the strand is "-", additional steps are added:
#' The positions of both segments are ordered from the last position to the
#' first one.
#' all positions are merged in one column and subtracted from the maximum
#' position.
#' the column is split in 2. The first and second correspond to the positions
#' of the first and second segments respectively.
#' Both segments are subjected to lm fit and the positions predicted are used
#' on the same way as the opposite strand.
#' if the difference between the positions predicted is lower than negative
#' threshold, ps is assigned otherwise, and if the difference is higher than
#' the positive threshold, ITSS is assigned.
#'
#' @param data dataframe: the probe based data frame.
#' @param maxDis integer: the maximal distance allowed between two successive
#' fragments.
#' 
#' @return the probe data frame with the columns regarding statistics:
#' \describe{
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{position:}{The bin/probe specific position}
#'   \item{strand:}{The bin/probe specific strand}
#'   \item{delay:}{The delay value of the bin/probe}
#'   \item{delay_fragment:}{The delay fragment the bin belongs to}
#'   \item{intercept:}{The vintercept of fit through the respective delay
#'   fragment}
#'   \item{slope:}{The slope of the fit through the respective delay fragment}
#'   \item{TU:}{The overarching transcription unit}
#'   \item{pausing_site:}{}
#'   \item{iTSS_I:}{}
#'   \item{ps_ts_fragment:}{}
#' }
#' 
#' @examples
#' data(fit_minimal)
#' predict_ps_itss(data = fit_minimal, maxDis = 300)
#' 
#' @export

predict_ps_itss <- function(data, maxDis = 300) {
  data$pausing_site <- "-"
  data$iTSS_I <- "-"
  data$ps_ts_fragment <- NA

  # select unique TUs
  uniqueTU <- unique(data$TU)
  uniqueTU <- uniqueTU[grep("_NA|_T", uniqueTU, invert = TRUE)]
  uniqueTU <- na.omit(uniqueTU)

  for (i in seq_along(uniqueTU)) {
    # select ID, position, delay, delay fragments, coordinates of the slope
    # and TU
    tu <- data[
      which(data$TU %in% uniqueTU[i]),
      c(
        "ID",
        "position",
        "strand",
        "delay",
        "delay_fragment",
        "slope",
        "intercept",
        "TU"
      )
    ]
    if (unique(tu$strand) == "-") {
      tu <- tu[order(tu$position, decreasing = FALSE), ]
    }

    # select delay segments in the TU
    del_segs <-
      unique(tu[grep(paste0("\\D_\\d+", "$"), tu$delay_fragment),
                "delay_fragment"])

    if (length(del_segs) > 1) {
      for (j in seq_len(length(del_segs) - 1)) {
        if (unique(tu$strand) == "+") {
          del.1 <- tu[which(tu$delay_fragment == del_segs[j]), ]
          y1 <- last(del.1$position) * unique(del.1$slope) +
            unique(del.1$intercept)
          del.2 <- tu[which(tu$delay_fragment == del_segs[j + 1]), ]
          y2 <- last(del.1$position) * unique(del.2$slope) +
            unique(del.2$intercept)
        } else {
          del.1 <- tu[which(tu$delay_fragment == del_segs[j + 1]), ]
          del.1 <- del.1[nrow(del.1):1, ]
          del.2 <- tu[which(tu$delay_fragment == del_segs[j]), ]
          del.2 <- del.2[nrow(del.2):1, ]
          del <- rbind(del.1, del.2)
          del$position <- abs(del$position - max(del$position)) + 1
          del.1 <- del[seq_len(nrow(del.1)), ]
          del.2 <- del[nrow(del.1) + seq_len(nrow(del.2)), ]
          coef.del.1 <- coef(lm(del.1$delay ~ del.1$position))
          y1 <- last(del.1$position) * coef.del.1[2] + coef.del.1[1]
          coef.del.2 <- coef(lm(del.2$delay ~ del.2$position))
          y2 <- last(del.1$position) * coef.del.2[2] + coef.del.2[1]
        }
        dis <- abs(last(del.1$position) - del.2$position[1])
        if (dis > maxDis) {
          next()
        } else {
          y.dif <- y1 - y2
          if (y.dif <= 0) {
            if (unique(tu$strand) == "+") {
              data[which(data$ID %in% last(del.1$ID)), "pausing_site"] <- "+"
              data[which(data$delay_fragment %in% del.1$delay_fragment[1]),
                   "ps_ts_fragment"] <-
                paste0(
                  del.1$delay_fragment[1],
                  ":",
                  del.2$delay_fragment[2]
                )
            } else {
              data[which(data$ID %in% del.2$ID[1]), "pausing_site"] <- "+"
              data[which(data$delay_fragment %in% del.2$delay_fragment[1]),
                   "ps_ts_fragment"] <-
                paste0(
                  del.2$delay_fragment[1],
                  ":",
                  del.1$delay_fragment[2]
                )
            }
          } else if (y.dif > 0) {
            if (unique(tu$strand) == "+") {
              data[which(data$ID %in% last(del.1$ID)), "iTSS_I"] <- "+"
              data[which(data$delay_fragment %in% del.1$delay_fragment[1]),
                   "ps_ts_fragment"] <-
                paste0(
                  del.1$delay_fragment[1],
                  ":",
                  del.2$delay_fragment[2]
                )
            } else {
              data[which(data$ID %in% del.2$ID[1]), "iTSS_I"] <- "+"
              data[which(data$delay_fragment %in% del.2$delay_fragment[1]),
                   "ps_ts_fragment"] <-
                paste0(
                  del.2$delay_fragment[1],
                  ":",
                  del.1$delay_fragment[2]
                )
            }
          }
        }
      }
    }
  }
  return(data)
}
