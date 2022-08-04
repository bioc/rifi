#' =========================================================================
#' rifi_penalties          
#' -------------------------------------------------------------------------
#' rifi_penalties wraps conveniently all penalty steps
#'
#' rifi_penalties wraps the functions:

#' 1. make_pen,

#' 2. viz_pen_obj
#'
#' @param inp SummarizedExperiment: the input data frame with correct format.
#' @param details logical: whether to return the penalty objects or just the
#' logbook.
#' @param viz logical: whether to visualize the output or not. Default is FALSE
#' @param top_i integer: the number of top results visualized. Default is all.
#' @param cores integer: the number of assigned cores for the task.
#' @param dpt integer: the number of times a full iteration cycle is repeated
#' with a more narrow range based on the previous cycle. Default is 2.
#' @param smpl_min integer: the smaller end of the sampling size. Default is 10.
#' @param smpl_max integer: the larger end of the sampling size. Default is 100.
#' @param sta_pen numeric: the lower starting penalty. Default is 0.5.
#' @param end_pen numeric: the higher starting penalty. Default is 4.5.
#' @param rez_pen numeric: the number of penalties iterated within the penalty
#' range. Default is 9.
#' @param sta_pen_out numeric: the lower starting outlier penalty. Default is
#' 0.5.
#' @param end_pen_out numeric: the higher starting outlier penalty. Default is
#' 3.5.
#' @param rez_pen_out numeric: the number of outlier penalties iterated within
#' the outlier penalty range. Default is 7.
#' 
#' @return The SummarizedExperiment object: with the penalties in the
#' logbook added to the metadata. Also adds logbook_details if details is TRUE,
#' and plots the penalties if viz is TRUE.
#'
#' @seealso `make_pen`
#' @seealso `viz_pen_obj`
#'
#' @examples
#' data(fit_minimal)
#' rifi_penalties(
#'   inp = fit_minimal, details = FALSE, viz = FALSE,
#'   top_i = 25, cores = 2, dpt = 1, smpl_min = 10, smpl_max = 100,
#'   sta_pen = 0.5, end_pen = 4.5, rez_pen = 9, sta_pen_out = 0.5,
#'   end_pen_out = 4.5, rez_pen_out = 9
#' )
#' 
#' @export

rifi_penalties <- function(inp,
                           details = FALSE,
                           viz = FALSE,
                           top_i = 25,
                           cores = 1,
                           dpt = 1,
                           smpl_min = 10,
                           smpl_max = 100,
                           sta_pen = 0.5,
                           end_pen = 4.5,
                           rez_pen = 9,
                           sta_pen_out = 0.5,
                           end_pen_out = 4.5,
                           rez_pen_out = 9) {
  
  logbook <- as.numeric(rep(NA, 8))

  names(logbook) <- c(
    "delay_penalty",
    "delay_outlier_penalty",
    "half_life_penalty",
    "half_life_outlier_penalty",
    "intensity_penalty",
    "intensity_outlier_penalty",
    "TI_penalty",
    "TI_outlier_penalty"
  )

  message("running make_pen on delay...")
  pen_obj_delay <-
    make_pen(
      inp = inp,
      FUN = fragment_delay_pen,
      cores = cores,
      logs = logbook,
      dpt = dpt,
      smpl_min = smpl_min,
      smpl_max = smpl_max,
      sta_pen = sta_pen,
      end_pen = end_pen,
      rez_pen = rez_pen,
      sta_pen_out = sta_pen_out,
      end_pen_out = end_pen_out,
      rez_pen_out = rez_pen_out
    )
  logbook <- pen_obj_delay[[1]]

  message("running make_pen on half-life...")
  pen_obj_HL <-
    make_pen(
      inp = inp,
      FUN = fragment_HL_pen,
      cores = cores,
      logs = logbook,
      dpt = dpt,
      smpl_min = smpl_min,
      smpl_max = smpl_max,
      sta_pen = sta_pen,
      end_pen = end_pen,
      rez_pen = rez_pen,
      sta_pen_out = sta_pen_out,
      end_pen_out = end_pen_out,
      rez_pen_out = rez_pen_out
    )
  logbook <- pen_obj_HL[[1]]

  message("running make_pen on intensity...")
  pen_obj_inty <-
    make_pen(
      inp = inp,
      FUN = fragment_inty_pen,
      cores = cores,
      logs = logbook,
      dpt = dpt,
      smpl_min = smpl_min,
      smpl_max = smpl_max,
      sta_pen = sta_pen,
      end_pen = end_pen,
      rez_pen = rez_pen,
      sta_pen_out = sta_pen_out,
      end_pen_out = end_pen_out,
      rez_pen_out = rez_pen_out
    )
  logbook <- pen_obj_inty[[1]]

  message("running make_pen on TI...")
  pen_obj_TI <-
    make_pen(
      inp = inp,
      FUN = fragment_TI_pen,
      cores = cores,
      logs = logbook,
      dpt = dpt,
      smpl_min = smpl_min,
      smpl_max = smpl_max,
      sta_pen = sta_pen,
      end_pen = end_pen,
      rez_pen = rez_pen,
      sta_pen_out = sta_pen_out,
      end_pen_out = end_pen_out,
      rez_pen_out = rez_pen_out
    )
  logbook <- pen_obj_TI[[1]]

  if (viz == TRUE) {
    message("running visualization...")
    viz_pen_obj(obj = pen_obj_delay, top_i = top_i)
    viz_pen_obj(obj = pen_obj_HL, top_i = top_i)
    viz_pen_obj(obj = pen_obj_inty, top_i = top_i)
    viz_pen_obj(obj = pen_obj_TI, top_i = top_i)
  }
  
  tmp_df <- inp_df(inp, "ID", "delay", "half_life", "intensity",
                   "TI_termination_factor", "position_segment", "flag")
  
  spltd <- split(tmp_df, tmp_df$position_segment)

  A <- lapply(spltd, function(x) {
    x <- na.omit(x$delay)
  })

  B <- lapply(spltd, function(x) {
    x <- na.omit(x$half_life)
  })

  C <- lapply(spltd, function(x) {
    x <- na.omit(x$intensity)
  })

  D <- lapply(spltd, function(x) {
    x <- na.omit(x$TI)
  })

  if (!any(between(lapply(A, length), smpl_min, smpl_max))) {
    warning(
      "There is no position segment with enough delay values in the given
      sample range! Default penalties for delay fragmentation will be returned!"
    )
    logbook[c("delay_penalty", "delay_outlier_penalty")] <- c(1.8, 0.7)
  }

  if (!any(between(lapply(B, length), smpl_min, smpl_max))) {
    warning(
      "There is no position segment with enough half_life values in the given
      sample range! Default penalties for half_life fragmentation will be
      returned!"
    )
    logbook[c("half_life_penalty", "half_life_outlier_penalty")] <-
      c(1.5, 1)
  }

  if (!any(between(lapply(C, length), smpl_min, smpl_max))) {
    warning(
      "There is no position segment with enough intensity values in the given
      sample range! Default penalties for intensity fragmentation will be
      returned!"
    )
    logbook[c("intensity_penalty", "intensity_outlier_penalty")] <-
      c(2, 1)
  }

  if (!any(between(lapply(D, length), smpl_min, smpl_max))) {
    warning(
      "There is no position segment with enough TI_termination values in the
      given sample range! Default penalties for TI fragmentation will be
      returned!"
    )
    logbook[c("TI_penalty", "TI_outlier_penalty")] <- c(1.5, 1)
  }

  if (details == TRUE) {
    res <-
      list(
        logbook,
        pen_obj_delay,
        pen_obj_HL,
        pen_obj_inty,
        pen_obj_TI
      )
    names(res) <-
      c(
        "logbook",
        "pen_obj_delay",
        "pen_obj_HL",
        "pen_obj_inty",
        "pen_obj_TI"
      )
    metadata(inp)$logbook_details <- res
  }
  
  metadata(inp)$logbook <- logbook
  
  inp
}
