#' rifi_penalties: conveniently wraps all penalty steps
#'
#' wraps the functions: make_pen and viz_pen_obj.
#'
#' @param probe data frame: the probe based data frame.
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
#' @return If details is set to TRUE a list with 5 items otherwise a vector of
#' length 8:
#' \describe{
#'   \item{logbook:}{The logbook vector containing all penalty information}
#'   \item{pen_obj_delay:}{A list with 4 items:
#'     \describe{
#'       \item{logbook:}{The logbook vector containing all penalty information}
#'       \item{delay_penalties:}{a vetor with the delay penalty and delay
#'       outlier penalty}
#'       \item{correct:}{a matrix of the correct splits}
#'       \item{wrong:}{a matrix of the incorrect splits}
#'     }
#'   }
#'   \item{pen_obj_HL:}{A list with 4 items:
#'     \describe{
#'       \item{logbook:}{The logbook vector containing all penalty information}
#'       \item{HL_penalties:}{a vetor with the half-life penalty and half-life
#'       outlier penalty}
#'       \item{correct:}{a matrix of the correct splits}
#'       \item{wrong:}{a matrix of the incorrect splits}
#'     }
#'   }
#'   \item{pen_obj_inty:}{A list with 4 items:
#'     \describe{
#'       \item{logbook:}{The logbook vector containing all penalty information}
#'       \item{inty_penalties:}{a vetor with the intensity penalty and intensity
#'       outlier penalty}
#'       \item{correct:}{a matrix of the correct splits}
#'       \item{wrong:}{a matrix of the incorrect splits}
#'     }
#'   }
#'   \item{pen_obj_TI:}{A list with 4 items:
#'     \describe{
#'       \item{logbook:}{The logbook vector containing all penalty information}
#'       \item{TI_penalties:}{a vetor with the TI penalty and TI outlier
#'       penalty}
#'       \item{correct:}{a matrix of the correct splits}
#'       \item{wrong:}{a matrix of the incorrect splits}
#'     }
#'   }
#' }
#' 
#'
#' @seealso `make_pen`
#' @seealso `viz_pen_obj`
#'
#' @examples
#' data(fit_minimal)
#' rifi_penalties(
#'   probe = fit_minimal, details = FALSE, viz = FALSE,
#'   top_i = 25, cores = 2, dpt = 1, smpl_min = 10, smpl_max = 100,
#'   sta_pen = 0.5,
#'   end_pen = 4.5, rez_pen = 9, sta_pen_out = 0.5, end_pen_out = 4.5,
#'   rez_pen_out = 9
#' )
#' 
#' @export

rifi_penalties <- function(probe,
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
  num_args <-
    list(
      top_i,
      cores,
      dpt,
      smpl_min,
      smpl_max,
      sta_pen,
      end_pen,
      rez_pen,
      sta_pen_out,
      end_pen_out,
      rez_pen_out
    )
  names(num_args) <-
    c(
      "top_i",
      "cores",
      "dpt",
      "smpl_min",
      "smpl_max",
      "sta_pen",
      "end_pen",
      "rez_pen",
      "sta_pen_out",
      "end_pen_out",
      "rez_pen_out"
    )
  assert(
    all(unlist(lapply(
      num_args,
      FUN = function(x) {
        (is.numeric(x) &
          length(x) == 1)
      }
    ))),
    paste0("'", names(which(
      unlist(lapply(
        num_args,
        FUN = function(x) {
          (is.numeric(x) &
            length(x) == 1)
        }
      )) == FALSE
    ))[1], "' must be numeric of length one")
  )
  assert(cores > 0, "'cores' must be a positive integer")
  assert(is.logical(details), "'details' must be a logical")
  assert(is.logical(viz), "'viz' must be a logical")
  req_cols_probe <-
    c(
      "ID",
      "position",
      "strand",
      "intensity",
      "position_segment",
      "delay",
      "half_life",
      "TI_termination_factor"
    )
  assert(
    all(req_cols_probe %in% colnames(probe)),
    paste0("'", req_cols_probe[which(!req_cols_probe %in% colnames(probe))],
           "' must be a column in 'probe'!")
  )

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
      probe = probe,
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
      probe = probe,
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
      probe = probe,
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
      probe = probe,
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

  tmp_df <-
    data.frame(
      ID = probe$ID,
      delay = probe$delay,
      half_life = probe$half_life,
      intensity = probe$intensity,
      TI = probe$TI_termination_factor,
      seg = probe$position_segment,
      flag = probe$flag
    )

  spltd <- split(tmp_df, tmp_df$seg)

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

  res <- logbook

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
  }

  res
}
