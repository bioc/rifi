#' rifi_fragmentation: conveniently wraps all fragmentation steps
#'
#' rifi_fragmentation wraps the following functions:
#' 1. fragment_delay
#' 2. fragment_HL
#' 3. fragment_inty
#' 4. TUgether
#' 5. fragment_TI
#'
#' @param probe probe data frame: the probe based data frame.
#' @param cores integer: the number of assigned cores for the task.
#' @param logbook numeric vector: the logbook vector, if it exists.
#' @param pen_delay numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer fragments. Default is the auto generated value.
#' @param pen_out_delay numeric: an internal parameter for the dynamic
#' programming. Higher values result in fewer allowed outliers. Default is the
#' auto generated value.
#' @param pen_HL numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer fragments. Default is the auto generated value.
#' @param pen_out_HL numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer allowed outliers. Default is the auto generated
#' value.
#' @param pen_inty numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer fragments. Default is the auto generated value.
#' @param pen_out_inty numeric: an internal parameter for the dynamic
#' programming. Higher values result in fewer allowed outliers. Default is
#' the auto generated value.
#' @param pen_TU numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer fragments. Default -0.75.
#' @param pen_TI numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer fragments. Default is the auto generated value.
#' @param pen_out_TI numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer allowed outliers. Default is the auto generated
#' value.
#'
#' @return the probe data frame with the columns regarding the fragmentation:
#' \describe{
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{position:}{The bin/probe specific position}
#'   \item{strand:}{The bin/probe specific strand}
#'   \item{intensity:}{The relative intensity at time point 0}
#'   \item{probe_TI:}{An internal value to determine which fitting model is
#'   applied}
#'   \item{flag:}{Information on which fitting model is applied}
#'   \item{position_segment:}{The position based segment}
#'   \item{delay:}{The delay value of the bin/probe}
#'   \item{half_life:}{The half-life of the bin/probe}
#'   \item{TI_termination_factor:}{The termination factor of the bin/probe}
#'   \item{delay_fragment:}{The delay fragment the bin belongs to}
#'   \item{velocity_fragment:}{The velocity value of the respective delay
#'   fragment}
#'   \item{intercept:}{The vintercept of fit through the respective delay
#'   fragment}
#'   \item{slope:}{The slope of the fit through the respective delay fragment}
#'   \item{HL_fragment:}{The half-life fragment the bin belongs to}
#'   \item{HL_mean_fragment:}{The mean half-life value of the respective
#'   half-life fragment}
#'   \item{intensity_fragment:}{The intensity fragment the bin belongs to}
#'   \item{intensity_mean_fragment:}{The mean intensity value of the respective
#'   intensity fragment}
#'   \item{TU:}{The overarching transcription unit}
#'   \item{TI_termination_fragment:}{The TI fragment the bin belongs to}
#'   \item{TI_mean_termination_factor:}{The mean termination factor of the
#'   respective TI fragment}
#'   \item{seg_ID:}{The combined ID of the fragment}
#' }
#' 
#' @seealso `fragment_delay`
#' @seealso `fragment_HL`
#' @seealso `fragment_inty`
#' @seealso `TUgether`
#' @seealso `fragment_TI`
#' 
#' @examples
#' data(fit_minimal)
#' data(penalties_minimal)
#' rifi_fragmentation(probe = fit_minimal, cores = 2,
#' logbook = penalties_minimal)
#' 
#' @export


rifi_fragmentation <-
  function(probe,
           cores = 1,
           logbook,
           pen_delay = logbook["delay_penalty"],
           pen_out_delay = logbook["delay_outlier_penalty"],
           pen_HL = logbook["half_life_penalty"],
           pen_out_HL = logbook["half_life_outlier_penalty"],
           pen_inty = logbook["intensity_penalty"],
           pen_out_inty = logbook["intensity_outlier_penalty"],
           pen_TU = -0.75,
           pen_TI = logbook["TI_penalty"],
           pen_out_TI = logbook["TI_outlier_penalty"]) {
    num_args <-
      list(
        cores,
        pen_delay,
        pen_out_delay,
        pen_HL,
        pen_out_HL,
        pen_inty,
        pen_out_inty,
        pen_TU,
        pen_TI,
        pen_out_TI
      )
    names(num_args) <-
      c(
        "cores",
        "pen_delay",
        "pen_out_delay",
        "pen_HL",
        "pen_out_HL",
        "pen_inty",
        "pen_out_inty",
        "pen_TU",
        "pen_TI",
        "pen_out_TI"
      )
    assert(
      all(unlist(lapply(
        num_args,
        FUN = function(x) {
          (is.numeric(x) &
            length(x) == 1)
        }
      ))),
      paste0(
        "'",
        names(which(unlist(
          lapply(
            num_args,
            FUN = function(x) {
              (is.numeric(x) &
                length(x) == 1)
            }
          )
        ) == FALSE))[1],
        "' must be numeric of length one or given by the logbook"
      )
    )
    assert(cores > 0, "'cores' must be a positive integer")
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

    message("running fragment_delay...")
    probe <-
      fragment_delay(
        probe = probe,
        cores = cores,
        pen = pen_delay,
        pen_out = pen_out_delay
      )

    message("running fragment_HL...")
    probe <-
      fragment_HL(
        probe = probe,
        cores = cores,
        pen = pen_HL,
        pen_out = pen_out_HL
      )

    message("running fragment_inty...")
    probe <-
      fragment_inty(
        probe = probe,
        cores = cores,
        pen = pen_inty,
        pen_out = pen_out_inty
      )

    message("running TUgether...")
    probe <- TUgether(
      probe = probe,
      cores = cores,
      pen = pen_TU
    )

    message("running fragment_TI...")
    probe <-
      fragment_TI(
        probe = probe,
        cores = cores,
        pen = pen_TI,
        pen_out = pen_out_TI
      )

    probe$seg_ID <-
      paste(
        probe$position_segment,
        probe$TU,
        probe$delay_fragment,
        probe$HL_fragment,
        probe$intensity_fragment,
        sep = "|"
      )
    probe$seg_ID <- gsub("_O", "", probe$seg_ID)

    probe
  }
