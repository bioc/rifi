#' rifi_fragmentation: conveniently wraps all fragmentation steps
#'
#' rifi_fragmentation wraps the following functions:
#' 1. fragment_delay
#' 2. fragment_HL
#' 3. fragment_inty
#' 4. TUgether
#' 5. fragment_TI
#'
#' @param inp SummarizedExperiment: the input data frame with correct format.
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
#' @return the SummarizedExperiment object: with delay_fragment, HL_fragment,
#' intensity_fragment, TI_termination_fragment and TU, and the respective values
#' added to the rowRanges.
#' 
#' @seealso `fragment_delay`
#' @seealso `fragment_HL`
#' @seealso `fragment_inty`
#' @seealso `TUgether`
#' @seealso `fragment_TI`
#' 
#' @examples
#' data(penalties_minimal)
#' rifi_fragmentation(inp = penalties_minimal, cores = 2,
#' logbook = metadata(penalties_minimal)$logbook)
#' 
#' @export


rifi_fragmentation <-
  function(inp,
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
    
    message("running fragment_delay...")
    inp <-
      fragment_delay(
        inp = inp,
        cores = cores,
        pen = pen_delay,
        pen_out = pen_out_delay
      )

    message("running fragment_HL...")
    inp <-
      fragment_HL(
        inp = inp,
        cores = cores,
        pen = pen_HL,
        pen_out = pen_out_HL
      )

    message("running fragment_inty...")
    inp <-
      fragment_inty(
        inp = inp,
        cores = cores,
        pen = pen_inty,
        pen_out = pen_out_inty
      )

    message("running TUgether...")
    inp <- TUgether(
      inp = inp,
      cores = cores,
      pen = pen_TU
    )

    message("running fragment_TI...")
    inp <-
      fragment_TI(
        inp = inp,
        cores = cores,
        pen = pen_TI,
        pen_out = pen_out_TI
      )

    rowRanges(inp)$seg_ID <-
      paste(
        rowRanges(inp)$position_segment,
        rowRanges(inp)$TU,
        rowRanges(inp)$delay_fragment,
        rowRanges(inp)$HL_fragment,
        rowRanges(inp)$intensity_fragment,
        sep = "|"
      )
    rowRanges(inp)$seg_ID <- gsub("_O", "", rowRanges(inp)$seg_ID)

    inp
  }
