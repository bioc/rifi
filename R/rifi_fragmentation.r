#' =========================================================================
#' rifi_fragmentation 
#' -------------------------------------------------------------------------
#' rifi_fragmentation wraps conveniently all fragmentation steps
#' 
#' rifi_fragmentation is wrapper of the following functions: 
#' 1. fragment_delay
 
#' 2. fragment_HL
 
#' 3. fragment_inty
 
#' 4. TUgether
 
#' 5. fragment_TI

#'
#' @param inp SummarizedExperiment: the input data frame with correct format.
#' @param cores integer: the number of assigned cores for the task.
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
#' rifi_fragmentation(inp = penalties_minimal, cores = 2)
#' 
#' @export


rifi_fragmentation <-
  function(inp,
           cores = 1,
           pen_delay = NULL,
           pen_out_delay = NULL,
           pen_HL = NULL,
           pen_out_HL = NULL,
           pen_inty = NULL,
           pen_out_inty = NULL,
           pen_TU = NULL,
           pen_TI = NULL,
           pen_out_TI = NULL) {
		
    lo <- metadata(inp)$logbook
    if(is.null(pen_delay)){
      pen_delay <- lo["delay_penalty"]
    }
    if(is.null(pen_out_delay)){
      pen_out_delay <- lo["delay_outlier_penalty"]
    }
    if(is.null(pen_HL)){
      pen_HL <- lo["half_life_penalty"]
    }
    if(is.null(pen_out_HL)){
      pen_out_HL <- lo["half_life_outlier_penalty"]
    }
    if(is.null(pen_inty)){
      pen_inty <- lo["intensity_penalty"]
    }
    if(is.null(pen_out_inty)){
      pen_delay <- lo["intensity_outlier_penalty"]
    }
    if(is.null(pen_TU)){
      pen_TU <- -0.75
    }
    if(is.null(pen_TI)){
      pen_TI <- lo["TI_penalty"]
    }
    if(is.null(pen_out_TI)){
      pen_out_TI <- lo["TI_outlier_penalty"]
    }
	
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
