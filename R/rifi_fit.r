#' rifi_fit: conveniently wraps all fitting steps
#'
#' Wraps the functions: nls2_fit, TI_fit, plot_nls2_function and
#' plot_singleProbe_function.
#'
#' @param inp SummarizedExperiment: the input with correct format.
#' @param cores integer: the number of assigned cores for the task.
#' @param viz logical: whether to visualize the output.
#' @param restr numeric: a parameter that restricts the freedom of the fit to
#' avoid wrong TI-term_factors, ranges from 0 to 0.2
#' @param decay numeric vector: A sequence of starting values for the decay.
#' Default is seq(.08, 0.11, by=.02)
#' @param delay numeric vector: A sequence of starting values for the delay.
#' Default is seq(0,10, by=.1)
#' @param k numeric vector: A sequence of starting values for the synthesis
#' rate.Default is seq(0.1,1,0.2)
#' @param bg numeric vector: A sequence of starting values. Default is 0.2.
#' @param TI_k numeric vector: A sequence of starting values for the synthesis
#' rate. Default is seq(0, 1, by = 0.5).
#' @param TI_decay numeric vector: A sequence of starting values for the decay.
#' Default is c(0.05, 0.1, 0.2, 0.5, 0.6).
#' @param TI numeric vector: A sequence of starting values for the TI. Default
#' is  seq(0, 1, by = 0.5).
#' @param TI_delay numeric vector: A sequence of starting values for the delay.
#' Default is seq(0, 2, by = 0.5).
#' @param TI_rest_delay numeric vector: A sequence of starting values. Default
#' is seq(0, 2, by = 0.5).
#' @param TI_bg numeric vector: A sequence of starting values. Default is 0.
#' 
#' @return the SummarizedExperiment object: with delay, decay  and
#' TI_termination_factor added to the rowRanges. The full fit data is saved in
#' the metadata as "fit_STD" and "fit_TI". A plot is given if viz = TRUE.
#'
#' @seealso `nls2_fit`
#' @seealso `TI_fit`
#' @seealso `plot_nls2`
#' @seealso `plot_singleProbe`
#'
#' @examples
#' data(preprocess_minimal)
#' rifi_fit(
#'   inp = preprocess_minimal,
#'   cores = 1, viz = FALSE, restr = 0.1,
#'   decay = seq(.08, 0.11, by = .02),
#'   delay = seq(0, 10, by = .1), k = seq(0.1, 1, 0.2), bg = 0.2,
#'   TI_k = seq(0, 1, by = 0.5), TI_decay = c(0.05, 0.1, 0.2, 0.5, 0.6),
#'   TI = seq(0, 1, by = 0.5), TI_delay = seq(0, 2, by = 0.5),
#'   TI_rest_delay = seq(0, 2, by = 0.5), TI_bg = 0,
#' )
#' 
#' @export

rifi_fit <-
  function(inp,
           cores = 1,
           viz = FALSE,
           restr = 0.2,
           decay = seq(.08, 0.11, by = .02),
           delay = seq(0, 10, by = .1),
           k = seq(0.1, 1, 0.2),
           bg = 0.2,
           TI_k = seq(0, 1, by = 0.5),
           TI_decay = c(0.05, 0.1, 0.2, 0.5, 0.6),
           TI = seq(0, 1, by = 0.5),
           TI_delay = seq(0, 2, by = 0.5),
           TI_rest_delay = seq(0, 2, by = 0.5),
           TI_bg = 0) {
    
    message("running nls2_fit...")
    inp <- tryCatch({
      inp <-
        nls2_fit(
          inp = inp,
          cores = cores,
          decay = decay,
          delay = delay,
          k = k,
          bg = bg
        )
      inp
    },
    error = function(e) {
      writeLines(
        paste(
          "An unknown error has appeared!\n",
          e,
          "An emergency output was returned!\n Please rerun nls2_fit manually!"
        )
      )
      return(inp)
    })
    
    message("running TI_fit...")
    inp <- tryCatch({
      inp <-
        TI_fit(
          inp = inp,
          cores = cores,
          restr = restr,
          k = TI_k,
          decay = TI_decay,
          ti = TI,
          ti_delay = TI_delay,
          rest_delay = TI_rest_delay,
          bg = TI_bg
        )
      inp
    },
    error = function(e) {
      writeLines(
        paste(
          "An unknown error has appeared!\n",
          e,
          "An emergency output was returned!\n Please rerun TI_fit manually!"
        )
      )
      return(inp)
    })
    
    tryCatch({
      if (viz == TRUE) {
        message("running visualization...")
        plot_nls2_function(inp = inp)
      }
    },
    error = function(e) {
      writeLines(
        paste(
          "An unknown error has appeared!\n",
          e,
          "The visualization could not be plotted!\n Please rerun
          plot_nls2_function or plot_singleProbe_function manually!"
        )
      )
    })
    inp
  }
