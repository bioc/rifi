#' rifi_fit: conveniently wraps all fitting steps
#'
#' Wraps the functions: nls2_fit, TI_fit, plot_nls2 and
#' plot_singleProbe_function.
#'
#' @param inp data frame: the input data frame with correct format.
#' @param probe data frame: the probe based data frame.
#' @param cores integer: the number of assigned cores for the task.
#' @param details logical: whether to return the fit objects or just the probe.
#' @param viz logical: whether to visualize the output.
#' @param restr numeric: a parameter that restricts the freedom of the fit to
#' avoid wrong TI-term_factors, ranges from 0 to 0.2
#' @param decay numeric vector: A sequence of starting values for the decay.
#' Default is seq(.08, 0.11, by=.02)
#' @param delay numeric vector: A sequence of starting values for the delay.
#' Default is seq(0,10, by=.1)
#' @param k numeric vector: A sequence of starting values for the synthesis
#' rate.Default is seq(0.1,1,0.2)
#' @param intyf numeric vector: A sequence of starting values. Default is 0.2.
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
#' @param color string vector: A sequence of color for plot.
#' 
#' @return A list of 3 data frames :
#' \describe{
#'   \item{probe_df:}{the probe data frame with the columns delay, half_life
#'   and TI_termination_factor added:
#'   \describe{
#'     \item{ID:}{The bin/probe specific ID}
#'     \item{position:}{The bin/probe specific position}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{intensity:}{The relative intensity at time point 0}
#'     \item{probe_TI:}{An internal value to determine which fitting model is
#'     applied}
#'     \item{flag:}{Information on which fitting model is applied}
#'     \item{postion_segment:}{The position based segment}
#'     \item{delay:}{The delay value of the bin/probe}
#'     \item{half_life:}{The half-life of the bin/probe}
#'     \item{TI_termination_factor:}{The termination factor of the bin/probe}
#'     }
#'   }
#'   \item{fit_obj_STD:}{the fit object for the standard fit:
#'   \describe{
#'     \item{ID:}{The bin/probe specific ID}
#'     \item{delay:}{The delay value of the bin/probe}
#'     \item{half_life:}{The half-life of the bin/probe}
#'     \item{inty_S0:}{The relative intensity at time point 0}
#'     \item{intyf:}{The background value of the fit}
#'     }
#'   }
#'   \item{fit_obj_TI:}{the fit object for the TI fit:
#'   \describe{
#'     \item{delay:}{The delay value of the bin/probe}
#'     \item{ti_delay:}{The ti-delay value of the bin/probe}
#'     \item{half_life:}{The half-life of the bin/probe}
#'     \item{ti_value:}{The ti-value of the bin/probe}
#'     \item{TI_termination_factor:}{The termination factor of the bin/probe}
#'     \item{synthesis_rate:}{The synthesis rate of the bin/probe}
#'     \item{TI_background:}{The background value of the fit}
#'     \item{position:}{The bin/probe specific position}
#'     \item{ID:}{The bin/probe specific ID}
#'     }
#'   }
#' }
#'
#' @seealso `nls2_fit`
#' @seealso `TI_fit`
#' @seealso `plot_nls2`
#' @seealso `plot_singleProbe`
#'
#' @examples
#' data(preprocess_minimal)
#' rifi_fit(
#'   inp = preprocess_minimal$input_df, probe = preprocess_minimal$probe_df,
#'   cores = 2, viz = FALSE, details = FALSE, restr = 0.2,
#'   decay = seq(.08, 0.11, by = .02),
#'   delay = seq(0, 10, by = .1), k = seq(0.1, 1, 0.2), intyf = 0.2,
#'   TI_k = seq(0, 1, by = 0.5), TI_decay = c(0.05, 0.1, 0.2, 0.5, 0.6),
#'   TI = seq(0, 1, by = 0.5), TI_delay = seq(0, 2, by = 0.5),
#'   TI_rest_delay = seq(0, 2, by = 0.5)
#' )
#' 
#' @export

rifi_fit <-
  function(inp,
           probe,
           cores = 1,
           details = FALSE,
           viz = FALSE,
           restr = 0.2,
           decay = seq(.08, 0.11, by = .02),
           delay = seq(0, 10, by = .1),
           k = seq(0.1, 1, 0.2),
           intyf = 0.2,
           TI_k = seq(0, 1, by = 0.5),
           TI_decay = c(0.05, 0.1, 0.2, 0.5, 0.6),
           TI = seq(0, 1, by = 0.5),
           TI_delay = seq(0, 2, by = 0.5),
           TI_rest_delay = seq(0, 2, by = 0.5),
           color = c(
             "blue",
             "green",
             "yellow",
             "grey",
             "orange",
             "cyan",
             "pink",
             "grey64",
             "grey90",
             "black",
             "grey45"
           )) {
    num_args <- list(cores, restr)
    names(num_args) <- c("cores", "restr")
    assert(all(unlist(lapply(
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
    ))[1], "' must be numeric of length one"))
    vec_args <-
      list(decay,
           delay,
           k,
           intyf,
           TI_k,
           TI_decay,
           TI,
           TI_delay,
           TI_rest_delay)
    names(vec_args) <-
      c(
        "decay",
        "delay",
        "k",
        "intyf",
        "TI_k",
        "TI_decay",
        "TI",
        "TI_delay",
        "TI_rest_delay"
      )
    assert(all(unlist(lapply(
      vec_args,
      FUN = function(x) {
        (is.numeric(x) &
           is.vector(x))
      }
    ))), paste0("'", names(which(
      unlist(lapply(
        vec_args,
        FUN = function(x) {
          (is.numeric(x) &
             is.vector(x))
        }
      )) == FALSE
    ))[1], "' must be numeric vector"))
    assert(is.logical(details), "'details' must be a logical")
    assert(is.logical(viz), "'viz' must be a logical")
    assert(cores > 0, "'cores' must be a positive integer")
    req_cols_inp <- c("0", "ID", "position", "strand", "filtration")
    assert(
      all(req_cols_inp %in% colnames(inp)),
      paste0("'", req_cols_inp[which(!req_cols_inp %in% colnames(inp))],
             "' must be a column in 'inp'!")
    )
    req_cols_probe <- c("ID", "position", "strand", "flag")
    assert(
      all(req_cols_probe %in% colnames(probe)),
      paste0("'", req_cols_probe[which(!req_cols_probe %in% colnames(probe))],
             "' must be a column in 'probe'!")
    )
    
    message("running nls2_fit...")
    fit_nls <- tryCatch({
      fit_nls <-
        nls2_fit(
          inp = inp,
          probe = probe,
          cores = cores,
          decay = decay,
          delay = delay,
          k = k,
          intyf = intyf
        )
      fit_nls
    },
    error = function(e) {
      writeLines(
        paste(
          "An unknown error has appeared!\n",
          e,
          "An emergency output was returned!\n Please rerun nls2_fit manually!"
        )
      )
      fit_nls2 <- data.frame(matrix(nrow = 0, ncol = 5))
      colnames(fit_nls2) <-
        c("ID", "delay", "half_life", "inty_S0", "intyf")
      fit_nls <- list(probe, fit_nls2)
      return(fit_nls)
    })
    probe <- fit_nls[[1]]
    
    message("running TI_fit...")
    res4 <- tryCatch({
      res4 <-
        TI_fit(
          inp = inp,
          probe = probe,
          cores = cores,
          restr = restr,
          k = TI_k,
          decay = TI_decay,
          ti = TI,
          ti_delay = TI_delay,
          rest_delay = TI_rest_delay
        )
      res4
    },
    error = function(e) {
      writeLines(
        paste(
          "An unknown error has appeared!\n",
          e,
          "An emergency output was returned!\n Please rerun TI_fit manually!"
        )
      )
      res3 <- data.frame(matrix(nrow = 0, ncol = 9))
      colnames(res3) <-
        c(
          "delay",
          "ti_delay",
          "half_life",
          "ti_value",
          "TI_termination_factor",
          "synthesis_rate",
          "TI_background",
          "position",
          "ID"
        )
      res4 <- list(probe, res3)
      return(res4)
    })
    probe <- res4[[1]]
    
    tryCatch({
      if (viz == TRUE) {
        message("running visualization...")
        plot_nls2_function(
          data = fit_nls[[2]],
          inp = inp,
          cores = cores,
          color = color
        )
        plot_singleProbe_function(
          data = res4[[2]],
          inp = inp,
          cores = cores,
          color = color
        )
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
    
    res <- probe
    
    if (details == TRUE) {
      res <- list(probe, fit_nls[[2]], res4[[2]])
      names(res) <- c("probe_df", "fit_obj_STD", "fit_obj_TI")
    }
    res
  }
