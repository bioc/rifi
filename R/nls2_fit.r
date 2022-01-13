#' nls2_fit: estimates decay for each probe or bin.
#' nls2_fit uses nls2 function to fit a probe or bin using intensities of the
#' time series data from different time point.
#' nls2 uses different starting values through expand grid and selects the best
#' fit. Different filters could be applied prior fitting to the model.
#' To apply nls2_fit function, prior filtration could applied.
#' 1. generic_filter_BG: filter probes with intensities below background using
#' threshold. Those probes are flagged with "_FLT_".
#' 2. filtration_below_backg: additional functions exclusive to microarrays
#' could be applied. Its very strict to the background (not recommended in
#' usual case).
#' 3. filtration_above_backg: selects probes with a very high intensity and
#' above the background (recommended for special transcripts). Probes are
#' flagged with "_ABG_".
#' Those transcripts are usually related to a specific function in bacteria.
#' This filter selects all probes with the same ID,
#' the mean is applied, the last time point is selected and compared to the
#' threshold.
#' The function used is:
#' relativity function: normalizes intensities of all time points to 1.
#'
#' the model used estimates the delay, decay, intensity of the first time
#' point (synthesis rate/decay) and the background.
#' The coefficients are gathered in vectors with the corresponding IDs.
#' Absence of the fit or a very bad fit are assigned with NA.
#' In case of probes with very high intensities and above the background,
#' the model used makes abstinence of background coefficient.
#' The output of all coefficients is saved as object list with dataframe an
#' the corresponding coefficients.
#' The fits are plotted using the function_plot_fit.r through rifi_fit.
#'
#' @param inp data frame: the input data frame with correct format.
#' @param probe data frame: the input dataframe containing the relative
#' intensities for each timepoint as well as information about ID, position,
#' strand, probe_TI and flag.
#' @param cores integer: the number of assigned cores for the task.
#' @param decay numeric vector: A sequence of starting values for the decay.
#' Default is seq(.08, 0.11, by=.02)
#' @param delay numeric vector: A sequence of starting values for the delay.
#' Default is seq(0,10, by=.1)
#' @param k numeric vector: A sequence of starting values for the synthesis
#' rate. Default is seq(0.1,1,0.2)
#' @param intyf numeric vector: A sequence of starting values. Default is 0.2.
#'
#' @return A list of 3 data frames :
#' \describe{
#'   \item{probe_df:}{the probe data frame with the columns delay and half_life
#'   added:
#'   \describe{
#'     \item{ID:}{The bin/probe specific ID}
#'     \item{position:}{The bin/probe specific position}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{intensity:}{The relative intensity at time point 0}
#'     \item{probe_TI:}{An internal value to determine which fitting model
#'     is applied}
#'     \item{flag:}{Information on which fitting model is applied}
#'     \item{postion_segment:}{The position based segment}
#'     \item{delay:}{The delay value of the bin/probe}
#'     \item{half_life:}{The half-life of the bin/probe}
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
#' }
#' 
#' @examples
#' data(preprocess_minimal)
#' nls2_fit(
#'   inp = preprocess_minimal$input_df, probe = preprocess_minimal$probe_df,
#'   cores = 2, decay = seq(.08, 0.11, by = .02), delay = seq(0, 10, by = .1),
#'   k = seq(0.1, 1, 0.2), intyf = 0.2
#' )
#' 
#' @export

nls2_fit <-
  function(inp,
           probe,
           cores,
           decay = seq(.08, 0.11, by = .02),
           delay = seq(0, 10, by = .1),
           k = seq(0.1, 1, 0.2),
           intyf = 0.2) {
    num_args <- list(cores)
    names(num_args) <- c("cores")
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
    vec_args <- list(decay, delay, k)
    names(vec_args) <- c("decay", "delay", "k")
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
    
    time <-
      as.numeric(colnames(inp)[seq_len(which(colnames(inp) %in% "ID") - 1)])
    data_i.1 <- inp[!grepl("FLT", inp$filtration), ]
    data_f.1 <- probe[!grepl("_TI_", probe$flag), ]
    ids_ABG <- data_f.1[grepl("ABG", data_f.1$flag), "ID"]
    Data <- data_i.1[data_i.1$ID %in% data_f.1$ID, ]
    probe$delay[probe$ID %in% Data$ID] <- NA
    probe$half_life[probe$ID %in% Data$ID] <- NA
    n_fit <- mclapply(c(seq_along(unique(Data$ID))), function(i) {
      # assign objects and empty vectors for the fit and the coefficients.
      decay_v <- c()
      delay_v <- c()
      ID_v <- c()
      pos_v <- c()
      inty_S0_v <- c()
      intyf_v <- c()
      Data_fit <- data.frame()
      data_c <- data.frame()
      
      # assemble all IDs in a vector
      ID_v <- c(ID_v, unique(Data$ID)[i])
      ID_v <- unique(ID_v)
      pos_v <- c(pos_v,unique(Data$position)[i])
      # I. input preparation:
      # 1. Time and intensity values are assembled in dataframe and used as
      # the first input for the fit.
      # 2. All probes with the same IDs are collected as in input to generate
      # the dataframe.
      tmp <- Data[which(Data$ID %in% unique(Data$ID)[i]), ]
      
      Data_fit_function <- function(data) {
        # intensities values are transformed to generate the dataframe.
        inty <-
          as.vector(t(data[seq_len(which(colnames(inp) %in% "ID") - 1)]))
        # time vector is adapted to the length of number of replicates and
        # number of intensities time point.
        Time <- rep(time, times = length(inty) / length(time))
        Data_fit <- data.frame(as.vector(t(Time)), as.vector(t(inty)))
        # NAs are eliminated from the dataframe.
        Data_fit <- na.omit(Data_fit)
        # intensities values are normalized to 1
        Data_fit[, 2] <- relativity(Data_fit[, 2])
        colnames(Data_fit) <- c("time", "inty")
        return(Data_fit)
      }
      
      # II. fit using nls2 function with different starting values assigned
      # as st1. The fit object is assigned as halfLE2.
      # nls2 fit object is assigned to NA
      halfLE2 <- NA
      
      # probes with flag different from "_" are selected for the model with
      # background coefficient,
      # otherwise the model without background coefficient is applied.
      if (unique(Data$ID)[i] %in% ids_ABG) {
        Data_fit <- Data_fit_function(tmp)
        st1 <- expand.grid(
          decay = decay,
          delay = delay,
          k = k
        )
        cc <- capture.output(type="message",
                             tryCatch({
                               halfLE2 <- nls2(
                                 inty ~ I(time < delay) * k / decay +
                                   (time >= delay) * I(k / decay * (exp(
                                     -decay * (time - delay)
                                   ))),
                                 data = Data_fit,
                                 algorithm = "port",
                                 control = list(warnOnly = TRUE),
                                 start = st1,
                                 lower = list(decay = 0.01, delay = 0.001)
                               )},
                               error = function(e) {}
                             ))
      } else {
        Data_fit <- Data_fit_function(tmp)
        st1 <- expand.grid(
          decay = decay,
          delay = delay,
          k = k,
          intyf = intyf
        )
        cc <- capture.output(type="message",
                             tryCatch({
                               halfLE2 <- nls2(
                                 inty ~ I(time < delay) * k / decay +
                                   (time >= delay) * I(intyf + (k / decay - intyf) * (exp(
                                     -decay * (time - delay)
                                   ))),
                                 data = Data_fit,
                                 algorithm = "port",
                                 control = list(warnOnly = TRUE),
                                 start = st1,
                                 lower = list(decay = 0.01, delay = 0.001)
                               )},
                               error = function(e) {}
                             ))
      }
      
      # III. gathering the coefficients in vectors. NA is assigned to empty
      # fit or a very bad fit where the coefficients can not be extracted.
      tryCatch({
        if (is.null(halfLE2)[1] | is.na(halfLE2)[1]) {
          decay_v <- c(decay_v, NA)
          delay_v <- c(delay_v, NA)
          intyf_v <- c(intyf_v, NA)
          inty_S0_v <- c(inty_S0_v, NA)
        } else {
          decay_v <- c(decay_v, log(2) / coef(halfLE2)[1])
          delay_v <- c(delay_v, coef(halfLE2)[2])
          inty_S0_v <-
            c(inty_S0_v, (coef(halfLE2)[3] / coef(halfLE2)[1]))
          if (length(coef(halfLE2)) == 4) {
            intyf_v <- c(intyf_v, coef(halfLE2)[4])
          } else {
            intyf_v <- c(intyf_v, NA)
          }
        }
      },
      warning = function(war) {
        print(paste("my warning in processing HalfLE2:", i, war))
      },
      error = function(err) {
        print(paste("my error in processing HalfLE2:", i, err))
      }
      )
      data_c <- cbind(ID_v, pos_v, delay_v, decay_v, inty_S0_v, intyf_v)
      data_c <- unique(data_c)
      colnames(data_c) <-
        c("ID", "position", "delay", "half_life", "inty_S0", "intyf")
      return(data_c)
    }, mc.preschedule = FALSE, mc.cores = cores)
    fit_nls2 <- as.data.frame(do.call(rbind, n_fit))
    if (length(n_fit) == 0) {
      fit_nls2 <- data.frame(matrix(nrow = 0, ncol = 6))
      colnames(fit_nls2) <-
        c("ID", "position", "delay", "half_life", "inty_S0", "intyf")
    }
    probe <- probe[with(probe, order(probe$ID)), ]
    fit_nls2 <-
      fit_nls2[with(fit_nls2, order(fit_nls2$ID)), ]
    
    probe$delay[probe$ID %in% Data$ID] <- fit_nls2$delay
    probe$half_life[probe$ID %in% Data$ID] <- fit_nls2$half_life
    
    probe <- probe[with(probe, order(-xtfrm(probe$strand), probe$position)), ]
    
    res <- list(probe, fit_nls2)
    return(res)
  }