#' TI_fit: estimates transcription interference and termination factor using
#' nls function for probe or bin flagged as "TI".
#' TI_fit uses nls function to fit the flagged probes or bins with "TI" found
#' using finding_TI.r.
#' It estimates the transcription interference level (referred later to TI) as
#' well as the transcription factor fitting the probes/bins with nls function
#' looping into several starting values.
#' To determine TI and termination factor, TI_fit function is applied to the
#' flagged probes and to the probes localized 200 nucleotides upstream.
#' Before applying TI_fit function, some probes/bins are filtered out if they
#' are below the background using generic_filter_BG.
#' The model loops into a dataframe containing sequences of starting values and
#' the coefficients are extracted from the fit with the lowest
#' residuals. When many residuals are equal to 0, the lowest residual can not
#' be determined and the coefficients extracted could be wrong.
#' Therefore, a second filter was developed. First we loop into all starting
#' values, we collect nls objects in tmp_v vector and the corresponding
#' residuals in tmp_r vector. They are sorted and residuals non equal to 0
#' are collected in a vector. If the first residuals are not equal to 0,
#' 20 % of the best residuals are collected in tmp_r_min vector and the
#' minimum termination factor is selected. In case the first residuals are
#' equal to 0 then values between 0 to 20% of the values collected in tmp_r_min
#' vector are gathered. The minimum termination factor coefficient is
#' determined and saved. The coefficients are gathered in res vector and saved
#' as an object.
#'
#' @param inp data frame: the input data frame with correct format.
#' @param probe data frame: the probe based data frame.
#' @param cores integer: the number of assigned cores for the task.
#' @param restr numeric: a parameter that restricts the freedom of the fit
#' to avoid wrong TI-term_factors, ranges from 0 to 0.2.
#' @param k numeric vector: A sequence of starting values for the synthesis
#' rate. Default is seq(0, 1, by = 0.5).
#' @param decay numeric vector: A sequence of starting values for the decay
#' Default is c(0.05, 0.1, 0.2, 0.5, 0.6).
#' @param ti numeric vector: A sequence of starting values for the delay.
#' Default is  seq(0, 1, by = 0.5).
#' @param ti_delay numeric vector: A sequence of starting values for the
#' delay.
#' Default is seq(0, 2, by = 0.5).
#' @param rest_delay numeric vector: A sequence of starting values. Default
#' is seq(0, 2, by = 0.5).
#'
#' @return A list of 2 data frames :
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
#' @examples 
#' data(preprocess_minimal)
#' TI_fit(inp = preprocess_minimal$input_df, probe = preprocess_minimal$probe_df,
#'  cores=2, restr=0.01)
#'
#' @export

TI_fit <-
  function(inp,
           probe,
           cores = 1,
           restr = 0.2,
           k = seq(0, 1, by = 0.5),
           decay = c(0.05, 0.1, 0.2, 0.5, 0.6),
           ti = seq(0, 1, by = 0.5),
           ti_delay = seq(0, 2, by = 0.5),
           rest_delay = seq(0, 2, by = 0.5)) {
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
    vec_args <- list(k, decay, ti, ti_delay, rest_delay)
    names(vec_args) <-
      c("k", "decay", "ti", "ti_delay", "rest_delay")
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
          (is.numeric(x) & is.vector(x))
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
    # I. inputs
    # vector of different timing point
    time <-
      as.numeric(colnames(inp)[seq_len(which(colnames(inp) %in% "ID") - 1)])

    # II. extract probes IDs
    # extract the IDs of the flagged probes
    corr_IDs <- probe[grepl("TI", probe$flag), "ID"]
    ABG_IDs <-
      probe[grepl("TI", probe$flag) &
              grepl("ABG", probe$flag), "ID"]
    section <- inp[inp$ID %in% corr_IDs, ]

    # III. filtration
    section <- section[!grepl("FLT", section$filtration), ]

    probe$delay[probe$ID %in% corr_IDs] <- NA
    probe$half_life[probe$ID %in% corr_IDs] <- NA
    probe$TI_termination_factor[probe$ID %in% corr_IDs] <- NA

    # IV. generate the dataframe input with the corresponding intensity
    # values
    inty <-
      as.vector(t(section[seq_len(which(colnames(section) %in% "ID") - 1)]))
    # time vector is adapted to the length of number of replicates and number
    # of intensities time point.
    Time <- rep(time, times = length(inty) / length(time))
    id <-
      as.vector(t(rep(
        section$ID,
        each = which(colnames(section) %in% "ID") - 1
      )))
    position <-
      as.vector(t(rep(
        section$position,
        each = which(colnames(section) %in% "ID") - 1
      )))
    input_data <-
      data.frame(as.vector(t(Time)), as.vector(t(inty)), as.vector(t(id)),
                 as.vector(t(position)))
    # NAs are eliminated from the dataframe.
    input_data <- na.omit(input_data)
    colnames(input_data) <- c("time", "inty", "id", "position")

    # V. TI model
    model <-
      inty ~ I(time < ti_delay) * I(k / decay - ti / decay + bg) +
      I(time < ti_delay + rest_delay & time >= ti_delay) *
      I(k / decay - ti / decay * exp(-decay * (time - ti_delay)) + bg) +
      I(time >= ti_delay + rest_delay) *
      I(bg + (k / decay - ti / decay * exp(-decay * rest_delay)) *
          exp(-decay * (time - (ti_delay + rest_delay))))

    model2 <-
      inty ~ I(time < ti_delay) * I(k / decay - ti / decay) +
      I(time < ti_delay + rest_delay & time >= ti_delay) *
      I(k / decay - ti / decay * exp(-decay * (time - ti_delay))) +
      I(time >= ti_delay + rest_delay) *
      I((k / decay - ti / decay * exp(-decay * rest_delay)) *
          exp(-decay * (time - (ti_delay + rest_delay))))

    # VI. list of the starting values generated from a loop
    start_list <-
      "list(ti_delay = pars[i,4], k = rep(pars[i,1],length(gr)),
                    decay=pars[i,2], ti = pars[i,3],
                    rest_delay = pars[i,5], bg = rep(0,length(gr)))"
    start_list2 <-
      "list(ti_delay = pars[i,4], k = rep(pars[i,1],length(gr)),
                    decay=pars[i,2], ti = pars[i,3],
                    rest_delay = pars[i,5])"
    lo <- "rep(0, ncol(pars) + 1)"
    lo2 <- "rep(0, ncol(pars))"
    k <- k
    decay <- decay
    ti <- ti
    ti_delay <- ti_delay
    rest_delay <- rest_delay
    pars <- expand.grid(k, decay, ti, ti_delay, rest_delay)

    temp_res <- data.frame()
    delay_Ve <- c()
    k_Ve <- c()
    decay_Ve <- c()
    decay2_Ve <- c()
    ti_Ve <- c()
    ti2_Ve <- c()
    bgs_Ve <- c()
    positions_Ve <- c()
    ids_Ve <- c()
    ti_delay_Ve <- c()
    term_prob_Ve <- c()

    gr1 <- unique(input_data$id)

    r <- mclapply(gr1, function(j) {
      gr <- j
      out <- NA
      temp_data <- input_data[which(input_data$id == gr), ]
      # normalize intensities values to 1
      temp_data$inty <- temp_data$inty / temp_data$inty[1]

      tmp_r <- c()
      tmp_v <- c()
      # vector gathering values above 0
      value <- c()

      # VII. fitting the probes
      for (i in seq_len(nrow(pars))) {
        tmp <- NA
        if (gr %in% ABG_IDs) {
          tmp <- tryCatch({
            nls(
              model2,
              data = temp_data,
              algorithm = "port",
              control = list(warnOnly = TRUE),
              start = eval(parse(text = start_list2)),
              lower = eval(parse(text = lo2))
            )
          },
          error = function(e) {},
		  warning = function(e){}
		  )
        } else {
          tmp <- tryCatch({
            nls(
              model,
              data = temp_data,
              algorithm = "port",
              control = list(warnOnly = TRUE),
              start = eval(parse(text = start_list)),
              lower = eval(parse(text = lo))
            )
          },
          error = function(e) {},
		  warning = function(e){}
		  )
        }
        if (any(!is.na(tmp))) {
          if (!is(tmp, "simpleError")) {
            if (any(is.na(out))) {
              out <- list(tmp)
            } else {
              # gather all nls objects and select the 20% of those with best
              # coefficient and lowest ti
              # since we did find a high ti with the lowest coefficients
              # which
              # is wrong
              tmp_v <- c(tmp_v, list(tmp))
              tmp_r <- c(tmp_r, sum(residuals(tmp) ** 2))
              names(tmp_r) <- seq_along(tmp_r)
              tmp_r_s <- sort(tmp_r, decreasing = FALSE)
              for (l in seq_along(tmp_r_s)) {
                if (tmp_r_s[l] != 0) {
                  value <- c(value, tmp_r_s[l])
                }
              }
              value <- unique(value)
              if (tmp_r_s[1] != 0) {
                tmp_r_min <-
                  tmp_r_s[between(tmp_r_s, tmp_r_s[1], (tmp_r_s[1] * restr +
                                                          tmp_r_s[1]))]
              } else {
                tmp_r_min <-
                  tmp_r_s[between(tmp_r_s, tmp_r_s[1], (value[1] * restr +
                                                          value[1]))]
              }
              if (is_empty(tmp_r_min)) {
                next ()
              }
              out <-
                tmp_v[as.numeric(as.character(names(tmp_r_min)))]
            }
          }
        }
      }

      # VIII. Collect the best residuals from nls object from out and extract
      # the coefficients
      ti_delay_v <- c()
      k_v <- c()
      decay_v <- c()
      ti_v <- c()
      bgs_v <- c()
      rest_delay_v <- c()

      # loop into out and extract the coefficients
      for (k in seq_along(out)) {
        if (any(!is.na(out[[k]]))) {
          fitted_pars <- coef(out[[k]])
          ti_delay_v <- c(ti_delay_v, fitted_pars[1])
          k_v <- c(k_v, fitted_pars[2])
          decay_v <- c(decay_v, fitted_pars[3])
          ti_v <- c(ti_v, fitted_pars[4])
          bgs_v <- c(bgs_v, fitted_pars[6])
          rest_delay_v <- c(rest_delay_v, fitted_pars[5])
        }
      }

      # determine the lowest termination factor
      minV <- which.min(ti_v)
      ti_delay <- ti_delay_v[minV]
      k <- k_v[minV]
      decay <- decay_v[minV]
      ti <- ti_v[minV]
      bgs <- bgs_v[minV]
      rest_delay <- rest_delay_v[minV]
      delay <- ti_delay + rest_delay
      positions <- unique(temp_data$position)
      ids <- unique(as.numeric(as.character(temp_data$id)))
      term_prob <- ti / k

      if (all(is.na(out))) {
        ti_delay <- NA
        k <- NA
        decay <- NA
        ti <- NA
        bgs <- NA
        rest_delay <- NA
        delay <- NA
        positions <- unique(temp_data$position)
        ids <- unique(as.numeric(as.character(temp_data$id)))
        term_prob <- NA
      }

      delay_Ve <- c(delay_Ve, delay)
      ti_delay_Ve <- c(ti_delay_Ve, ti_delay)
      decay_Ve <- c(decay_Ve, log(2) / decay)
      ti_Ve <- c(ti_Ve, ti)
      term_prob_Ve <- c(term_prob_Ve, term_prob)
      k_Ve <- c(k_Ve, k)
      bgs_Ve <- c(bgs_Ve, bgs)
      positions_Ve <- c(positions_Ve, positions)
      ids_Ve <- c(ids_Ve, ids)

      temp_res <- cbind(
        delay_Ve,
        ti_delay_Ve,
        decay_Ve,
        ti_Ve,
        term_prob_Ve,
        k_Ve,
        bgs_Ve,
        positions_Ve,
        ids_Ve
      )
      colnames(temp_res) <-
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
      return(temp_res)
    }, mc.preschedule = FALSE, mc.cores = cores)
    res3 <- as.data.frame(do.call(rbind, r))

    res3[apply(res3, 2, is.infinite)] <- NA

    if (length(r) == 0) {
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
    }

    probe <- probe[with(probe, order(probe$ID)), ]
    res3 <-
      res3[with(res3, order(res3$ID)), ]

    probe$delay[probe$ID %in% corr_IDs] <- NA
    probe$delay[probe$ID %in% corr_IDs] <- res3$delay
    probe$half_life[probe$ID %in% corr_IDs] <- NA
    probe$half_life[probe$ID %in% corr_IDs] <- res3$half_life
    probe$TI_termination_factor <- NA
    probe$TI_termination_factor[probe$ID %in% corr_IDs] <-
      res3$TI_termination_factor

    probe <-
      probe[with(probe, order(-xtfrm(probe$strand), probe$position)), ]

    res4 <- list(probe, res3)
    return(res4)
  }
