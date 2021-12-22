#' rifi_preprocess: conveniently wraps all pre-processing steps.
#' Wraps the functions:
#' 1. check_input
#' 2. make_df
#' 3. function_seg
#' 4. finding_PDD
#' 5. finding_TI
#' Allows for the optional integration of filter functions.
#' Filter functions mark replicates with "FLT". Those are then not considered
#' in the fit!
#' Three different filter functions can be applied. The first (FUN_filter_BG)
#' is a general filter usually to exclude probes with low expression or "bad"
#' patterns.
#' The second (FUN_filter_STD) only filters IDs that will be fitted with the
#' standard fit. The third (FUN_filter_TI) only filters IDs that will be fitted
#' with the TI fit
#'
#' @param inp data frame: the input data frame with correct format.
#' @param cores integer: the number of assigned cores for the task.
#' @param FUN_filter_BG function: A function of x, returning a character
#' string containing "FLT" (e.g ("FLT_BG_3_7")). x is the numeric vector of the
#' intensity from all time points.
#' @param FUN_filter_STD function: A function of x, returning a character string
#' containing "FLT" (e.g ("FLT_STD_alpha")). x is the numeric vector of the
#' intensity from all time points.
#' @param FUN_filter_TI function:A function of x, returning a character string
#' containing "FLT" (e.g ("FLT_TI_Z")). x is the numeric vector of the intensity
#' from all time points.
#' @param bg numeric: threshold over which the last time point has to be to be
#' fitted with the above background mode.
#' @param rm_FLT logical: remove IDs where all replicates are marked as filtered
#' by the background check. Default is FALSE.
#' @param thrsh_check numeric: the minimal allowed intensity for time point "0".
#' Advised to be kept at 0! Default is 0.
#' @param dista integer: the amount of nucleotides defining the gap. Default is
#' 300.
#' @param run_PDD logical: running the PDD flag function
#' @param pen_PDD numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer fragments. Advised to be kept at 2. Default
#' is 2.
#' @param pen_out_PDD numeric: an internal parameter for the dynamic
#' programming. Higher values result in fewer possible outliers. Advised to
#' be kept at 1. Default is 1.
#' @param thrsh_PDD numeric: an internal parameter that allows fragments with
#' slopes steeper than the threshold to be flagged with "_PDD_". Higher values
#' result in fewer candidates . Advised to be kept at 0.001. Default is 0.001.
#' @param pen_TI numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer fragments. Advised to be kept at 10. Default
#' is 10.
#' @param thrsh_TI numeric: an internal parameter that allows fragments with a
#' certain amount of IDs with higher relative intensities at time points later
#' than "0" to be flagged as "_TI_". Higher values result in fewer candidates.
#' -0.5 is 25 %, 0 is 50%, 0.5 is 75%. Advised to be kept at 0.5. Default
#' is 0.5.
#' @param add integer: range of nucleotides before a potential TI event where
#' in IDs are fitted with the TI fit.
#' 
#' @return A list of 2 data frames:
#' \describe{
#'   \item{probe_df:}{the probe dataframe:
#'   \describe{
#'     \item{ID:}{The bin/probe specific ID}
#'     \item{position:}{The bin/probe specific position}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{intensity:}{The relative intensity at time point 0}
#'     \item{probe_TI:}{An internal value to determine which fitting model is
#'     applied}
#'     \item{flag:}{Information on which fitting model is applied}
#'     \item{postion_segment:}{The position based segment}
#'       }
#'     }
#'   \item{fit_obj_TI:}{the fit object for the TI fit:
#'     \describe{
#'     \item{...}{all timepoints}
#'     \item{ID:}{unique IDs}
#'     \item{position:}{genome positions}
#'     \item{strand:}{strand information}
#'     \item{filtration:}{indicator wether the replicate is filtered or not}
#'     }
#'   }
#' }
#'
#' @seealso `check_input`
#' @seealso `make_df`
#' @seealso `segment_pos`
#' @seealso `finding_PDD`
#' @seealso `finding_TI`
#'
#' @examples
#' data(example_input_minimal)
#' rifi_preprocess(
#'   inp = example_input_minimal, cores = 2, bg = 0, rm_FLT = FALSE,
#'   thrsh_check = 0, dista = 300, run_PDD = FALSE
#'   )
#'   
#' @export

rifi_preprocess <-
  function(inp,
           cores,
           FUN_filter_BG = function(x) {
             "_"
           },
           FUN_filter_STD = function(x) {
             "_"
           },
           FUN_filter_TI = function(x) {
             "_"
           },
           bg = 0,
           rm_FLT = FALSE,
           thrsh_check = 0,
           dista = 300,
           run_PDD = FALSE,
           pen_PDD = 2,
           pen_out_PDD = 1,
           thrsh_PDD = 0.001,
           pen_TI = 10,
           thrsh_TI = 0.5,
           add = 1000) {
    num_args <-
      list(
        cores,
        bg,
        thrsh_check,
        dista,
        pen_PDD,
        pen_out_PDD,
        thrsh_PDD,
        pen_TI,
        thrsh_TI,
        add
      )
    names(num_args) <-
      c(
        "cores",
        "bg",
        "thrsh_check",
        "dista",
        "pen_PDD",
        "pen_out_PDD",
        "thrsh_PDD",
        "pen_TI",
        "thrsh_TI",
        "add"
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
    FUN_args <- list(FUN_filter_BG, FUN_filter_STD, FUN_filter_TI)
    names(FUN_args) <-
      c("FUN_filter_BG", "FUN_filter_STD", "FUN_filter_TI")
    assert(all(unlist(lapply(
      FUN_args, is.function
    ))), paste0("'", names(which(
      unlist(lapply(FUN_args, is.numeric)) == FALSE
    )[1]), "' must be a function"))
    assert(is.logical(rm_FLT), "'rm_FLT' must be a logical")
    assert(cores > 0, "'cores' must be a positive integer")
    req_cols_inp <- c("0", "ID", "position", "strand")
    assert(
      all(req_cols_inp %in% colnames(inp)),
      paste0("'", req_cols_inp[which(!req_cols_inp %in% colnames(inp))],
             "' must be a column in 'inp'!")
    )
    message("running check_input...")
    tmp <- check_input(inp = inp, thrsh = thrsh_check)
    inp <- tmp[[1]]
    inp_save <- inp
    time <-
      as.numeric(colnames(inp)[seq_len(which(colnames(inp) %in% "ID") - 1)])
    inp <- tryCatch({
        for (i in seq_len(nrow(inp))) {
          inp[i, "filtration"] <-
            FUN_filter_BG(x = inp[i,
                                  seq_len(which(colnames(inp) %in% "ID") - 1)])
        }
        inp
      },
      error = function(e) {
        writeLines(
          paste(
            "An unknown error has appeared!\n",
            e,
            "An emergency output was returned!\n The given filtration
            malfunctioned!"
          )
        )
        inp <- inp_save
        return(inp)
      }
    )

    message("running make_df...")
    tryCatch({
        probe <- make_df(
          inp = inp,
          cores = cores,
          bg = bg,
          rm_FLT = rm_FLT
        )
      },
      error = function(e) {
        stop(e, call. = FALSE)
      }
    )

    message("running segment_pos...")
    tryCatch({
        probe <- segment_pos(probe = probe, dista = dista)
      },
      error = function(e) {
        stop(e, call. = FALSE)
      }
    )
    if (run_PDD == TRUE) {
      message("running finding_PDD...")
      tryCatch({
          probe <-
            finding_PDD(
              probe = probe,
              pen = pen_PDD,
              pen_out = pen_out_PDD,
              thrsh = thrsh_PDD,
              cores = cores
            )
        },
        error = function(e) {
          writeLines(
            paste(
              "An unknown error has appeared!\n",
              e,
              "An emergency output was returned!\n Please rerun finding_PDD
              manually!"
            )
          )
        }
      )
    }

    message("running finding_TI...")
    tryCatch({
        probe <-
          finding_TI(probe,
            pen = pen_TI,
            thrsh = thrsh_TI,
            cores,
            add = add
          )
      },
      error = function(e) {
        writeLines(
          paste(
            "An unknown error has appeared!\n",
            e,
            "An emergency output was returned!\n Please rerun finding_TI
            manually!"
          )
        )
      }
    )
    inp_save <- inp
    inp <- tryCatch({
        ID_STD <- probe$ID[!grepl("TI", probe$flag)]
        ID_TI <- probe$ID[grepl("TI", probe$flag)]
        inp_STD <- inp[which(inp$ID %in% ID_STD), ]
        inp_TI <- inp[which(inp$ID %in% ID_TI), ]

        if (nrow(inp_STD) > 0) {
          for (i in seq_len(nrow(inp_STD))) {
            tmp <- FUN_filter_STD(
              x = inp_STD[i, seq_len(which(colnames(inp_STD) %in% "ID") - 1)])
            inp_STD[i, "filtration"] <-
              paste0(inp_STD[i, "filtration"], tmp)
          }
        }
        if (nrow(inp_TI) > 0) {
          for (i in seq_len(nrow(inp_TI))) {
            inp_TI[i, "filtration"] <-
              paste0(inp_TI[i, "filtration"], FUN_filter_TI(
                x = inp_TI[i, seq_len(which(colnames(inp_TI) %in% "ID") - 1)]))
          }
        }
        inp <- rbind(inp_STD, inp_TI)
        inp
      },
      error = function(e) {
        writeLines(
          paste(
            "An unknown error has appeared!\n",
            e,
            "An emergency output was returned!\n The given filtration
            malfunctioned!"
          )
        )
        inp <- inp_save
        return(inp)
      }
    )
    res <- list(probe, inp)
    names(res) <- c("probe_df", "input_df")
    res
  }
