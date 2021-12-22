#' check_input: reviews the input given by the user.
#' 'check_input' stops the operation if the input data frame has severe faults.
#' Less severe faults lead to the removal of wrong IDs and a warnings describing
#' the problem.
#' The input data frame with all possible replicates need to have the following
#' columns:
#' "0","'time point x'",...,"'time point n'","ID","position","strand"!
#' It is required that the names referring to the time points can be converted
#' to numeric.
#' IDs can be given as character (except ",","|","_"), or numeric, but need to
#' refer to a unique position/strand combination.
#' Strand information need to be given by "+" and "-".
#' The columns "ID","position","strand" must be the last three columns in this
#' order!
#' The relative intensity for all time points should be numeric.
#' The relative intensity for the first time point cannot be 0 or NA.
#'
#' @param inp data frame: the input data frame with correct format.
#' @param thrsh numeric: the minimal allowed intensity for time point "0"...
#' ...advised to be kept at 0! Default is 0.
#'
#' @return the input data frame, with the filtration column added:
#'  \describe{
#'     \item{...}{all timepoints}
#'     \item{ID:}{unique IDs}
#'     \item{position:}{genome positions}
#'     \item{strand:}{strand information}
#'     \item{filtration:}{indicator wether the replicate is filtered or not}
#'     }
#'     
#' @examples
#' data(example_input_minimal)
#' check_input(inp = example_input_minimal, thrsh = 0)
#' @export

check_input <- function(inp, thrsh = 0) {
  num_args <- list(thrsh)
  names(num_args) <- c("trsh")
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
  req_cols_inp <- c("0", "ID", "position", "strand")
  assert(
    all(req_cols_inp %in% colnames(inp)),
    paste0("'", req_cols_inp[which(!req_cols_inp %in% colnames(inp))],
           "' must be a column in 'inp'!")
  )
  time <- names(inp[seq_len(which(colnames(inp) %in% "ID") - 1)])
  # checks if time is given as numeric
  if (any(is.na(as.numeric(time)))) {
    inp <- inp[, -c(which(is.na(as.numeric(time))))]
    stop(
      "A timeslot in the column names is not given as a numeric.\n
      Please review the input dataframe!"
    )
  }
  # check if "ID","position","strand" are there correctly
  if (!all(names(inp[(length(time) + 1):(length(time) + 3)]) ==
           c("ID", "position", "strand"))) {
    stop(
      "A column for ID, position or strand is not given correctly.\n
      Please review the input dataframe!"
    )
  }
  # the columns are ordered in a way, that t0 is in first position
  time <- as.numeric(time)
  names(time) <- seq_along(time)
  ord <- as.numeric(names(sort(time)))
  ord <- c(ord, seq((max(ord) + 1), ncol(inp), 1))
  inp <- inp[, c(ord)]
  # checks if the first time point is 0
  if (names(inp)[1] != "0") {
    stop("No timeslot for timepoint 0 is given.\nPlease review the input
         dataframe!")
  }
  # checks the strand
  if (!(all(unique(inp$strand) %in% c("+", "-")) |
        all(unique(inp$strand) %in% c("-", "+")) |
        all(unique(inp$strand) %in% c("-")) |
        all(unique(inp$strand) %in% c("+")))) {
    symb <-
      paste(unique(inp$strand)[grep("\\+|\\-", unique(inp$strand),
                                    invert = TRUE)])
    inp <-
      inp[-c(which(!(inp$strand == "+" | inp$strand == "-"))), ]
    warning(
      "Your strand information contained the symbol '",
      symb, "'. Rows with this symbol will be removed!\n"
    )
    if (nrow(inp) == 0) {
      stop("All strand information was given in a wrong way!")
    }
  }
  a <- length(unique(paste(inp$position, inp$strand)))
  b <- length(unique(inp$ID))
  if (a < b) {
    warning(
      "There are more unique IDs than positions!\nThe same position will be
      evaluated multiple times!\nPotentially, replicates are not considered
      correctly!"
    )
  }
  if (a > b) {
    stop("There are more positions than unique IDs!\nTwo positions cannot
         share the same ID!")
  }
  if (length(grep("\\||_|\\,", inp$ID)) != 0) {
    bad_IDs <- unique(inp$ID[grep("\\||_|\\,", inp$ID)])
    warning(paste0(
      "IDs cannot contain |, _ or ,.The following IDs have been removed: ",
      paste0(bad_IDs, collapse = ", ")
    ))
  }
  # makes mean for t0
  zero <- split(inp, inp$ID)
  zero_tmp <-
    unlist(lapply(zero, function(x) {
      mean(unlist(x[1]), na.rm = TRUE)
    }))
  # checks if t0 is 0 or NA
  if (any(is.na(zero_tmp)) | any(zero_tmp <= thrsh)) {
    bad_IDs <- names(zero_tmp)[which(is.na(zero_tmp) |
                                       zero_tmp <= thrsh)]
    inp <- inp[-c(which(inp$ID %in% bad_IDs)), ]
    warning(
      paste0(
        "The input contains values below or equal to ",
        thrsh,
        " or NAs at time point 0! The following IDs have been removed: ",
        paste0(bad_IDs, collapse = ", ")
      )
    )
    if (nrow(inp) == 0) {
      stop("The whole input consists of values below or equal to ",
           thrsh,
           " or NAs at time point 0!")
    }
  }
  inp$filtration <- "_"
  res <- list(inp)
  names(res) <- c("input_df")
  res
}
