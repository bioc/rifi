# =========================================================================
# Check_input           reviews the input given by the user
# -------------------------------------------------------------------------
#'
#'
#' 'check_input' stops the operation if the input data frame has severe faults.
#' Less severe faults lead to the removal of wrong IDs and a warnings describing
#' the problem.

#' The Summarized Experiment colData must have the columns "timepoint" with
#' the timepoints convertible to numeric and containing the timepoint 0.

#' If replicates are used the column in colData must be called "replicate".
#' The replicate must be convertible to numeric.

#' In the RowRanges, optionally, IDs can be given as character
#' (except ",","|","_"),but need to refer to a unique position/strand
#' combination.

#' Strand information needs to be given.

#' The relative intensity in the assay must be numeric.

#' The relative intensity for the first time point cannot be 0 or NA.
#'
#' @param inp SummarizedExperiment: the input data frame with correct format.
#' @param thrsh numeric: the minimal allowed intensity for time point "0".
#'
#' @return the SummarizedExperiment object: checked, and with position, ID and
#' filtration added to the rowRanges.
#'     
#' @examples
#' data(example_input_minimal)
#' check_input(inp = example_input_minimal, thrsh = 0)
#' @export

check_input <- function(inp, thrsh = 0) {
  assert(is.numeric(thrsh) & length(thrsh) == 1,
  "thrsh must be numeric of length one")
  #checks if time is given as numeric
  assert(is.numeric(colData(inp)$timepoint),
         "The timepoints must be given as numeric.\n
          Please review the input dataframe!")
  assert(0 %in% colData(inp)$timepoint,
         "A timepoint 0 must be given.\n
          Please review the input dataframe!")
  #check replicates
  if(!"replicate" %in% names(colData(inp))){
    colData(inp)$replicate<-1
  }
  repli <- unique(colData(inp)$replicate)
  metadata(inp)$replicates <- repli
  message("Number of replicates: ",paste(max(repli),collapse = " "))
  #create position column
  if(!"position" %in% names(mcols(rowRanges(inp)))){
    rowRanges(inp)$position <- end(resize(inp,width = 1,fix = "end"))
  }
  #order the input
  inp <- inp_order(inp)
  #get the time
  time <- colData(inp)$timepoint
  metadata(inp)$timepoints <- time
  message("Timepoints: ",paste(unique(time),collapse = " "))
  inp <- inp_order(inp)
  if(!"ID" %in% names(mcols(rowRanges(inp)))){
    rowRanges(inp)$ID<-as.character(seq(nrow(inp)))
    message("No IDs were given in the input. Default IDs were assigned.")
  }
  if(!"FLT" %in% names(mcols(rowRanges(inp)))){
    rowRanges(inp)$FLT <- 0
  }
  if("*" %in% decode(strand(inp))){
    rows <- which(decode(strand(inp)) %in% "*")
    rep <- metadata(inp)$replicate
    inp <- encode_FLT(obj = inp, rows = rows, rep = rep)
    warning("Probes without strand information cannot be considered.")
  }
  # checks if t0 is 0 or NA
  zero_tmp <- rowMeans(assay(inp[,inp$timepoint == "0"]), na.rm = TRUE)
  if (any(is.na(zero_tmp) | any(zero_tmp <= thrsh))) {
    rows <- which(is.na(zero_tmp) | zero_tmp <= thrsh)
    rep <- metadata(inp)$replicate
    inp <- encode_FLT(obj = inp, rows = rows, rep = rep)
    warning("Probes below ",thrsh," or NA at timepoint 0 connot be considered.")
  }
  comb<-cbind(decode(strand(inp)),rowRanges(inp)$position)
  if(any(duplicated(comb))){
    inp <- inp[which(!duplicated(comb)),]
    warning("Probes without unique position and strand cannot be considered.")
  }
  inp
}
