#' segment_pos: divides all IDs by position into position_segments.
#' segment_pos adds the column "position_segment" to the probe based dataframe.
#' To reduce run time, the data is divided by regions of no expression larger
#' than "dist" nucleotides.
#' The data frame needs to contain at least "strand" and "position".
#' 
#' @param probe data frame: the probe based data frame.
#' @param dista integer: the amount of nucleotides defining the gap. Default
#' is 300.
#' 
#' @return the probe based data frame with the position_segment column added:
#'  \describe{
#'     \item{ID:}{The bin/probe specific ID}
#'     \item{position:}{The bin/probe specific position}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{intensity:}{The relative intensity at time point 0}
#'     \item{probe_TI:}{An internal value to determine which fitting model is
#'     applied}
#'     \item{flag:}{Information on which fitting model is applied}
#'     \item{postion_segment:}{The position based segment}
#'       }
#'       
#' @examples
#' data(preprocess_minimal)
#' segment_pos(probe = preprocess_minimal$probe_df, dista = 300)
#' 
#' @export

segment_pos <- function(probe, dista = 300) {
  num_args <- list(dista)
  names(num_args) <- c("dista")
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
  req_cols_probe <- c("strand", "position")
  assert(
    all(req_cols_probe %in% colnames(probe)),
    paste0("'", req_cols_probe[which(!req_cols_probe %in% colnames(probe))],
           "' must be a column in 'probe'!")
  )

  count <- 1
  tmp <- probe
  tmp$position_segment <- NA
  strand_info <- c("+", "-", "NA")
  for (k in seq_along(strand_info)) {
    x2 <- probe[which(probe$strand == strand_info[k]), ]
    x2 <- x2[order(x2$position), ]
    pos <- unique(x2$position)
    if (length(pos) != 0) {
      names(pos) <- seq_along(pos)
      pos_d <- diff(pos)
      pos_dv <- c()
      for (i in seq_along(pos_d)) {
        if (pos_d[i] >= dista) {
          pos_dv <- c(pos_dv, pos_d[i])
        }
      }
      pos_dvn <- names(pos_dv)
      pos_dvn <- as.numeric(as.character(pos_dvn))
      splitAt <- function(x, pos) {
        (split(x, cumsum(seq_along(x) %in% pos)))
      }
      fragmentsToAnalyze <- splitAt(pos, pos = pos_dvn)
      fragmentsToAnalyze <- unname(fragmentsToAnalyze)

      for (i in seq_along(fragmentsToAnalyze)) {
        tmp$position_segment[which(tmp$position %in% fragmentsToAnalyze[[i]] &
          tmp$strand == strand_info[k])] <- paste0("S_", count)
        count <- count + 1
      }
    }
  }
  return(tmp)
}
