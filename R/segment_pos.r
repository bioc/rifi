#' =========================================================================
#' segment_pos
#' -------------------------------------------------------------------------
#'segment_pos divides all IDs by position into position_segments
#'
#' segment_pos adds the column "position_segment" to the rowRanges.
#' To reduce run time, the data is divided by regions of no expression larger
#' than "dist" nucleotides.
#' 
#' @param inp SummarizedExperiment: the input.
#' @param dista integer: the amount of nucleotides defining the gap. Default
#' is 300.
#' 
#' @return The SummarizedExperiment object: with position_segment added to the
#' rowRanges.
#'       
#' @examples
#' data(preprocess_minimal)
#' segment_pos(inp = preprocess_minimal, dista = 300)
#' 
#' @export

segment_pos <- function(inp, dista = 300) {
  assert(is.numeric(dista) & length(dista) == 1,
         "dista must be numeric of length one")
  #order the input
  inp <- inp_order(inp)
  #make the tmp_df
  tmp_df <- inp_df(inp, "ID")
  #calculate difference
  tmp_diff <- abs(diff(tmp_df$position))
  diff_logi <- c(FALSE, tmp_diff > dista)
  #check strand difference
  strand_logi <- c(FALSE,
                   tmp_df$strand[-1] != tmp_df$strand[-length(tmp_df$strand)])
  #each TRUE increases the number
  number <- cumsum(strand_logi | diff_logi)+1
  position_segment <- paste0("S_", number)
  #fuse with input
  rowRanges(inp)$position_segment <- as.character(NA)
  rowRanges(inp)$position_segment[rowRanges(inp)$ID %in% tmp_df$ID] <-
    position_segment
  inp
}
