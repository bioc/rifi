synthesis_r_Function <- function(x, data) {
  for (k in seq_len(nrow(data))) {
    if (is.na(as.data.frame(rowRanges(data))[k, x])) {
      next ()
    }
    if (as.data.frame(rowRanges(data))[k, x] >= 0) {
      rowRanges(data)$synthesis_ratio_event[k] <- "iTSS_II"
    } else {
      rowRanges(data)$synthesis_ratio_event[k] <- "Termination"
    }
  }
  data
}
