synthesis_r_Function <- function(x, data) {
  for (k in seq_len(nrow(data))) {
    if (is.na(data[k, x])) {
      next ()
    }
    if (data[k, x] >= 1) {
      data[k, "synthesis_ratio_event"] <- "iTSS_II"
    } else{
      data[k, "synthesis_ratio_event"] <- "Termination"
    }
  }
  data
}
