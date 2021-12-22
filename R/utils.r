char_TI <- function(y, hig = 0.075) {
  y <- y / y[1] # normalizes the input
  tmp <-
    max(y, na.rm = TRUE) - y[1] # subtracts t0 point from the highest point
  res <- -1 # not TI exept..
  if (tmp > hig) {
    # ..the difference is higher than 7.5%
    res <- 1 # if yes, the output is 1
  }
  res
}

relativity <- function(x) {
  (x / max(x))
}

finding_above_backg <- function(tmp, bg_lastTime, time) {
  tmp[, "fl_above_backg"] <- "_ABG_"
  for (j in seq_len(nrow(tmp))) {
    tmp_ID <- tmp[j, which(colnames(tmp) %in% "ID") - 1]
    if (!is.na(tmp_ID)) {
      if (tmp_ID < bg_lastTime) {
        tmp[j, "fl_above_backg"] <- "_"
      }
    }
  }
  return(tmp$fl_above_backg)
}

assert <- function(expr, error) {
  if (!expr) {
    stop(error, call. = FALSE)
  }
}
