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

finding_above_backg <- function(df, bg) {
  res <- rep("_ABG_",nrow(df))
  for (j in seq_len(nrow(df))) {
    tmp_last <- df[j, ncol(df)]
    if (!is.na(tmp_last)) {
      if (tmp_last < bg) {
        res[j] <- "_"
      }
    }
  }
  res
}

assert <- function(expr, error) {
  if (!expr) {
    stop(error, call. = FALSE)
  }
}

encode_FLT <- function(obj, rows, rep){
  bi <- intToBits(0)
  bi[rep] <- as.raw(1)
  int<-packBits(bi,"integer")
  rowRanges(obj)$FLT[rows] <- bitwOr(rowRanges(obj)$FLT[rows], int)
  obj
}


decode_FLT <- function(obj){
  #get int vector
  int_vec <- rowRanges(obj)$FLT
  #get logical from bits
  logi <- lapply(int_vec, function(x){as.logical(intToBits(x))})
  #representation fo replicates to filter
  repli <- lapply(logi, function(x){metadata(obj)$replicates[x]})
  #logical of correct length
  logi2 <- lapply(repli, function(x){colData(obj)$replicate %in% x})
  #logical matrix
  log_mat <- do.call(rbind,logi2)
  log_mat
}

inp_order <- function(inp){
  #the columns are ordered in a way, that t0 is in first position
  #the rows are order by strand and position increasing(plus) decreasing(minus)
  assert(all(c("timepoint", "replicate") %in% names(colData(inp))),
       "timepoint and replicate must be columns in the colData")
  assert("position" %in% names(mcols(rowRanges(inp))),
         "position must be columns in the rowRanges")
  ord_col <- order(colData(inp)$replicate,colData(inp)$timepoint)
  ord_row <- order(strand(inp))
  inp <- inp[ord_row, ord_col]
  ord_plus <- order(rowRanges(inp[strand(inp) == "+",])$position)
  ord_minus <- order(rowRanges(inp[strand(inp) == "-",])$position,
                     decreasing = TRUE)
  inp[strand(inp) == "+",] <- inp[strand(inp) == "+",][ord_plus,]
  inp[strand(inp) == "-",] <- inp[strand(inp) == "-",][ord_minus,]
  inp
}

inp_df <- function(inp, ...){
  cols <- c(..., "FLT")
  assert(all(cols %in% names(mcols(rowRanges(inp)))),
         paste0("one of the necessary columns in rowRanges is missing: ",
                paste0(cols,collapse = ", ")))
  #get the integer representing all replicates
  bi<-intToBits(0)
  bi[metadata(inp)$replicate] <- as.raw(1)
  int<-packBits(bi,"integer")
  #get the data frame
  tmp_df <- rowRanges(inp)[,cols]
  tmp_df <- as.data.frame(mcols(tmp_df))
  tmp_df$strand <- decode(strand(inp))
  tmp_df$position <- rowRanges(inp)$position
  #filter the data frame
  tmp_df <- tmp_df[tmp_df$FLT != int,]
  tmp_df
}

tmp_df_rev <- function(tmp_df, stra){
  tmp_df[tmp_df$strand == stra, "position"] <-
    ((tmp_df[tmp_df$strand == stra,
             "position"]) - (tmp_df[tmp_df$strand == stra, ]
                             [nrow(tmp_df[tmp_df$strand == stra, ]),
                               "position"])) * -1
  tmp_df
}
