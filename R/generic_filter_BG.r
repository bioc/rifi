generic_filter_BG <- function(tmp, threshold=4500, bg_thrsh=250){
  time <- as.numeric(colnames(tmp)[seq_along(which(colnames(tmp) %in%
                                                     "ID")-1)])
  for(j in seq_len(nrow(tmp))){
  tmp_r <- tmp[j, seq_along(which(colnames(tmp) %in% "ID")-1)]
  Exp <- abs(tmp_r[1]-tmp_r[7])
  Max <- which.is.max(tmp_r)
  Min <- which.min(tmp_r)
  res<-"_"
  if((Exp < bg_thrsh &  tmp_r[1] < threshold) | Min == which.min(time) |
     Max == which.is.max(time)){
    res<-"_FLT_"
    }
  }
  res
}
