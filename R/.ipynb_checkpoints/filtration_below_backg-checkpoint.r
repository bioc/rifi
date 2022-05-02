filtration_below_backg <-
  function(tmp,
           firstTimeP = 1200,
           bg_thrsh,
           time,
           verbose = FALSE) {
    tmp[, "filtration"] <- "_PASS_"
    threshold <- mean(tmp[, 1])
    #extracted from a previous run of filtration and mean of time 0 was applied
    for (j in seq_len(nrow(tmp))) {
      print(j)
      tmp_r <- tmp[j, seq_along(which(colnames(tmp) %in% "ID") - 1)]
      #set a difference between first and last intensity from time serie
      Exp <- abs(tmp_r[1] - tmp_r[7])
      tmp_i <- tmp[j, 1]
      Max <- which.is.max(tmp_r)
      Min <- which.min(tmp_r)
      #@assign an empty vector to assemble row numbers, subtracted from tmp
      try(if (Exp < bg_thrsh & tmp_i < threshold) {
        tmp[j, "filtration"] <- "_1_"
      } else if (Min == which.min(time)) {
        tmp[j, "filtration"] <- "_2_"
      } else if (Max == which.is.max(time)) {
        tmp[j, "filtration"] <- "_3_"
      } else if (Max == which.is.max(time) - 1 & tmp_i < threshold) {
        tmp[j, "filtration"] <- "_4_"
      } else if (Max == which.is.max(time) - 2 & tmp_i < threshold) {
        tmp[j, "filtration"] <- "_5_"
      })
      if (verbose == TRUE) {
        print(paste("Exp", j, "is less than the threshold"))
      }
    }
    return(tmp)
  }
