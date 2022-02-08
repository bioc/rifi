fragment_function <- function(seg) {
  segs <- c()
  unique_Input <- unique(seg)
  for (j in seq_along(unique_Input)) {
    ind_input <- which(unique_Input[j] == seg)
    if (length(ind_input) == 1) {
      next ()
    } else{
      segs <- c(segs, unique_Input[j])
    }
  }
  return(segs)
}

t_test_function <-
  function(data, seg, frag, param, o, tu, threshold) {
    if (length(seg) > 1) {
      for (j in seq_along(seg) - 1) {
        frg.1.h <- tu[which(tu[, frag] == seg[j]), ]
        frg.2.h <- tu[which(tu[, frag] == seg[j + 1]), ]
        #if the dataframe has less than 3 probes/bins, will not be subjected to
        # t-test
        if (nrow(frg.1.h) < 3 | nrow(frg.2.h) < 3) {
          next ()
        } else if (abs(last(frg.1.h$position) - frg.2.h$position[1]) <=
                   threshold) {
          frg.1.h <- tu[which(tu[, frag] == seg[j]), param]
          frg.2.h <- tu[which(tu[, frag] == seg[j + 1]), param]
          t_h <-
            t.test(frg.1.h,
                   frg.2.h,
                   alternative = "two.sided",
                   var.equal = FALSE)
          #write fragments name separated by ":" subjected to t-test
          frag_hl <- paste0(seg[j], ":", seg[j + 1])
          #extract the p_value from t-test
          t_h <- t_h[[3]]
          #Fold change of the neighboring fragments, the mean of the
          #fragment2/the mean of the fragment1
          quot.hl <- log2(mean(frg.2.h)/mean(frg.1.h))
          #assign the fragments name on the corresponding column
          data[which(data[, frag] == seg[j]), paste0("FC_fragment_", o)] <-
            rep(frag_hl, times = length(which(data[, frag] == seg[j])))
          #assign the ratio of FC on the corresponding column
          data[which(data[, frag] == seg[j]), paste0("FC_", o)] <-
            rep(quot.hl, times = length(which(data[, frag] == seg[j])))
          #assign the p_value on the corresponding column
          data[which(data[, frag] == seg[j]), paste0("p_value_", o)] <-
            rep(t_h, times = length(which(data[, frag] == seg[j])))
        } else{
          next ()
        }
      }
    }
    return(data)
  }

manova_function <- function(x, y, z, data) {
  res.man <- manova(cbind(x, y) ~ z, data = data)
  p_value_Manova <- summary(res.man)[[4]]
  p_value_Manova <- p_value_Manova[1, 6]
  return(p_value_Manova)
}
