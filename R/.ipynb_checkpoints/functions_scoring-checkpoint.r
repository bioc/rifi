# score_fun_linear scores the values of y on how close they are to a linear fit.
# The script also returns the respective slope of the respective fit.
# The output is "IDs_slope_intercept_outliers" as name and score as value.
# y is a vector of values, x is a vector of respective position,
# z is a vector of identifiers (IDs in our case) (is the position by default)
# stran should be TRUE if strand is not (all) NA. n_out is the number
# of allowed outliers.
score_fun_linear <-
  function(y,
           x,
           z = x,
           pen,
           stran,
           n_out = min(10, max(1, 0.2 * length(x)))) {
    if (length(unique(x)) > 1) {
      # a linear regression can not be used on just one position
      n_out <-
        min(n_out, length(x) - 3) # the number of allowed outliers is decided
      mo <- lm(y ~ x) # the linear fit is performed
      mo_save <- mo
      tmp <- abs(residuals(mo)) # the residuals are cached
      tmp_velo <- coef(mo)[[2]] # velocity...
      tmp_inter <- coef(mo)[[1]] # .. and intercept are cached
      if (coef(mo_save)[[2]] < 0 &
        stran == TRUE) {
        # if the stranded option is active and the slope is negative, the
        # residuals are calculated as if the slope was 0.
        tmp <- abs(y - (mean(y))) # the score is overwritten
        tmp_velo <- 0 # velocity...
        tmp_inter <- mean(y) # .. and intercept are overwritten
      }
      if (coef(mo_save)[[2]] < (- (1 / 60)) &
        stran == FALSE) {
        # if the stranded option is inactive and ...
        # ...the slope is smaller - 1/60 the residuals are calculated as if
        # the slope was -1/60
        mo <- lm(y - I((-1 / 60)) * x ~ 1)
        tmp <- abs(residuals(mo))
        tmp_velo <- (-1 / 60) # velocity...
        tmp_inter <- coef(mo)[[1]] # .. and intercept are overwritten
      }
      if (coef(mo_save)[[2]] > (1 / 60)) {
        # if the slope is bigger 1/60, the residuals are calculated as
        # if the slope was 1/60
        mo <- lm(y - I(1 / 60) * x ~ 1)
        tmp <- abs(residuals(mo))
        tmp_velo <- 1 / 60 # velocity...
        tmp_inter <- coef(mo)[[1]] # .. and intercept are overwritten
      }
      tmp_score <- sum(tmp) # the sum of residuals is the temporary score
      names(tmp) <-
        z # tmp, y and x are named by x to decide for outliers
      names(y) <- z
      names(x) <- z
      out <-
        sort(tmp, decreasing = TRUE)[1] # out is the sorted vector of residuals
      if (n_out >= 1) {
        # checks if more than 0 outliers are allowed
        for (i in seq_len(n_out)) {
          # the loop iterates over the number of allowed outliers
          tmp_n <- names(out) # all but the i worst are selected
          tmp_y <-
            y[!names(y) %in% tmp_n] # new y is chosen (without outliers)
          tmp_x <-
            x[!names(y) %in% tmp_n] # new x is chosen (without outliers)
          mo <- lm(tmp_y ~ tmp_x) # new linear fit
          mo_save <- mo
          tmp_velo_o <- coef(mo)[[2]] # velocity is cached
          if (coef(mo_save)[[2]] < (- (1 / 60)) &
            stran == FALSE) {
            # if the stranded option is inactive and ...
            # ...the slope is smaller - 1/60 the residuals are calculated as
            # if the slope was -1/60
            mo <- lm(tmp_y - I((-1 / 60)) * tmp_x ~ 1)
            tmp_velo_o <- (-1 / 60) # velocity is overwritten
          }
          if (coef(mo_save)[[2]] > (1 / 60)) {
            # if the slope is bigger 1/60, the residuals are calculated as
            # if the slope was 1/60
            mo <- lm(tmp_y - I(1 / 60) * tmp_x ~ 1)
            tmp_velo_o <- 1 / 60 # velocity is overwritten
          }
          tmp_inter_o <- coef(mo)[[1]] # intercept is cached
          # the sum of residuals is the temporary score with outlier
          # penalty times i
          tmp_score_o <- sum(abs(residuals(mo))) + pen * i
          tmp_tmp <- abs(residuals(mo))
          if (coef(mo_save)[[2]] < 0 &
            stran == TRUE) {
            # if the stranded option is active and ...
            # ...the slope is negative, the residuals are calculated as
            # if the slope was 0
            tmp_score_o <-
              sum(abs(tmp_y - (mean(tmp_y)))) + pen * i # the score is
            #overwritten
            tmp_velo_o <- 0 # velocity...
            tmp_inter_o <- mean(tmp_y) # .. and intercept are overwritten
            tmp_tmp <- abs(tmp_y - (mean(tmp_y)))
          }
          tmp_velo <-
            c(tmp_velo, tmp_velo_o) # velocity, intercept and score...
          tmp_inter <-
            c(tmp_inter, tmp_inter_o) # ...are put into one long vector
          tmp_score <- c(tmp_score, tmp_score_o)
          out <- c(out, sort(tmp_tmp, decreasing = TRUE)[1])
        }
      }
      mi <-
        which(tmp_score == min(tmp_score))[1] # the lowest score is chosen
      nam <- paste0(z, collapse = ",") # the name is all IDs divided by ","
      nam <- paste0(nam, "_", tmp_velo[mi]) # the velocity is pasted behind that
      #by "_"
      nam <- paste0(nam, "_", tmp_inter[mi]) # same with the intercept
      if (mi > 1) {
        # if multiple tmp_scores exist, outliers were found
        outlier <- paste(names(out)[seq_len(mi - 1)], collapse = ",")
        nam <- paste0(nam, "_", outlier) # the outliers are pasted behind the
        # name with "_"
      }
      res <- tmp_score[mi] # the final result is cached
      names(res) <- nam # and the name is added
    } else {
      # if only one value is given..
      nam <- paste0(z, collapse = ",") # ... the names are pasted together..
      nam <- paste0(nam, "_", "0") # ...the velocity is 0..
      nam <-
        paste0(nam, "_", mean(y)) # ... intercept is the mean of the input
      #values
      res <- 0 # and the score is 0 (perfect score)
      names(res) <- nam
    }
    res # the score is returned with additional info in the name
  }

# score_fun_linear_PDD scores the values of y on how close they are to a linear
# fit.
# The script also returns the respective slope of the respective fit.
# The output is "IDs_slope_intercept_outliers" as name and score as value.
# y is a vector of values, x is a vector of respective position,
# z is a vector of identifiers (IDs in our case) (is the position by default)
# stran should be TRUE if strand is not (all) NA. n_out is the number
# of allowed outliers.Spezial version for PDD

score_fun_linear_PDD <-
  function(y,
           x,
           z = x,
           pen,
           stran,
           n_out = min(10, max(1, 0.2 * length(x)))) {
    if (length(unique(x)) > 1) {
      # a linear regression can not be used on just one position
      n_out <-
        min(n_out, length(x) - 3) # the number of allowed outliers is decided
      mo <- lm(y ~ x) # the linear fit is performed
      mo_save <- mo
      tmp <- abs(residuals(mo)) # the residuals are cached
      tmp_velo <- coef(mo)[[2]] # velocity...
      tmp_inter <- coef(mo)[[1]] # .. and intercept are cached
      if (coef(mo_save)[[2]] < 0 &
        stran == TRUE) {
        # if the stranded option is active and ...
        # ...the slope is negative, the residuals are calculated as if the
        # slope was 0
        tmp <- abs(y - (mean(y))) # the score is overwritten
        tmp_velo <- 0 # velocity...
        tmp_inter <- mean(y) # .. and intercept are overwritten
      }
      tmp_score <- sum(tmp) # the sum of residuals is the temporary score
      names(tmp) <-
        z # tmp, y and x are named by x to decide for outliers
      names(y) <- z
      names(x) <- z
      out <-
        sort(tmp, decreasing = TRUE)[1] # out is the sorted vector of residuals
      if (n_out >= 1) {
        # checks if more than 0 outliers are allowed
        for (i in seq_len(n_out)) {
          # the loop iterates over the number of allowed outliers
          tmp_n <- names(out) # all but the i worst are selected
          tmp_y <-
            y[!names(y) %in% tmp_n] # new y is chosen (without outliers)
          tmp_x <-
            x[!names(y) %in% tmp_n] # new x is chosen (without outliers)
          mo <- lm(tmp_y ~ tmp_x) # new linear fit
          mo_save <- mo
          tmp_velo_o <- coef(mo)[[2]] # velocity is cached
          tmp_inter_o <- coef(mo)[[1]] # intercept is cached
          # the sum of residuals is the temporary score with outlier penalty
          # times i
          tmp_score_o <- sum(abs(residuals(mo))) + pen * i
          tmp_tmp <- abs(residuals(mo))
          if (coef(mo_save)[[2]] < 0 &
            stran == TRUE) {
            # if the stranded option is active and ...
            # ...the slope is negative, the residuals are calculated as if the
            # slope was 0
            tmp_score_o <-
              sum(abs(tmp_y - (mean(tmp_y)))) + pen * i # the score is
            #overwritten
            tmp_velo_o <- 0 # velocity...
            tmp_inter_o <- mean(tmp_y) # .. and intercept are overwritten
            tmp_tmp <- abs(tmp_y - (mean(tmp_y)))
          }
          tmp_velo <-
            c(tmp_velo, tmp_velo_o) # velocity, intercept and score...
          tmp_inter <-
            c(tmp_inter, tmp_inter_o) # ...are put into one long vector
          tmp_score <- c(tmp_score, tmp_score_o)
          out <- c(out, sort(tmp_tmp, decreasing = TRUE)[1])
        }
      }
      mi <-
        which(tmp_score == min(tmp_score))[1] # the lowest score is chosen
      nam <- paste0(z, collapse = ",") # the name is all IDs divided by ","
      nam <-
        paste0(nam, "_", tmp_velo[mi]) # the velocity is pasted behind that
      #by "_"
      nam <- paste0(nam, "_", tmp_inter[mi]) # same with the intercept
      if (mi > 1) {
        # if multiple tmp_scores exist, outliers were found
        outlier <- paste(names(out)[seq_len(mi - 1)], collapse = ",")
        nam <-
          paste0(nam, "_", outlier) # the outliers are pasted behind the
        #name with "_"
      }
      res <- tmp_score[mi] # the final result is cached
      names(res) <- nam # and the name is added
    } else {
      # if only one value is given..
      nam <- paste0(z, collapse = ",") # ... the names are pasted together..
      nam <- paste0(nam, "_", "0") # ...the velocity is 0..
      nam <-
        paste0(nam, "_", mean(y)) # ... intercept is the mean of the input
      # values
      res <- 0 # and the score is 0 (perfect score)
      names(res) <- nam
    }
    res # the score is returned with additional info in the name
  }


# score_fun_ave scores the values of y on how close they are to their
# respective mean.
# also returns the respective mean.
# the output is IDs_mean_outliers as name, and score as value.
# y is a vector of values, z is a vector of identifiers (in our case IDs)
score_fun_ave <-
  function(y, z, pen, n_out = min(10, max(1, 0.2 * length(z)))) {
    n_out <-
      min(n_out, length(z) - 3) # the number of allowed outliers is decided
    tmp_ave <- mean(y) # the mean is calculated
    tmp_dif <- abs(y - tmp_ave) # the difference to the mean is calculated
    names(tmp_dif) <-
      z # tmp_dif and y are named by z to decide for outliers
    names(y) <- z
    out <-
      sort(tmp_dif, decreasing = TRUE)[1] # out is is the sorted vector of
    #differences to the mean
    tmp_score <-
      sum(tmp_dif) # the sum of the differences to the mean is used as score
    tmp_mean <- tmp_ave # the mean is cached
    if (n_out >= 1) {
      # checks if more than 0 outliers are allowed
      for (i in seq_len(n_out)) {
        # the loop iterates over the number of allowed outliers
        tmp_out <- names(out) # all but the i worst are selected
        tmp_y <- y[!names(y) %in% tmp_out]
        tmp_ave_o <- mean(tmp_y) # new mean
        tmp_dif_o <- abs(tmp_y - tmp_ave_o) # new difference to mean
        out <- c(out, sort(tmp_dif_o, decreasing = TRUE)[1])
        tmp_score_o <-
          sum(tmp_dif_o) + pen * i # new score with outliers penalty
        tmp_score <- c(tmp_score, tmp_score_o) # the score(s) are cached
        tmp_mean <- c(tmp_mean, tmp_ave_o) # the mean(s) are cached
      }
    }
    mi <-
      which(tmp_score == min(tmp_score))[1] # the lowest score is selected
    nam <- paste0(z, collapse = ",") # all names (z) are collected
    nam <- paste0(nam, "_", tmp_mean[mi]) # the mean is pasted to the name
    if (mi > 1) {
      # checks if outliers were found
      outlier <-
        paste(names(out)[seq_len(mi - 1)], collapse = ",") # pastes all
      #outliers together
      nam <- paste0(nam, "_", outlier) # pastes the outliers to the mean
    }
    res <- tmp_score[mi] # the lowest score
    names(res) <- nam # the name
    res
  }

score_fun_increasing <- function(x, y) {
  x <-
    as.matrix(x) # x need to be given as a matrix, only relevant if only two...
  # ...values are given as an input.
  tmp <- diff(x) # the difference is calculated
  tmp_score <-
    sum(tmp) # the sum of the differences is the first part of the scoring
  tmp_x <- x[2, ] # only the second row is taken
  tmp_x[which(tmp_x < 0)] <- 0 # all values below 0 are made 0
  tmp_x[which(tmp_x == 0)] <- NA # all 0 are made NA
  tmp_x <- (log(tmp_x) + 0.5) # the ln(x)+0.5 is calculated
  tmp_x[is.na(tmp_x)] <-
    -Inf # all former 0 (now NA) ate -inf as ln(0) is NaN
  tmp_start_pen <-
    sum(tmp_x) # the sum of the ln is the second part of the scoring
  score <-
    sum(c(tmp_start_pen, tmp_score)) # the two scorings are combined
  names(score) <-
    paste0(y, collapse = ",") # the name is all given IDs divided by ,
  score # returns the score
}
