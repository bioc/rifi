fragment_inty_pen <- function(probe, pen, pen_out, from, to, cores) {
  # I.Preparations: the dataframe is configured and some other variables are
  # assigned

  registerDoMC(cores)

  probe <- probe[with(probe, order(-xtfrm(probe$strand), probe$position)), ]

  probe[probe$strand == "-", ] <- probe[probe$strand == "-", ][
    order(probe[probe$strand == "-", ]$position, decreasing = TRUE), ]

  # In this special version we look at position_segments
  tmp_df <-
    data.frame(
      ID = probe$ID,
      val = probe$intensity,
      seg = probe$position_segment
    )

  # this is just for safety, nothing should be omitted here
  tmp_df <- na.omit(tmp_df)

  # take log2 to adress noise
  tmp_df$val <- log2(tmp_df$val)

  unique_seg <- unlist(unique(tmp_df$seg))

  # II. Dynamic Programming: the scoring function is interpreted

  frags <- foreach(k = seq_along(unique_seg)) %dopar% {
    section <- tmp_df[which(tmp_df$seg == unique_seg[k]), ]

    best_frags <- c()
    best_names <- c()

    if (nrow(section) >= from &
      nrow(section) <= to) {
      # the sampling is between 10 an 200

      for (i in 2:nrow(section)) {
        tmp_score <-
          score_fun_ave(section[seq_len(i), "val"],
                        section[seq_len(i), "ID"], pen_out)
        tmp_name <- names(tmp_score)
        if (i > 3) {
          for (j in (i - 1):3) {
            tmp_val <- section[j:i, "val"]
            tmp_ID <- section[j:i, "ID"]
            tmp <-
              score_fun_ave(tmp_val, tmp_ID, pen_out) + pen + best_frags[j - 2]
            tmp_score <- c(tmp_score, tmp)
            tmp_n <- paste0(best_names[j - 2], "|", names(tmp))
            tmp_name <- c(tmp_name, tmp_n)
          }
        }
        pos <- which(tmp_score == min(tmp_score))[1]
        tmp_score <- tmp_score[pos]
        tmp_name <- tmp_name[pos]
        best_frags <- c(best_frags, tmp_score)
        best_names <- c(best_names, tmp_name)
      }
      best_names[length(best_names)]
    }
  }

  # III. Fill the dataframe

  frags <- frags[!(vapply(frags, is.null, logical(1)))]

  if (length(frags) == 0) {
    res1 <- 0
    res2 <- 0
    res <- list(res1, res2)
    names(res) <- c("intensity", "intensity_out")
    return(res)
  }

  comp <- foreach(k = seq_along(frags)) %dopar% {
    na <- strsplit(frags[[k]], "\\|")[[1]]

    parts <- vector("list", length = length(na))

    for (i in seq_along(na)) {
      tmp_trgt <- strsplit(na[i], "_")[[1]][1]
      trgt <- strsplit(tmp_trgt, ",")[[1]]
      if (length(strsplit(na[i], "_")[[1]]) == 3) {
        tmp_outl <- strsplit(na[i], "_")[[1]][3]
        outl <- strsplit(tmp_outl, ",")[[1]]
        trgt <- trgt[-which(trgt %in% outl)]
      }
      rows <- match(trgt, probe[, "ID"])

      parts[[i]] <- probe[rows, "intensity"]
    }
    parts
  }

  p <- c()
  all <- 0

  for (k in seq_along(comp)) {
    all <- all + (length(comp[[k]]) - 1)
    if (length(comp[[k]]) > 1) {
      for (i in seq_len(length(comp[[k]]) - 1)) {
        if (!identical(unique(comp[[k]][[i]]), unique(comp[[k]][[i + 1]]))) {
          tryCatch({
              p_val <- t.test(comp[[k]][[i]], comp[[k]][[i + 1]])$p.value
              p <- c(p, p_val)
            },
            error = function(e) {
            }
          )
        }
      }
    }
  }
  p <- p.adjust(p)
  count <- length(p[which(p < 0.01)])
  res1 <- count
  res2 <- all - count
  res <- list(res1, res2)
  names(res) <- c("intensity", "intensity_out")
  res
}
