fragment_delay_pen <- function(inp, pen, pen_out, from, to, cores) {
  stranded <- 1

  # I. Preparations: the dataframe is configured and some other variables are
  # assigned

  registerDoMC(cores) # cores for DoMC

  # the dataframe is sorted by strand and position.
  inp <- inp_order(inp)
  
  tmp_df <- inp_df(inp, "ID", "delay", "position_segment")
  
  tmp_df <- tmp_df_rev(tmp_df, "-")

  # All lines with any NA are ignored at this point and are taken care of later
  # in the function.
  tmp_df <- na.omit(tmp_df)
  # makes a vector of all position segments (S_1,S_2,...)
  unique_seg <- unlist(unique(tmp_df$position_segment))
  # II. Dynamic Programming: the scoring function is interpreted
  # the foreach loop iterates over each unique segment
  frags <- foreach(k = seq_along(unique_seg)) %dopar% {
    # only the part of the tmp_df that responds to the respective segment is
    # picked
    section <- tmp_df[which(tmp_df$position_segment == unique_seg[k]), ]

    # best_frags collects all scores that the dp is referring to
    best_frags <- c()
    # best names collects the names, and its last element is returned as result
    best_names <- c()

    if (nrow(section) >= from &
      nrow(section) <= to) {
      # the sampling is between 10 an 200

      for (i in 3:nrow(section)) {
        # the loop iterates over each value in the segment
        # this part always goes from position 1 to the referred position
        # 1:3,1:4...
        tmp_score <-
          score_fun_linear(section[seq_len(i), "delay"],
                           section[seq_len(i), "position"],
                           section[seq_len(i), "ID"], pen_out, stranded)
        tmp_name <- names(tmp_score)
        # in this loop all smaller parts are scored eg (i = 6) 6:6,5:6,4:6...
        # they are then combined with the former score e.g 1,2,3,4,5|6,
        # 1,2,3,4|5,6...
        if (i > 5) {
          # only parts bigger than 6 are accepted as three is the smallest
          # possible fragment size
          for (j in (i - 2):4) {
            tmp_val <- section[j:i, "delay"]
            tmp_position <- section[j:i, "position"]
            tmp_ID <- section[j:i, "ID"]
            # penalty for a new fragment and former scores are added
            tmp <-
              score_fun_linear(tmp_val, tmp_position, tmp_ID, pen_out,
                               stranded) + pen +
              best_frags[j - 3]
            tmp_score <- c(tmp_score, tmp)
            # the new fragment is pasted to its corresponding former fragment
            tmp_n <- paste0(best_names[j - 3], "|", names(tmp))
            tmp_name <- c(tmp_name, tmp_n) # the names are cached
          }
        }
        # from the first score eg 1:6 and the smaller scores from the loop
        # 1,2,3,4,5|6, 1,2,3,4|5,6... the smallest is chosen and passed to
        # best_frags and best_names for the next iteration
        pos <-
          which(tmp_score == min(tmp_score))[1] # lowest score is collected
        tmp_score <- tmp_score[pos]
        tmp_name <- tmp_name[pos]
        best_frags <- c(best_frags, tmp_score)
        best_names <- c(best_names, tmp_name)
      }
    }
    # the final result put into a list called frags
    best_names[length(best_names)]
  }

  frags <- frags[!(vapply(frags, is.null, logical(1)))]

  if (length(frags) == 0) {
    res1 <- 0
    res2 <- 0
    res <- list(res1, res2)
    names(res) <- c("delay", "delay_out")
    return(res)
  }

  # III. Fill the dataframe

  # this loop iterates over the segments to fill the dataframe
  comp <- foreach(k = seq_along(frags)) %dopar% {
    na <- strsplit(frags[[k]], "\\|")[[1]]
    group <- c()
    position <- c()
    delay <- c()
    # for the first one
    A <- 0
    B <- 0
    # for all between
    if (length(na) > 1) {
      for (i in seq_len(length(na) - 1)) {
        tmp_trgt <- strsplit(na[i], "_")[[1]][1]
        trgt <- strsplit(tmp_trgt, ",")[[1]]
        if (length(strsplit(na[i], "_")[[1]]) == 4) {
          tmp_outl <- strsplit(na[i], "_")[[1]][4]
          outl <- strsplit(tmp_outl, ",")[[1]]
          trgt <- trgt[-which(trgt %in% outl)]
        }
        rows <- match(trgt, rowRanges(inp)$ID)
        tmp_group <- rep(as.character(i), length(rows))
        tmp_position <- rowRanges(inp)$position[rows]
        stra <- unique(decode(strand(inp))[rows])
        if (stranded == TRUE & stra == "-") {
          tmp_position <- (tmp_position - (tmp_df[tmp_df$strand == "-", ][
              nrow(tmp_df[tmp_df$strand == "-", ]), "position"])) * -1
        }
        tmp_position <- tmp_position - min(tmp_position, na.rm = TRUE)
        tmp_delay <- rowRanges(inp)$delay[rows] - A + B
        sta_point <-
          as.numeric(strsplit(na[i], "_")[[1]][2]) *
          rowRanges(inp)$position[rows[1]] + as.numeric(strsplit(na[i], "_")[[1]][3])
        end_point <-
          as.numeric(strsplit(na[i], "_")[[1]][2]) *
          rowRanges(inp)$position[rows[length(rows)]] +
          as.numeric(strsplit(na[i], "_")[[1]][3])
        new_point <-
          as.numeric(strsplit(na[i + 1], "_")[[1]][2]) *
          rowRanges(inp)$position[rows[length(rows)]] +
          as.numeric(strsplit(na[i + 1], "_")[[1]][3])
        A <- sum(A, (new_point - sta_point))
        B <- sum(B, (new_point - end_point))
        group <- c(group, tmp_group)
        position <- c(position, tmp_position)
        delay <- c(delay, tmp_delay)
      }
    }
    # for the last one
    tmp_trgt <- strsplit(na[length(na)], "_")[[1]][1]
    trgt <- strsplit(tmp_trgt, ",")[[1]]
    if (length(strsplit(na[length(na)], "_")[[1]]) == 4) {
      tmp_outl <- strsplit(na[length(na)], "_")[[1]][4]
      outl <- strsplit(tmp_outl, ",")[[1]]
      trgt <- trgt[-which(trgt %in% outl)]
    }
    rows <- match(trgt, rowRanges(inp)$ID)
    tmp_group <- rep(as.character(length(na)), length(rows))
    tmp_position <- rowRanges(inp)$position[rows]
    stra <- unique(decode(strand(inp))[rows])
    if (stranded == TRUE & stra == "-") {
      tmp_position <- (tmp_position -
        (tmp_df[tmp_df$strand == "-", ][
          nrow(tmp_df[tmp_df$strand == "-", ]), "position"])) * -1
    }
    tmp_position <- tmp_position - min(tmp_position, na.rm = TRUE)
    tmp_delay <- rowRanges(inp)$delay[rows] - A + B
    group <- c(group, tmp_group)
    position <- c(position, tmp_position)
    delay <- c(delay, tmp_delay)

    parts <- data.frame(
      group = group,
      position = position,
      delay = delay
    )
    parts
  }

  p <- c()
  all <- 0

  for (k in seq_along(comp)) {
    unique_group <- unique(comp[[k]]$group)
    all <- all + (length(unique_group) - 1)
    if (length(unique_group) > 1) {
      for (i in seq_len(length(unique_group) - 1)) {
        tmp_data <-
          comp[[k]][c(
            which(comp[[k]]$group == unique_group[i]),
            which(comp[[k]]$group == unique_group[i + 1])
          ), ]

        model.1 <- lm(delay ~ position + group + position:group,
                      data = tmp_data)
        tryCatch({
            p_val <- Anova(model.1, type = "II")$"Pr(>F)"[3]
            names(p_val) <- paste0(k, ",", i)
            p <- c(p, p_val)
          },
          error = function(e) {
          }
        )
      }
    }
  }

  p <- p.adjust(p)
  count <- length(p[which(p < 0.01)])
  rest <- names(p[which(p >= 0.01)])
  p <- c()
  for (j in rest) {
    k <- as.numeric(strsplit(j, ",")[[1]][1])
    i <- as.numeric(strsplit(j, ",")[[1]][2])
    unique_group <- unique(comp[[k]]$group)

    tmp_data <-
      comp[[k]][c(
        which(comp[[k]]$group == unique_group[i]),
        which(comp[[k]]$group == unique_group[i + 1])
      ), ]
    tryCatch({
        model.2 <- lm(delay ~ position + group, data = tmp_data)
        p_val <- Anova(model.2, type = "II")$"Pr(>F)"[2]
        p <- c(p, p_val)
      },
      error = function(e) {
      }
    )
  }
  p <- p.adjust(p)
  count <- count + length(p[which(p < 0.01)])
  res1 <- count
  res2 <- all - count
  res <- list(res1, res2)
  names(res) <- c("delay", "delay_out")
  res
}
