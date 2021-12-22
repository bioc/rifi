#' fragment_delay: performs the delay fragmentation.
#' fragment_delay makes delay_fragments based on position_segments and assigns
#' all gathered information to the probe based data frame.
#' The columns "delay_fragment", "velocity_fragment", "intercept" and "slope"
#' are added.
#' fragment_delay makes delay_fragments, assigns slopes, which are 1/velocity
#' at the same time, and intercepts for the TU calculation.
#' The function used is:
#' .score_fun_linear
#' the input is the probe, a data frame with ID, delay, position,
#' position_segment and strand.
#' pen is the penalty for new fragments in the dynamic programming, pen_out is
#' the outlier penalty.
#'
#' @param probe data frame: the probe based data frame.
#' @param cores cores: integer: the number of assigned cores for the task.
#' @param pen numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer fragments. Default is the auto generated value.
#' @param pen_out numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer allowed outliers. Default is the auto generated
#' value.
#' 
#' @return the probe data frame with the columns regarding the delay:
#'  delay_fragment, velocity_fragment, intercept and slope:
#' \describe{
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{position:}{The bin/probe specific position}
#'   \item{strand:}{The bin/probe specific strand}
#'   \item{position_segment:}{The position based segment}
#'   \item{delay:}{The delay value of the bin/probe}
#'   \item{delay_fragment:}{The delay fragment the bin belongs to}
#'   \item{velocity_fragment:}{The velocity value of the respective delay
#'   fragment}
#'   \item{intercept:}{The vintercept of fit through the respective delay
#'   fragment}
#'   \item{slope:}{The slope of the fit through the respective delay fragment}
#' }
#' 
#' @examples
#' data(fragmentation_minimal)
#' data(penalties_minimal)
#' fragment_delay(
#'   probe = fragmentation_minimal, cores = 2,
#'   pen = penalties_minimal["delay_penalty"],
#'   pen_out = penalties_minimal["delay_outlier_penalty"]
#' )
#' 
#' @export


fragment_delay <- function(probe, cores = 1, pen, pen_out) {
  stranded <- 1

  num_args <- list(pen, pen_out)
  names(num_args) <- c("pen", "pen_out")
  assert(
    all(unlist(lapply(
      num_args,
      FUN = function(x) {
        (is.numeric(x) &
          length(x) == 1)
      }
    ))),
    paste0(
      "'",
      names(which(unlist(
        lapply(
          num_args,
          FUN = function(x) {
            (is.numeric(x) &
              length(x) == 1)
          }
        )
      ) == FALSE))[1],
      "' must be numeric of length one or given by the logs"
    )
  )
  assert(cores > 0, "'cores' must be a positive integer")
  req_cols_probe <-
    c("ID", "position", "strand", "position_segment", "delay")
  assert(all(req_cols_probe %in% colnames(probe)),
    paste0("'", req_cols_probe[which(!req_cols_probe %in% colnames(probe))],
           "' must be a column in 'probe'!")
  )
  # I.Preparations: the dataframe is configured and some other variables are
  # assigned
  registerDoMC(cores) # cores for DoMC

  # the dataframe is sorted by strand and position.
  probe <- probe[with(probe, order(-xtfrm(probe$strand), probe$position)), ]

  probe[probe$strand == "-", ] <- probe[probe$strand == "-", ][
    order(probe[probe$strand == "-", ]$position, decreasing = TRUE), ]

  probe[, "delay_fragment"] <- NA
  probe[, "velocity_fragment"] <- NA
  probe[, "intercept"] <- NA
  probe[, "slope"] <- NA

  # tmp_df takes all relevant variables from probe
  tmp_df <-
    data.frame(
      ID = probe$ID,
      val = probe$delay,
      position = probe$position,
      seg = probe$position_segment
    )
  # if strand information is given and stranded is TRUE, all "-" are inverted...
  # ...and strand is added to tmp_df.
  if (stranded == TRUE) {
    tmp_df$strand <- probe$strand
    # the positions are inverted for the "-" strand
    # so that the slope is seen as increasing by the scoring function.
    tmp_df[tmp_df$strand == "-", "position"] <-
      ((tmp_df[tmp_df$strand == "-", "position"]) -
      (tmp_df[tmp_df$strand == "-", ][1, "position"])) * -1
  }
  # All lines with any NA are ignored at this point and are taken care of later
  # in the function.
  tmp_df <- na.omit(tmp_df)

  # makes a vector of all position segments (S_1,S_2,...)
  unique_seg <- unlist(unique(tmp_df$seg))

  # count is needed to name the fragments.
  count <- 1

  # II. Dynamic Programming: the scoring function is interpreted

  # the foreach loop iterates over each unique segment
  frags <- foreach(k = seq_along(unique_seg)) %dopar% {
    # only the part of the tmp_df that responds to the respective segment is
    # picked
    section <- tmp_df[which(tmp_df$seg == unique_seg[k]), ]

    # best_frags collects all scores that the dp is referring to
    best_frags <- c()
    # best names collects the names, and its last element is returned as result
    best_names <- c()

    if (nrow(section) > 2) {
      # only segments with more than two values are grouped...*

      for (i in 3:nrow(section)) {
        # the loop iterates over each value in the segment
        # this part always goes from position 1 to the referred position
        # 1:3,1:4...
        tmp_score <-
          score_fun_linear(section[seq_len(i), "val"],
                           section[seq_len(i), "position"],
                           section[seq_len(i), "ID"], pen_out, stranded)
        tmp_name <- names(tmp_score)
        # in this loop all smaller parts are scored e.g (i = 6) 6:6,5:6,4:6...
        # they are then combined with the former score e.g 1,2,3,4,5|6,
        # 1,2,3,4|5,6...
        if (i > 5) {
          # only parts bigger than 6 are accepted as three is the smallest
          # possible fragment size
          for (j in (i - 2):4) {
            tmp_val <- section[j:i, "val"]
            tmp_position <- section[j:i, "position"]
            tmp_ID <- section[j:i, "ID"]
            # penalty for a new fragment and former scores are added
            tmp <-
              score_fun_linear(tmp_val, tmp_position, tmp_ID, pen_out,
                               stranded) + pen + best_frags[j - 3]
            tmp_score <- c(tmp_score, tmp)
            # the new fragment is pasted to its corresponding former fragment
            tmp_n <- paste0(best_names[j - 3], "|", names(tmp))
            tmp_name <- c(tmp_name, tmp_n) # the names is cached
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
    } else {
      #* ...all segments with less than three values are grouped automatically
      tmp_score <-
        score_fun_linear(section[, "val"], section[, "position"],
                         section[, "ID"], pen_out, stranded)
      tmp_name <- names(tmp_score)
      best_names <- c(best_names, tmp_name)
    }
    # the final result put into a list called frags
    best_names[length(best_names)]
  }

  # III. Fill the dataframe

  # this loop iterates over the segments to fill the dataframe
  for (k in seq_along(frags)) {
    # the single fragments are split by |
    na <- strsplit(frags[[k]], "\\|")[[1]]

    # the loop goes over each fragment
    for (i in seq_along(na)) {
      # trgt are all IDs in the fragment
      tmp_trgt <- strsplit(na[i], "_")[[1]][1]
      trgt <- strsplit(tmp_trgt, ",")[[1]]
      # This statement checks if there are outliers
      if (length(strsplit(na[i], "_")[[1]]) == 4) {
        # outl are the IDs of all outliers
        tmp_outl <- strsplit(na[i], "_")[[1]][4]
        outl <- strsplit(tmp_outl, ",")[[1]]
        # while the original trgt gets stripped of the outliers, so trgt are
        # now all valid IDs
        trgt <- trgt[-which(trgt %in% outl)]
        # now the real outliers are named
        rows <- match(outl, probe[, "ID"])
        nam <- paste0("D_", count, "_O")
        probe[rows, "delay_fragment"] <- nam
        # and the values are assigned (normal outliers are allowed to share...
        # ...the same values, while terminal outliers are not)
        probe[rows, "slope"] <-
          as.numeric(strsplit(na[i], "_")[[1]][2])
        probe[rows, "intercept"] <-
          as.numeric(strsplit(na[i], "_")[[1]][3])
        probe[rows, "velocity_fragment"] <-
          1 / (as.numeric(strsplit(na[i], "_")[[1]][2]))
      }
      # last thing is that the real fragment is named and the values are...
      # ...assigned the same way as for the outliers
      rows <- match(trgt, probe[, "ID"])
      nam <- paste0("D_", count)
      probe[rows, "delay_fragment"] <- nam
      probe[rows, "slope"] <- as.numeric(strsplit(na[i], "_")[[1]][2])
      probe[rows, "intercept"] <-
        as.numeric(strsplit(na[i], "_")[[1]][3])
      probe[rows, "velocity_fragment"] <-
        1 / (as.numeric(strsplit(na[i], "_")[[1]][2]))
      # the count increases by one for the next iteration
      count <- count + 1
    }
  }
  # in this part the NAs are dealt with initial check if the first lines are NA
  if (sum(cumprod(is.na(probe$delay_fragment))) > 0) {
    probe$delay_fragment[seq_len(sum(cumprod(is.na(probe$delay_fragment))))] <-
      "D_0.5_NA"
  }
  # all rows with NA in delay fragment are collected
  row_NA <- which(is.na(probe$delay_fragment))
  # check if there are any NAs
  if (length(row_NA) > 0) {
    # group is initiated with the first NA
    group <- c(row_NA[1])
    # the loop iterates over all NAs
    for (i in seq_along(row_NA)) {
      # if there are more NAs following the first NA, they will be put
      # together in group
      if (is.na(probe$delay_fragment[row_NA[i] + 1])) {
        group <- c(group, row_NA[i] + 1)
      }
      # if there are no more NAs to add, it will be checked if the fragment
      # above and below is the same or not. Normal outliers are ignored in
      # this case
      else if (gsub("_O", "", probe$delay_fragment[group[1] - 1]) ==
               gsub("_O", "", probe$delay_fragment[group[length(group)] + 1])) {
        # if they are the same, so the NAs are in the middle of one fragment
        # they get the name marked with an _NA and the values of the fragment
        probe$delay_fragment[group] <-
          paste0(gsub("_O", "", probe$delay_fragment[group[1] - 1]), "_NA")
        probe$velocity_fragment[group] <-
          probe$velocity_fragment[group[1] - 1]
        probe$intercept[group] <- probe$intercept[group[1] - 1]
        probe$slope[group] <- probe$slope[group[1] - 1]
        # group is initiated with the next NA after the ones that were just
        # dealt with
        group <- row_NA[i + 1]
      } else if (probe$delay_fragment[group[1] - 1] !=
                 probe$delay_fragment[group[length(group)] + 1]) {
        # if the NAs are between two fragments, it makes a new fragment
        # between the two
        probe$delay_fragment[group] <-
          paste0("D_", as.numeric(gsub(
            "D_|_O", "", probe$delay_fragment[group[1] - 1]
          )) + 0.5, "_NA")
        # group is initiated with the next NA after the ones that were just
        # dealt with
        group <- row_NA[i + 1]
      }
    }
    # this deals with the case, that NAs might be ate the very end of the
    # dataframe
    probe$delay_fragment[is.na(probe$delay_fragment)] <-
      paste0("D_", as.numeric(gsub("D_|_O", "",
               probe$delay_fragment[!is.na(probe$delay_fragment)]
               [length(probe$delay_fragment[!is.na(probe$delay_fragment)])])) +
               0.5, "_NA")
  }
  # the function ends with the return of the dataframe
  probe
}
