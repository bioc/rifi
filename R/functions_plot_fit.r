plot_nls2_function <-
  function(data,
           inp,
           color = c(
             "blue",
             "green",
             "yellow",
             "grey",
             "orange",
             "cyan",
             "pink",
             "grey64",
             "grey90",
             "black",
             "grey45"
           )) {
    time <-
      as.numeric(colnames(inp)[seq_len(which(colnames(inp) %in% "ID") - 1)])
    #assign a data frame
    Data_fit_function <- function(data) {
      #intensities values are transformed to generate the dataframe.
      inty <-
        as.vector(t(data[seq_len(which(colnames(data) %in% "ID") - 1)]))
      #time vector is adapted to the length of number of replicates and number
      # of intensities time point.
      Time <- rep(time, times = length(inty) / length(time))
      indice <-
        rep(seq_len(length(inty) / length(time)), each = length(time))
      Data_fit <-
        data.frame(as.vector(t(Time)), as.vector(t(inty)), as.vector(t(indice)))
      #NAs are eliminated from the dataframe.
      Data_fit <- na.omit(Data_fit)
      #intensities values are normalized to 1
      Data_fit[, 2] <- relativity(Data_fit[, 2])
      colnames(Data_fit) <- c("time", "inty", "indice")
      return(Data_fit)
    }
    pdf("fit_nls2.pdf")
    par(mfrow = c(3, 4))
    for (i in seq_len(nrow(data))) {
      tmp <- inp[which(inp$ID %in% data$ID[i]), ]
      Data_fit <- Data_fit_function(tmp)
      pos <- which(data$ID %in% tmp$ID[1])[1]
      decay <- log(2) / data[pos, "half_life"]
      tryCatch({
        #plot first the time and intensity for each probe/bin. The coefficients
        #are applied using the curve function in R.
        #to make sure all data points are plotted because some probes can miss 
        #one or more time point, a summary of intensity mean is calculated. 
        Data_fit_tmp <- as.data.frame(Data_fit[, c(1, 2)] %>%
                          group_by(time) %>%
                            summarise_each(funs(mean(inty))))
        ind <- unique(Data_fit$indice)
        plot(
          Data_fit_tmp[Data_fit$indice == ind[1], "time"],
          Data_fit_tmp[Data_fit$indice == ind[1], "inty"],
          pch = 16,
          xlab = "Time [min]",
          ylab = "Intensity [A.U.]",
          main = paste0("ID: ", tmp$ID[1], " ", tmp$strand[1]),
          ylim = c(min(Data_fit$inty), 1),
          xaxt = "n"
        )
        #plot the replicate in different colors
        if (length(ind) > 1) {
          for (k in seq_along(ind)) {
            points(
              Data_fit[Data_fit$indice == ind[k], "time"],
              Data_fit[Data_fit$indice == ind[k], "inty"],
              type = "p",
              xlab = "Time [min]",
              ylab = "Intensity [A.U.]",
              col = color[k]
            )
          }
        }
        if (!is.na(data$delay[pos])) {
          if (!is.na(data$intyf[pos])) {
            axis(
              side = 1,
              cex.axis = 1,
              at = c(
                round(as.numeric(
                gsub("\\,", ".", data$delay[pos])), digits = 2),
                time[length(time)-1])
            )
            curve((x < data$delay[pos]) *
                    data$inty_S0[pos] + (x >= data$delay[pos]) *
                    (data$intyf[pos] + (data$inty_S0[pos] - data$intyf[pos]) *
                       (exp(-decay * (x - data$delay[pos])))),
                  from = time[1],
                  to = data$delay[pos],
                  type = "l",
                  add = TRUE,
                  col = color[1],
                  pch = 7,
                  lty = 2
            )
            curve((x < data$delay[pos]) * data$inty_S0[pos] +
                    (x >= data$delay[pos]) *
                    (data$intyf[pos] +
                       (data$inty_S0[pos] - data$intyf[pos]) *
                       (exp(-decay * (x - data$delay[pos])))),
                  from = data$delay[pos],
                  to = 64,
                  type = "l",
                  add = TRUE,
                  col = color[8],
                  pch = 3,
                  lty = 1
            )
          } else{
            axis(
              side = 1,
              cex.axis = 1,
              at = c(
                round(as.numeric(
                  gsub("\\,", ".", data$delay[pos])), digits = 2),
                time[length(time)-1])
            )
            curve((x < data$delay[pos]) *
                    data$inty_S0[pos] + (x >= data$delay[pos]) *
                    (data$inty_S0[pos] * (exp(-decay * (x - data$delay[pos])
                    ))),
                  from = time[1],
                  to = data$delay[pos],
                  type = "l",
                  add = TRUE,
                  col = color[1],
                  pch = 7,
                  lty = 2
            )
            curve((x < data$delay[pos]) *
                    data$inty_S0[pos] + (x >= data$delay[pos]) *
                    (data$inty_S0[pos] * (exp(-decay * (x - data$delay[pos])))),
                  from = data$delay[pos],
                  to = last(time),
                  type = "l",
                  add = TRUE,
                  col = color[8],
                  pch = 3,
                  lty = 1
            )
          }
          legend(
            "topright",
            legend = paste("HL: ", formatC(
              data$half_life[pos],
              digits = 3,
              format = "fg"
            )),
            text.col = 4,
            bty = "n",
            cex = 0.8
          )
          abline(h = data$inty_S0[pos],
                 col = "grey",
                 lty = 3)
          abline(h = data$intyf[pos],
                 col = "grey",
                 lty = 2)
        }
      }, error = function(e) {
        e
      })
    }
    dev.off()
  }


plot_singleProbe_function <-
  function(data,
           inp,
           color = c(
             "blue",
             "green",
             "yellow",
             "grey",
             "orange",
             "cyan",
             "pink",
             "grey64",
             "grey90",
             "black",
             "grey45"
           )) {
    if (nrow(data) == 0) {
      return()
    }
    time <-
      as.numeric(colnames(inp)[seq_len(which(colnames(inp) %in% "ID") - 1)])

    # fit_nls2 is saved as df object to make it easier and shorter for the
    # plot
    unique_ID <- unique(inp$ID)
    #assign a data frame, the mean of each replicate is calculated
    tmp_df <-
      foreach(i = seq_along(unique_ID), .combine = rbind) %dopar% {
        tmp <- inp[which(inp$ID %in% unique_ID[i]), ]
        tmp[which(grepl("FLT", tmp$filtration)),
            seq_len(which(colnames(inp) %in% "ID") - 1)] <- NA
        tmp_mean <-
          t(colMeans(tmp[seq_len(which(colnames(inp) %in% "ID") - 1)],
                     na.rm = TRUE))
        tmp <-
          cbind(tmp_mean, tmp[1, (which(colnames(inp) %in% "ID")):
                                (which(colnames(inp) %in% "strand"))])
        tmp
      }
    tmp_df <- na.omit(tmp_df)
    ti_ids <- data[, "ID"]
    tmp_df <- tmp_df[which(tmp_df$ID %in% ti_ids), ]
    pdf("ti_fit.pdf")
    par(mfrow = c(2, 2))
    for (j in seq_along(ti_ids)) {
      print(j)
      gr <- ti_ids[j]
      temp_data <- tmp_df[which(tmp_df$ID == gr), ]
      temp_points <- inp[which(inp$ID %in% gr),
            seq_len(which(colnames(inp) %in% "ID") - 1)]
      inty <- na.omit(as.vector(t(
          temp_data[, seq_len(which(colnames(temp_data) %in% "ID") - 1)]))) /
        temp_points[1, 1]
      temp_points <- temp_points / temp_points[1, 1]
      temp_data <- cbind.data.frame(time, inty)
      if (is.finite(max(temp_data[, 2], na.rm = TRUE))) {
        plot(
          temp_data$time,
          temp_data$inty,
          col = color[1],
          pch = 1,
          type = "p",
          xlab = "Time [min]",
          ylab = "Intensity [A.U.]",
          ylim = c(0, max(temp_points, na.rm = TRUE)),
          main = paste0("ID: ", tmp_df$ID[j], " ", tmp_df$strand[j])
        )
        for (i in seq_len(nrow(temp_points))) {
          points(
            temp_data$time,
            temp_points[i, ],
            col = i,
            pch = 1,
            type = "p"
          )
        }
        points(
          temp_data$time,
          temp_data$inty,
          col = color[1],
          pch = 1,
          type = "p"
        )
        k <- data[j, "synthesis_rate"]
        decay <- log(2) / data[j, "half_life"]
        ti <- data[j, "ti_value"]
        bg <- data [j, "TI_background"]
        ti_delay <- data[j, "ti_delay"]
        rest_delay <- data[j, "delay"] - ti_delay
        term_prob <- ti / k
        curve((k / decay - ti / decay * x / x + bg),
              from = 0,
              to = ti_delay,
              add = TRUE,
              type = "l",
              lwd = 1,
              col = color[1],
              lty = 1
        )
        curve(
          k / decay + bg - ti / decay * exp(-decay * (x - ti_delay)),
          from = ti_delay,
          to = ti_delay + rest_delay,
          add = TRUE,
          type = "l",
          lwd = 1,
          col = color[1],
          lty = 1
        )
        curve(
          I(bg + (k / decay - ti / decay * exp(
            -decay * (rest_delay)
          )) * exp(-decay * (
            x - (ti_delay + rest_delay)
          ))),
          from = ti_delay + rest_delay,
          to = 64,
          lty = 2,
          add = TRUE,
          type = "l",
          lwd = 1,
          col = color[1]
        )
        legend(
          "topright",
          legend = c(
            paste0("pos=", unique(tmp_df$position[j])),
            paste0("ti_delay= ", round(ti_delay, digits = 1)),
            paste0("delay= ", round(ti_delay + rest_delay, digits = 1)),
            paste0("halflife= ", round(log(2) / decay, digits = 1)),
            paste0("ti= ", round(ti, digits = 1)),
            paste0("term_prob= ", round(term_prob, digits = 1))
          ),
          bty = "n",
          cex = 0.8,
          text.col = color[seq_len(6)]
        )
      }
    }
    dev.off()
  }
