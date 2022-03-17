plot_nls2_function <-
  function(inp) {
    inp <- inp_order(inp)
    assay(inp)[decode_FLT(inp)] <- NA
    time <- metadata(inp)$timepoints
    fit_STD <- metadata(inp)$fit_STD
    fit_TI <- metadata(inp)$fit_TI
    #assign a data frame
    pdf("fit_nls2.pdf")
    par(mfrow = c(2, 2))
    for (i in seq_len(nrow(inp))) {
      tmp_inp <- inp[i,]
      ID <- rowRanges(tmp_inp)$ID
      f<-data.frame(bg = 0)
      if(ID %in% fit_STD$ID){
        f <- fit_STD[fit_STD$ID == ID,]
      }
      if(ID %in% fit_TI$ID){
        f <- fit_TI[fit_TI$ID == ID,]
      }
      if(!all(is.na(assay(tmp_inp)))){
        plot(time, assay(tmp_inp), pch = 16, xlab = "Time [min]",
             ylab = "Intensity [A.U.]",
             main = paste0("ID: ", ID,
                           " position: ", rowRanges(tmp_inp)$position,
                           " ", decode(strand(tmp_inp))
             ),
             ylim = c(0.9*f$bg, max(max(assay(tmp_inp), na.rm = TRUE),0)),
             col = colData(inp)$replicate + 1
        )
        mean_r <- tapply(as.numeric(assay(tmp_inp)), colData(inp)$timepoint,
                         mean, na.rm = TRUE)
        points(as.numeric(names(mean_r)), mean_r, pch = 16, col = 1)
      }
      if(ID %in% fit_STD$ID){
        f <- fit_STD[fit_STD$ID == ID,]
        if(!any(is.na(f))){
          curve((f$k / f$decay * x / x + f$bg),
                from = 0,
                to = f$delay,
                type = "l",
                add = TRUE,
                pch = 7
          )
          curve(((f$k / f$decay) * exp(-f$decay * (x - f$delay)) + f$bg),
                from = f$delay,
                to = max(time),
                type = "l",
                add = TRUE,
                pch = 7
          )
          legend(
            "topright",
            legend = c(
              paste0("delay = ", round(f$delay, digits = 1)),
              paste0("halflife = ", round(log(2) / f$decay, digits = 1)),
              paste0("background = ", round(f$bg, digits = 1))
            ),
            bty = "n",
            cex = 0.8
          )
        }
      }
      if(ID %in% fit_TI$ID){
        f <- fit_TI[fit_TI$ID == ID,]
        if(!any(is.na(f))){
          curve(((f$k / f$decay - f$ti / f$decay) * x / x + f$bg),
                from = 0,
                to = f$ti_delay,
                type = "l",
                add = TRUE,
                pch = 7
          )
          curve(f$k / f$decay + f$bg - f$ti / f$decay *
                  exp(-f$decay * (x - f$ti_delay)),
                from = f$ti_delay,
                to = f$ti_delay + f$rest_delay,
                type = "l",
                add = TRUE,
                pch = 7
          )
          curve(
            f$bg + (f$k / f$decay - f$ti / f$decay * exp(
              -f$decay * f$rest_delay)) * exp(-f$decay * (
                x - (f$ti_delay + f$rest_delay))),
            from = f$ti_delay + f$rest_delay,
            to = max(time),
            type = "l",
            add = TRUE,
            pch = 7
          )
          legend(
            "topright",
            legend = c(
              paste0("ti_delay = ", round(f$ti_delay, digits = 1)),
              paste0("delay = ", round(f$ti_delay + f$rest_delay, digits = 1)),
              paste0("halflife = ", round(log(2) / f$decay, digits = 1)),
              paste0("ti = ", round(f$ti, digits = 1)),
              paste0("term_prob = ", round(f$ti / f$k, digits = 1)),
              paste0("background = ", round(f$bg, digits = 1))
            ),
            bty = "n",
            cex = 0.8
          )
        }
      }
    }
    dev.off()
  }
