#' viz_pen_obj: visualizes penalty objects.
#' An optional visualization of any penalty object created by make_pen.
#' Can be customized to show only the n = top_i top results.
#'
#' @param obj object: penalty object(make_pen output)
#' @param top_i integer: the number of top results visualized. Default is all.
#'
#' @return A visualization of the penalty object
#' 
#' @examples 
#' data(penalties_e_coli)
#' viz_pen_obj(penalties_e_coli$pen_obj_delay,25)
#' 
#' @export
#'
viz_pen_obj <-
  function(obj, top_i = nrow(obj[[3]][[1]]) * ncol(obj[[3]][[1]])) {
    for (i in seq_along(obj[[3]])) {
      name <- paste0("pen_plot_", names(obj[[2]]), i, ".pdf")
      pdf(file = name, width = 10)
      wrong3 <- obj[[4]][[i]]
      correct3 <- obj[[3]][[i]]
      diff <- correct3 - wrong3
      diff[is.na(diff)] <- min(diff, na.rm = TRUE)
      top <- max(rank(diff, ties.method = "last")) - top_i
      col <-
        matrix(rep(seq_len(nrow(diff)), ncol(diff)), nrow(diff), ncol(diff))
      col2 <-
        matrix(rep(as.numeric(colnames(diff)), each = nrow(diff)),
               ncol(diff), nrow(diff))
      ylim <- c(
        min(diff[which(rank(diff, ties.method = "last") >= top)],
            wrong3[which(rank(diff, ties.method = "last") >=
          top)], correct3[which(rank(diff, ties.method = "last") >= top)],
          na.rm = TRUE),
        max(diff[which(rank(diff, ties.method = "last") >= top)],
            wrong3[which(rank(diff, ties.method = "last") >=
          top)], correct3[which(rank(diff, ties.method = "last") >= top)],
          na.rm = TRUE)
      )
      xlim <-
        c(min((rank(
          diff,
          ties.method = "last"
        )[which(rank(diff, ties.method = "last") >= top)]), na.rm = TRUE),
        max((rank(
          diff,
          ties.method = "last"
        )[which(rank(diff, ties.method = "last") >= top)]), na.rm = TRUE)
        + (diff(c(
            min((rank(
              diff,
              ties.method = "last"
            )[which(rank(diff, ties.method = "last") >= top)]), na.rm = TRUE),
            max((rank(
              diff,
              ties.method = "last"
            )[which(rank(diff, ties.method = "last") >= top)]), na.rm = TRUE)
          )) * 0.1))
      fivety <- quantile(diff, 0.5)
      ninety <- quantile(diff, 0.9)
      ninetynine <- quantile(diff, 0.99)
      curve(
        0 * x + fivety,
        from = xlim[1] + diff(xlim) * 0.02,
        to = xlim[2] - diff(xlim) * 0.1,
        col = "lightgrey",
        ylim = ylim,
        xlim = xlim,
        main = paste0("top ", max(rank(
          diff,
          ties.method = "last"
        )) - top),
        xlab = "Rank",
        ylab = "Number"
      )
      curve(
        0 * x + ninety,
        add = TRUE,
        from = xlim[1] + diff(xlim) * 0.02,
        to = xlim[2] - diff(xlim) * 0.1,
        col = "lightgrey"
      )
      curve(
        0 * x + ninetynine,
        add = TRUE,
        from = xlim[1] + diff(xlim) * 0.02,
        to = xlim[2] - diff(xlim) * 0.1,
        col = "lightgrey"
      )
      text(xlim[1], fivety, "50%", cex = 0.5, col = "lightgrey")
      text(xlim[1], ninety, "90%", cex = 0.5, col = "lightgrey")
      text(xlim[1], ninetynine, "99%", cex = 0.5, col = "lightgrey")
      points(rank(diff, ties.method = "last")[
        which(rank(diff, ties.method = "last") >=  top)],
      wrong3[which(rank(diff, ties.method = "last") >= top)],
      col = "red",
      pch = 19
      )
      points(rank(diff, ties.method = "last")[
        which(rank(diff, ties.method = "last") >= top)],
      correct3[which(rank(diff, ties.method = "last") >= top)],
      col = "green",
      pch = 19
      )
      points(rank(diff, ties.method = "last")[
        which(rank(diff, ties.method = "last") >= top)],
      diff[which(rank(diff, ties.method = "last") >= top)],
      col = col[which(rank(diff, ties.method = "last") >= top)],
      pch = 19
      )
      text(
        rank(diff, ties.method = "last")[
          which(rank(diff, ties.method = "last") >= top)],
        diff[which(rank(diff, ties.method = "last") >= top)] +
          (ylim[2] * 0.02),
        label = as.character(col2[which(rank(diff, ties.method = "last")
                                        >= top)]),
        cex = 0.65,
        col = "darkgrey"
      )
      legend("right", c(row.names(diff)), fill = c(seq_len(nrow(diff))))
      dev.off()
    }
  }
