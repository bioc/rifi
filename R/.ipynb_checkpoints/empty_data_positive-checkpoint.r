#in case the dataframe from the positive strand is empty, a dataframe is
#created to be plotted for the genome annotation.
empty_data_positive <- function(data_p, data_n, frag, i, 
                                axis_title_y_size = axis_title_y_size, 
                                axis_text_y_size = axis_text_y_size,
                                Limit = Limit) {
    print(paste0(i, ": no data on positive strand"))
    df1_f <- data.frame(matrix(NA, nrow = 1, ncol = ncol(data_p)))
    colnames(df1_f) <- colnames(data_p)
    df1_f$ID <- "ID_fake"
    df1_f$position <- frag[i]
    df1_f$strand <- "+"
    df1_f$delay <- .01
    df1_f$half_life <- .01
    df1_f$intensity <- 1000
    Title <- NA

    #add arrow and virtual boxes to adjust the scale in ggplot
    p1 <-
      ggplot(df1_f, aes(x = get('position'), y = get('intensity'))) +
      scale_x_continuous(limits = c(frag[i], frag[i + 1])) +
      labs(y = "Intensity [A.U]") +
      scale_y_continuous(
        trans = 'log2',
        labels = label_log2_function,
        sec.axis = sec_axis(~ . * 1, name = "Coverage",
                            labels = label_square_function)
      ) +
      theme_bw() +
      background_grid(major = "xy", minor = "none") +
      theme(
        legend.text = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(colour = 5, size = axis_title_y_size),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(
          angle = 90,
          hjust = 1,
          size = axis_text_y_size
        ),
        panel.border = element_blank()
      )
    if (nrow(data_n) != 0) {
      p1 <- p1 +
        ggtitle(
          paste0(
            "ID: ",
            data_n$ID[1],
            "-",
            last(data_n$ID),
            "; ",
            "FC*: ",
            "significant t-test of two consecutive segments",
            "; Term: termination, NS: new start, PS: pausing site,",
            "iTSS_I: internal starting site,",
            "TI: transcription interferance."
          )
        ) +
        theme(plot.title = element_text(size = 6, hjust = .5))
    }
    p2 <-
      ggplot(df1_f, aes(x = get('position'), y = get('half_life'))) +
      scale_x_continuous(limits = c(frag[i], frag[i + 1])) +
      scale_y_continuous(
        limits = c(0, 10),
        breaks = seq(0, 10, by = 2),
        sec.axis = sec_axis(~ . * 1, name = "Half-life [min]", breaks =
                              seq(0, Limit, by = 2))
      ) +
      labs(y = "Half-life [min]") +
      theme_bw() +
      background_grid(major = "xy", minor = "none") +
      theme(
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(colour = 6, size = axis_title_y_size),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(
          angle = 90,
          hjust = 1,
          size = axis_text_y_size
        ),
        panel.border = element_blank()
      )
    p3 <-
      ggplot(df1_f, aes(x = get('position'), y = get('delay'))) +
      scale_x_continuous(limits = c(frag[i], frag[i + 1])) +
      scale_y_continuous(
        limits = c(0, Limit),
        breaks = seq(0, Limit, by = 2),
        sec.axis = sec_axis(~ . * 1, name = "Delay [min]", breaks =
                              seq(0, Limit, by = 2))
      ) +
      labs(y = "Delay [min]") +
      theme_bw() +
      background_grid(major = "xy", minor = "none") +
      theme(axis.title.x = element_blank()) +
      theme(
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(colour = 4, size = axis_title_y_size),
        axis.text.y = element_text(
          angle = 90,
          hjust = 1,
          size = axis_text_y_size
        ),
        axis.text.x = element_text(size = 6),
        panel.grid.major.x = element_blank(),
        plot.title = element_blank(),
        panel.border = element_blank()
      )
  p <- list(p1, p2, p3)
  return(p)
}