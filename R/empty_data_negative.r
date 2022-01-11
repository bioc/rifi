#in case the dataframe from the positive strand is empty, a dataframe is
#created to be plotted for the genome annotation.
empty_data_negative <- function(data_n,
                                frag,
                                i,
                                axis_title_y_size = axis_title_y_size,
                                axis_text_y_size = axis_text_y_size,
                                Limit = Limit) {
  print(paste0(i, ": no data on negative strand"))
  df2_f <- data.frame(matrix(NA, nrow = 1, ncol = ncol(data_n)))
  colnames(df2_f) <- colnames(data_n)
  df2_f$ID <- "ID_fake"
  df2_f$position <- frag[i]
  df2_f$strand <- "-"
  df2_f$delay <- .01
  df2_f$half_life <- .01
  df2_f$intensity <- 1000
  Title <- NA

  #generate empty data
  p4 <-
    ggplot(df2_f, aes(x = get('position'), y = get('intensity'))) +
    scale_x_continuous(limits = c(frag[i], frag[i + 1])) +
    scale_y_continuous(
      trans = 'log2',
      labels = label_log2_function,
      limits = c(NA, NA),
      sec.axis = sec_axis( ~ . * 1, name = "Coverage [A.U]")
    ) +
    labs(y = "Intensity [A.U]") +
    theme_bw() +
    background_grid(major = "xy", minor = "none") +
    theme(
      legend.title = element_blank(),
      legend.position = "none",
      axis.title.y = element_text(colour = 5, size = axis_title_y_size),
      axis.text.x = element_blank(),
      axis.text.y = element_text(
        angle = 90,
        hjust = 1,
        size = axis_text_y_size
      ),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.margin = margin(.1, .2, .1, .2, "cm"),
      panel.border = element_blank()
    )
  p5 <-
    ggplot(df2_f, aes(x = get('position'), y = get('half_life'))) +
    scale_x_continuous(limits = c(frag[i], frag[i + 1])) +
    scale_y_continuous(
      limits = c(0, 10),
      breaks = seq(0, 10, by = 2),
      sec.axis = sec_axis( ~ . * 1, name = "Half-life [min]", breaks =
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
      axis.text.y = element_text(
        angle = 90,
        hjust = 1,
        size = axis_text_y_size
      ),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.margin = margin(.1, .2, .1, .2, "cm"),
      panel.border = element_blank()
    )
  p6  <-
    ggplot(df2_f, aes(x = get('position'), y = get('delay'))) +
    scale_x_continuous(limits = c(frag[i], frag[i + 1])) +
    scale_y_continuous(
      limits = c(0, 10),
      breaks = seq(0, 10, by = 2),
      sec.axis = sec_axis( ~ . * 1, name = "Delay [min]", breaks =
                             seq(0, Limit, by = 2))
    ) +
    labs(y = "Delay [min]") +
    coord_trans(y = 'reverse') +
    theme_bw() +
    background_grid(major = "xy", minor = "none") +
    theme(
      legend.title = element_blank(),
      legend.position = "none",
      axis.title.y = element_text(colour = 4, size = axis_title_y_size),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.y = element_text(
        angle = 90,
        hjust = 1,
        size = axis_text_y_size
      ),
      plot.margin = margin(.1, .2, .1, .2, "cm"),
      plot.title = element_blank(),
      panel.border = element_blank()
    )
    Title <- NA
  p <- list(p6, p5, p4, Title)
  return(p)
}
