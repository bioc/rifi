##############################ggplot reverse strand#####################
negative_strand_function <- function(data_n,
                                     data,
                                     tmp.c2,
                                     frag,
                                     i,
                                     Limit = Limit,
                                     shape = shape,
                                     col_outiler = col_outiler,
                                     col_coverage = col_coverage,
                                     shape_outlier = shape_outlier,
                                     limit_intensity = limit_intensity,
                                     face = face,
                                     tick_length = tick_length,
                                     arrow.color = arrow.color,
                                     minVelocity = minVelocity,
                                     medianVelocity = medianVelocity,
                                     shape_above20 = shape_above20,
                                     col_above20 = col_above20,
                                     col_outlierabove10 = col_outlierabove10,
                                     shape_outlierabove10 = shape_outlierabove10,
                                     fontface = fontface,
                                     coverage = coverage,
                                     axis_text_y_size = axis_text_y_size,
                                     axis_title_y_size = axis_title_y_size,
                                     TI_threshold = TI_threshold,
                                     p_value_TI = p_value_TI,
                                     termination_threshold =
                                         termination_threshold,
                                     iTSS_threshold = iTSS_threshold,
                                     p_value_int = p_value_int,
                                     p_value_event = p_value_event,
                                     p_value_hl = p_value_hl,
                                     event_duration_ps = 
                                         event_duration_ps,
                                     event_duration_itss = 
                                         event_duration_itss,
                                     HL_threshold_1 = HL_threshold_1,
                                     HL_threshold_2 = HL_threshold_2,
                                     vel_threshold = vel_threshold,
                                     HL_threshold_color = 
                                         HL_threshold_color,
                                     vel_threshold_color = vel_threshold_color,
                                     ps_color = ps_color,
                                     iTSS_I_color = iTSS_I_color) {
    p4 <-
        ggplot(data, aes(x = get('position'),
                         y = get('intensity'))) +
        scale_x_continuous(limits = c(frag[i], frag[i + 1])) +
        scale_y_continuous(
            trans = 'log2',
            labels = label_log2_function,
            limits = c(NA, NA),
            sec.axis = sec_axis(~ . * 1, name = "Coverage",
                                labels = label_square_function)
        ) +
        coord_trans(y = 'reverse') +
        labs(y = "Intensity [A.U]") +
        theme_bw() +
        background_grid(major = "xy", minor = "none") +
        theme(
            legend.title = element_blank(),
            axis.title.y = element_text(colour = 5, size = axis_title_y_size),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(
                angle = 90,
                hjust = 1,
                size = axis_text_y_size
            ),
            panel.grid.major.x = element_blank(),
            legend.position = "none",
            plot.margin = margin(.1, .2, .1, .2, "cm"),
            panel.border = element_blank()
        )
    data_n <- indice_function(data_n, "HL_fragment")
    #increase the limit to 20 in case 3 or more probes/bins have a
    #delay above 10
    Limit_h_df2 <-
        limit_function(data_n, "half_life", ind = 1)
    if (Limit_h_df2 == 20) {
        Breaks_h2 <- seq(0, Limit_h_df2, by = 4)
    } else{
        Breaks_h2 <- seq(0, Limit_h_df2, by = 2)
    }
    df2.h <- secondaryAxis(data_n, "half_life", ind = 1)
    df2.h.o <- outlier_plot(data_n, "half_life", ind = 2, maxvalue = Limit_h_df2)
    #in case only one bin is available and the HL is above 20
    if (all(data_n$half_life > 20)) {
        data_n$half_life <- 20
    }
    p5 <-
        ggplot(data_n, aes(x = get('position'),
                           y = get('half_life'))) +
        scale_x_continuous(limits = c(frag[i], frag[c(i + 1)])) +
        scale_y_continuous(
            limits = c(0, Limit_h_df2),
            breaks = Breaks_h2,
            sec.axis = sec_axis(~ . * 1, name = "Half-life [min]", breaks =
                                    Breaks_h2)
        ) +
        labs(y = "Half-life [min]") +
        coord_trans(y = 'reverse') +
        theme_bw() +
        background_grid(major = "xy", minor = "none") +
        theme(
            legend.title = element_blank(),
            legend.position = "none",
            axis.title.y = element_text(colour = 6, size = axis_title_y_size),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(
                angle = 90,
                hjust = 1,
                size = axis_text_y_size
            ),
            panel.grid.major.x = element_blank(),
            plot.margin = margin(.1, .2, .1, .2, "cm"),
            panel.border = element_blank()
        )
    #in case of empty plot, a fake one is generated to avoid
    #disorganization.
    if (length(unique(data_n$half_life)) == 1) {
        if (is.na(unique(data_n$half_life))) {
            data_n$half_life[c(1, 2)] <- c(.1, .2)
            p5 <-
                ggplot(data_n, aes(
                    x = get('position'),
                    y = get('half_life')
                )) +
                scale_x_continuous(limits = c(frag[i], frag[i + 1])) +
                coord_trans(y = 'reverse') +
                labs(y = "Half-life [min]") +
                theme_bw() +
                background_grid(major = "xy", minor = "none") +
                theme(
                    legend.title = element_blank(),
                    legend.position = "none",
                    axis.title.y = element_text(colour = 6,
                                                size = axis_title_y_size),
                    axis.text.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    panel.grid.major.x = element_blank(),
                    axis.text.y = element_text(
                        angle = 90,
                        hjust = 1,
                        size = axis_text_y_size
                    ),
                    plot.margin = margin(.1, .2, .1, .2, "cm"),
                    panel.border = element_blank()
                )
        }
    }
    data_n <- indice_function(data_n, "delay_fragment")
    #increase the limit to 20 in case 3 or more probes/bins have a
    #delay above 10
    Limit_df2 <- limit_function(data_n, "delay", ind = 1)
    if (Limit_df2 == 20) {
        Breaks_d2 <- seq(0, Limit_df2, by = 4)
    } else{
        Breaks_d2 <- seq(0, Limit_df2, by = 2)
    }
    df2.d <- secondaryAxis(data_n, "delay", ind = 1)
    #in case only one bin is available and the delay is above 20
    if (all(data_n$delay > 20)) {
        data_n$delay <- 20
    }
    p6 <-
        ggplot(data_n, aes(x = get('position'), y = get('delay'))) +
        scale_x_continuous(limits = c(frag[i], frag[i + 1])) +
        scale_y_continuous(
            limits = c(0, Limit_df2),
            breaks = Breaks_d2,
            sec.axis = sec_axis(~ . * 1, name = "Delay [min]", breaks =
                                    Breaks_d2)
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
            axis.text.y = element_text(
                angle = 90,
                hjust = 1,
                size = axis_text_y_size
            ),
            panel.grid.major.x = element_blank(),
            plot.title = element_blank(),
            plot.margin = margin(.1, .2, .1, .2, "cm"),
            panel.border = element_blank()
        )
    #if delay data contains only one value or are all NAs, the plot is empty
    if (length(unique(data_n$delay)) == 1) {
        if (is.na(unique(data_n$delay))) {
            p6 <- p6 +
                scale_y_continuous(
                    limits = c(0, Limit_df2),
                    breaks = Breaks_d2,
                    sec.axis = sec_axis(~ . * 1, name = "Delay [min]", breaks =
                                            Breaks_d2)
                ) +
                coord_trans(y = 'reverse')
        }
    }
    #############################TI plot reverse strand###################
    #plot transcription interference on the background
    if (length(grep("TI", data_n$flag)) != 0) {
        df2_ti_a <- data_n[grep("TI", data_n$flag), ]
        #looping into TUs
        for (l in seq_along(unique(df2_ti_a$TU))) {
            df2_ti <- df2_ti_a[which(df2_ti_a$TU == unique(df2_ti_a$TU)[l]), ]
            df2_ti <-
                df2_ti[grep(paste0("\\TI_\\d+", "$"),
                            df2_ti$TI_termination_fragment), ]
            df2_ti <-
                df2_ti[!is.na(df2_ti$TI_mean_termination_factor), ]
            if (nrow(df2_ti) < 3 |
                nrow(df2_ti %>%
                     group_by(get('TI_termination_fragment'))) ==
                length(unique(df2_ti$TI_termination_fragment)) |
                all(df2_ti$TI_mean_termination_factor == 0)) {
                p4 <- p4
            } else{
                #plot the TI in parallel to the intensity lines. The values is
                #calculated as described below.
                df2_ti$value_d <-
                    df2_ti$intensity_mean_fragment /
                    (1 - df2_ti$TI_termination_factor)
                df2_ti$value_l <-
                    df2_ti$intensity_mean_fragment /
                    (1 - df2_ti$TI_mean_termination_factor)
                if (length(unique(df2_ti$TI_termination_fragment)) == 1) {
                    p4 <- p4 +
                        geom_line(data = df2_ti,
                                  aes(
                                      y = get('value_l'),
                                      group = get('TI_mean_termination_factor')
                                  ),
                                  size = .2) +
                        geom_point(
                            data = df2_ti,
                            aes(
                                y = get('value_d'),
                                group = get('TI_termination_fragment')
                            ),
                            size = .3,
                            shape = 1
                        )
                } else{
                    df2_ti <-
                        indice_function(df2_ti,
                                        "TI_termination_fragment")
                    df2_ti <- df2_ti %>%
                        filter(get('indice') == 1)
                    p4 <- p4 +
                        geom_step(data = df2_ti,
                                  aes(
                                      y = get('value_l'),
                                      group = get('TI_mean_termination_factor')
                                  ),
                                  size = .2) +
                        geom_point(
                            data = df2_ti,
                            aes(
                                y = get('value_d'),
                                group = get('TI_termination_fragment')
                            ),
                            size = .3,
                            shape = 1
                        )
                }
                #when strand is "-", the selection of the positions is opposite
                #as strand "+"
                df2.1_ti <-
                    arrange_byGroup(df2_ti, "TI_termination_fragment")
                if (nrow(df2.1_ti) > 1) {
                    TIs <- TI_frag_threshold(df2.1_ti, TI_threshold)
                    df2.1_ti_thr <-
                        df2.1_ti[which(df2.1_ti$TI_termination_fragment %in% TIs), ]
                    p4 <- p4 +
                        geom_vline(
                            df2.1_ti_thr,
                            xintercept = df2.1_ti_thr[, "position"],
                            colour = "green",
                            linetype = "dotted",
                            size = .2
                        )
                    if (length(na.omit(df2.1_ti$p_value_TI)) != 0) {
                        if (nrow(df2.1_ti %>%
                                 filter(get('p_value_TI')
                                        < p_value_TI)) != 0) {
                            p4 <- p4 +
                                geom_text(
                                    data = df2.1_ti,
                                    aes(
                                        x = get('position'),
                                        y = get('intensity_mean_fragment'),
                                        label = "Tinterf*"
                                    ),
                                    fontface = fontface,
                                    size = 1.3,
                                    check_overlap = TRUE,
                                    colour = "green"
                                )
                        }
                        else if (nrow(df2.1_ti %>%
                                      filter(get('p_value_TI')
                                             > p_value_TI)) != 0) {
                            p4 <- p4 +
                                geom_text(
                                    data = df2.1_ti,
                                    aes(
                                        x = get('position'),
                                        y = get('intensity_mean_fragment'),
                                        label = "Tinterf"
                                    ),
                                    fontface = fontface,
                                    size = 1.3,
                                    check_overlap = TRUE,
                                    colour = "green"
                                )
                        }
                    }
                }
            }
        }
    }
    #######################segments plot reverse strand###################
    if (length(na.omit(data_n$delay)) == 0) {
        p4 <- p4 +
            geom_point(size = .5)
        p5 <- p5
        p6 <- p6
    } else if (length(na.omit(data_n$delay)) < 3) {
        p4 <- p4 +
            geom_point(size = .5)
        p5 <- p5 +
            geom_point(size = .5)
        p6 <- p6 +
            geom_point(size = .5)
    } else{
        ######################intensity plot reverse strand################
        data_n <-
            indice_function(data_n, "intensity_fragment")
        if (nrow(data_n %>%
                 filter(get('indice') == 1)) != 0) {
            p4 <- p4 +
                geom_point(data = data_n %>%
                               filter(get('indice') == 1),
                           aes(col = get('intensity_fragment')),
                           size = .5)
            df2_wo <- data_n[which(data_n$indice == 1), ]
            df2_wo <-
                meanPosition(df2_wo, "intensity_fragment")
            if (length(df2_wo$intensity_fragment) != 0) {
                p4 <- p4 +
                    geom_line(data = data_n %>%
                                  filter(get('indice') == 1),
                              aes(
                                  x = get('position'),
                                  y = get('intensity_mean_fragment'),
                                  col = get('intensity_fragment')
                              )) +
                    geom_text(
                        data = df2_wo,
                        aes(
                            x = get('meanPosi'),
                            y = get('intensity_mean_fragment'),
                            label = get('intensity_fragment')
                        ),
                        size = 1.3,
                        check_overlap = TRUE
                    )
            }
        }
        if (nrow(data_n %>%
                 filter(get('indice') == 2)) != 0) {
            p4 <- p4 +
                geom_point(
                    data = data_n %>%
                        filter(get('indice') == 2),
                    col = col_outiler,
                    shape = shape_outlier,
                    size = .5
                )
        }
        if (coverage == 1) {
            p4 <- p4 +
                geom_line(data = tmp.c2,
                          aes(
                              x = get('position'),
                              y = get('coverage')
                          ),
                          col = get('col_coverage'))
        }
        ########################HL plot reverse strand#######################
        data_n <- indice_function(data_n, "HL_fragment")
        if (nrow(data_n %>%
                 filter(get('indice') == 1)) != 0) {
            p5 <- p5 +
                geom_point(data = data_n %>%
                               filter(get('indice') == 1),
                           aes(col = get('HL_fragment')),
                           size = .5)
            df2_wo <- data_n[which(data_n$indice == 1), ]
            df2_wo <-
                meanPosition(df2_wo, "HL_fragment")
            if (length(df2_wo$HL_fragment) != 0) {
                p5 <- p5 +
                    geom_line(data = data_n %>%
                                  filter(get('indice') == 1),
                              aes(
                                  x = get('position'),
                                  y = get('HL_mean_fragment'),
                                  col = get('HL_fragment')
                              )) +
                    geom_text(
                        data = df2_wo,
                        aes(
                            x = get('meanPosi'),
                            y = get('HL_mean_fragment'),
                            label = get('HL_fragment')
                        ),
                        size = 1.3,
                        check_overlap = TRUE
                    )
            }
        }
        if (nrow(data_n %>%
                 filter(get('indice') == 2)) != 0) {
            p5 <- p5 +
                geom_point(
                    data = data_n %>%
                        filter(get('indice') == 2),
                    col = col_outiler,
                    shape = shape_outlier,
                    size = .5
                )
        }
        if (nrow(df2.h %>%
                 filter(get('indice') == 1) %>%
                 filter(get('half_life') >= 20)) >= 1) {
            p5 <- p5 +
                geom_point(
                    data = df2.h %>%
                        filter(get('indice') == 1) %>%
                        filter(get('half_life') >= 20),
                    aes(y = Limit_h_df2),
                    col = col_above20,
                    shape = shape_above20,
                    size = .5
                )
        }
        #add outliers between 10 and 30
        if (nrow(df2.h.o) >= 1){
          p5 <- p5 +
            geom_point(
              data = df2.h.o,
          aes(y = Limit_h_df2),
          col = col_outlierabove10,
          shape = shape_outlierabove10,
          size = .5
          )
        }
        #######################delay plot reverse strand##################
        #plot delay segment and outliers
        data_n <- indice_function(data_n, "delay_fragment")
        if (nrow(data_n %>%
                 filter(get('indice') == 1)) != 0) {
            p6 <- p6 +
                geom_point(data = data_n %>%
                               filter(get('indice') == 1),
                           aes(col = get('delay_fragment')),
                           size = .5)
        }
        if (nrow(data_n %>%
                 filter(get('indice') == 2)) != 0) {
            p6 <- p6 +
                geom_point(
                    data = data_n %>%
                        filter(get('indice') == 2),
                    col = col_outiler,
                    shape = shape_outlier,
                    size = .5
                )
        }
        #plot the regression line from delay
        #plot of the delay or delay predicted is not always a line close to
        #the real values as delay data is very noisy,
        #the regression line could be far from the real data. In this cases,
        #an internal lm fit from ggplot is applied.
        if (length(na.omit(data_n$delay)) > 2) {
            df.c <- regr(data_n, ind = 1, data = data)
            #delay_mean column is added to draw a line in case of velocity
            #is NA or equal to 60.
            #the mean of the delay is calculated excluding outliers
            data_n$delay_mean <-
                delay_mean(data_n, "delay", ind = 1)
            #in case a fragment is split on two pages, each break get a
            #different velocity
            #upon the corresponding delay therefore a segment is checked
            #first if its flaking
            #the borders
            if (nrow(df.c) != 0) {
                p6 <- p6 +
                    geom_smooth(
                        data = df.c %>%
                            filter(get('indice') == 1),
                        formula = y ~ x,
                        aes(
                            y = get('delay.p_rev'),
                            col = get('delay_fragment')
                        ),
                        se = FALSE,
                        size = .4,
                        method = "lm"
                    )
            }
            if (nrow(data_n %>%
                     filter(get('indice') == 1) %>%
                     filter(get('slope') == 0)) > 3) {
                p6 <- p6 +
                    geom_line(
                        data = data_n %>%
                            filter(get('indice') == 1) %>%
                            filter(get('slope') == 0),
                        aes(
                            y = get('delay_mean'),
                            col = get('delay_fragment')
                        ),
                        size = .4
                    )
            }
        }
        if (nrow(df2.d %>%
                 filter(get('delay') >= 20)) >= 1) {
            p6 <- p6 +
                geom_point(
                    data = df2.d %>%
                        filter(get('delay') >= 20),
                    aes(y = Limit_df2),
                    col = col_above20,
                    shape = shape_above20,
                    size = .5
                )
        }
        df2_wo <-
            data_n[grep(paste0("\\D_\\d+", "$"), data_n$delay_fragment), ]
        if (nrow(df2_wo) != 0) {
            df2_wo <- meanPosition(df2_wo, "delay_fragment")
            df2_wo <-
                velocity_fun(df2_wo,
                             minVelocity = minVelocity,
                             medianVelocity = medianVelocity)
            if (length(unique(na.omit(df2_wo$velocity))) != 0) {
                p6 <- p6 +
                    geom_text(
                        data = df2_wo,
                        aes(
                            x = get('meanPosi'),
                            y = Limit - 1,
                            label = get('velocity')
                        ),
                        size = 1.3,
                        check_overlap = TRUE,
                        color = 2
                    ) +
                    geom_text(
                        data = df2_wo,
                        aes(
                            x = get('meanPosi'),
                            y = 5,
                            label = get('delay_fragment')
                        ),
                        size = 1.3,
                        check_overlap = TRUE
                    )
            }
        }
    }
    #######################Events plot reverse strand#####################
    #plot terminal and new start events
    forggtitle_ter <- 0
    forggtitle_NS <- 0
    df2_syR_T <- NA
    if (length(which(!is.na(data_n$synthesis_ratio))) != 0) {
        df2_syR <- data_n[which(!is.na(data_n$synthesis_ratio)), ]
        if (length(which(
            df2_syR$synthesis_ratio < termination_threshold &
            !is.na(df2_syR$FC_HL_intensity_fragment)
        )) != 0) {
            df2_syR_T <-
                df2_syR[which(
                    df2_syR$synthesis_ratio < termination_threshold &
                        !is.na(df2_syR$FC_HL_intensity_fragment)
                ), ]
            df2_syR_T <-
                df2_syR_T[!duplicated(df2_syR_T$FC_fragment_intensity), ]
            forggtitle_ter <- nrow(df2_syR_T)
        }
        df2_syR_T <- NA
        #plot New_start event from synthesis_ratio_event column
        if (length(which(
            df2_syR$synthesis_ratio > iTSS_threshold &
            !is.na(df2_syR$FC_HL_intensity_fragment)
        )) != 0) {
            df2_syR_T <-
                df2_syR[which(
                    df2_syR$synthesis_ratio > iTSS_threshold &
                        !is.na(df2_syR$FC_HL_intensity_fragment)
                ), ]
            df2_syR_T <-
                df2_syR_T[!duplicated(df2_syR_T$FC_fragment_intensity), ]
            forggtitle_NS <- nrow(df2_syR_T)
        }
    }
    ####################pausing site reverse strand##############
    #Select rows with pausing site event and draw a line
    if (nrow(data_n %>%
             filter(get('pausing_site') == "+")) != 0) {
        df2_ps <- data_n %>%
            filter(get('pausing_site') == "+")
        #select pausing site duration
        df2_ps <-
            df2_ps %>%
            filter(na.omit(get('event_duration')) >= event_duration_ps)
        #For those event with significant statistical test, are displayed
        #with a legend
        if (nrow(df2_ps %>%
                 filter(get('event_ps_itss_p_value_Ttest')
                        < p_value_event)) != 0) {
            df2_ps_s <-
                df2_ps %>%
                filter(get('event_ps_itss_p_value_Ttest')
                       < p_value_event)
            p6 <-
                my_segment_T(
                    p6,
                    data = df2_ps_s,
                    "PS*",
                    y = 0,
                    yend = 3,
                    dis = 50,
                    ytext = 3.8,
                    color = "orange",
                    linetype = "dotted",
                    df = "pausing",
                    fontface = fontface
                )
        }
        if (nrow(df2_ps %>%
                 filter(get('event_ps_itss_p_value_Ttest')
                        > p_value_event)) != 0) {
            df2_ps_b <-
                df2_ps %>%
                filter(get('event_ps_itss_p_value_Ttest') > p_value_event)
            p6 <-
                my_segment_T(
                    p6,
                    data = df2_ps_b,
                    "PS",
                    y = 0,
                    yend = 3,                    
                    dis = 50,
                    ytext = 3.8,
                    color = "orange",
                    linetype = "dotted",
                    df = "pausing",
                    fontface = fontface
                )
        }
        if (nrow(df2_ps %>%
                 filter(is.na(
                     get('event_ps_itss_p_value_Ttest')
                 ))) != 0) {
            df2_ps_b <- df2_ps[which(is.na(df2_ps$event_ps_itss_p_value_Ttest)), ]
            p6 <-
                my_segment_T(
                    p6,
                    data = df2_ps_b,
                    "PS",
                    y = 0,
                    yend = 3,                    
                    dis = 50,
                    ytext = 3.8,
                    color = "orange",
                    linetype = "dotted",
                    df = "pausing",
                    fontface = fontface
                )
        }
    }
    ####################iTSS_I reverse strand###################
    #Select rows with internal starting site event and draw a line
    if (nrow(data_n %>%
             filter(get('iTSS_I') == "+")) != 0) {
        #select iTSS
        df2_itss <- data_n %>%
            filter(get('iTSS_I') == "+")
        #select iTSS duration
        df2_itss <-
            df2_itss %>%
            filter(na.omit(get('event_duration')) <= event_duration_itss)
        if (nrow(df2_itss %>%
                 filter(get('event_ps_itss_p_value_Ttest')
                        < p_value_event)) != 0) {
            df2_itss_s <-
                df2_itss %>%
                filter(get('event_ps_itss_p_value_Ttest') < p_value_event)
            p6 <-
                my_segment_NS(
                    p6,
                    data = df2_itss_s,
                    "iTSS*",
                    y = 0,
                    yend = 3.2,
                    dis = 10,
                    ytext = 3.6,
                    color = 4,
                    linetype = "dotted",
                    fontface = fontface
                )
        }
        if (nrow(df2_itss %>%
                 filter(get('event_ps_itss_p_value_Ttest')
                        > p_value_event)) != 0) {
            df2_itss_b <-
                df2_itss %>%
                filter(get('event_ps_itss_p_value_Ttest') > p_value_event)
            p6 <-
                my_segment_NS(
                    p6,
                    data = df2_itss_b,
                    "iTSS",
                    y = 0,
                    yend = 3.2,
                    dis = 10,
                    ytext = 3.6,
                    color = 4,
                    linetype = "dotted",
                    fontface = fontface
                )
        }
        if (nrow(df2_itss %>%
                 filter(is.na(
                     get('event_ps_itss_p_value_Ttest')
                 ))) != 0) {
            df2_itss_b <- df2_itss[which(is.na(df2_itss$
                                                 event_ps_itss_p_value_Ttest)),]
            p6 <-
                my_segment_NS(
                    p6,
                    data = df2_itss_b,
                    "iTSS",
                    y = 0,
                    yend = 3.2,                    
                    dis = 10,
                    ytext = 3.6,
                    color = 4,
                    linetype = "dotted",
                    fontface = fontface
                )
        }
    }
    ####################FC reverse strand###################
    #add FC for intensity ratio test if p_value is significant
    data_n <- indice_function(data_n, "intensity_fragment")
    df2_wo <- data_n[which(data_n$indice == 1), ]
    #select the last row for each segment and add 40 nucleotides in case
    #of negative strand for a nice plot
    df2_wo_pvalue <-
        df2_wo[!duplicated(df2_wo$p_value_intensity), ]
    if (nrow(df2_wo_pvalue) != 0) {
        df2_p_val_int <-
            df2_wo_pvalue[which(df2_wo_pvalue$p_value_intensity
                                <= p_value_int), ]
        if (nrow(df2_p_val_int) != 0) {
            p4 <- p4 +
                geom_text(
                    data = df2_p_val_int,
                    aes(
                        x = get('position'),
                        y = get('intensity_mean_fragment')
                    ),
                    label = "FC*",
                    fontface = fontface,
                    size = 1,
                    check_overlap = TRUE
                )
        }
        df2_p_val_int <-
            df2_wo_pvalue[which(df2_wo_pvalue$p_value_intensity
                                > p_value_int), ]
        if (nrow(df2_p_val_int) != 0) {
            p4 <- p4 +
                geom_text(
                    data = df2_p_val_int,
                    aes(
                        x = get('position'),
                        y = get('intensity_mean_fragment')
                    ),
                    label = "FC",
                    fontface = fontface,
                    size = 1,
                    check_overlap = TRUE
                )
        }
    }
    #add FC for HL ratio higher than HL_threshold upon p_value significance
    data_n <- indice_function(data_n, "HL_fragment")
    df2_hl <- data_n[which(data_n$indice == 1), ]
    df2_hl <- df2_hl[!duplicated(df2_hl$FC_HL), ]
    df2_p_val_hl <-
        df2_hl[which(df2_hl$FC_HL >= HL_threshold_1), ]
    if (nrow(df2_p_val_hl) != 0) {
        df2_p_val_hl_p1 <-
            df2_p_val_hl[which(df2_p_val_hl$p_value_HL <= p_value_hl), ]
        if (nrow(df2_p_val_hl_p1) != 0) {
            p5 <- my_segment_NS(
                p5,
                data = df2_p_val_hl_p1,
                "HL*",
                y = 0,
                yend = 3,
                dis = 10,
                ytext = 3.4,
                color = HL_threshold_color,
                linetype = "dashed",
                fontface = fontface
            )
        }
        df2_p_val_hl_p2 <-
            df2_p_val_hl[which(df2_p_val_hl$p_value_HL > p_value_hl), ]
        if (nrow(df2_p_val_hl_p2) != 0) {
            p5 <- my_segment_NS(
                p5,
                data = df2_p_val_hl_p2,
                "HL",
                y = 0,
                yend = 3,
                dis = 10,
                ytext = 3.4,
                color = HL_threshold_color,
                linetype = "dashed",
                fontface = fontface
            )
        }
    }
    #add FC for HL ratio lower than HL_threshold upon p_value significance
    df2_p_val_hl <-
        df2_hl[which(df2_hl$FC_HL <= HL_threshold_2), ]
    if (nrow(df2_p_val_hl) != 0) {
        df2_p_val_hl_p1 <-
            df2_p_val_hl[which(df2_p_val_hl$p_value_HL <= p_value_hl), ]
        if (nrow(df2_p_val_hl_p1) != 0) {
            p5 <- 
                my_segment_NS(
                    p5,
                    data = df2_p_val_hl_p1,
                    "HL*",
                    y = 0,
                    yend = 3,
                    dis = 10,
                    ytext = 3.4,
                    color = HL_threshold_color,
                    linetype = "dashed",
                    fontface = fontface
                )
        }
        df2_p_val_hl_p2 <-
            df2_p_val_hl[which(df2_p_val_hl$p_value_HL > p_value_hl), ]
        if (nrow(df2_p_val_hl_p2) != 0) {
            p5 <- 
                my_segment_NS(
                    p5,
                    data = df2_p_val_hl_p2,
                    "HL",
                    y = 0,
                    yend = 3,
                    dis = 10,
                    ytext = 3.4,
                    color = HL_threshold_color,
                    linetype = "dashed",
                    fontface = fontface
                )
        }
    }
    #if FC_HL is between both threshold and has a significant p_value
    df2_p_val_hl <-
        df2_hl[which(df2_hl$FC_HL > HL_threshold_2 & 
                         df2_hl$FC_HL < HL_threshold_1),]
    if (nrow(df2_p_val_hl) != 0) {
        df2_p_val_hl_p1 <-
            df2_p_val_hl[which(df2_p_val_hl$p_value_HL <= p_value_hl), ]
        if (nrow(df2_p_val_hl_p1) != 0) {
            p5 <- 
                my_segment_NS(
                    p5,
                    data = df2_p_val_hl_p1,
                    "HL*",
                    y = 0,
                    yend = 3,
                    dis = 10,
                    ytext = 3.4,
                    color = HL_threshold_color,
                    linetype = "dashed",
                    fontface = fontface
                )
        }
    }
    #Select rows with velocity_ratio event and draw a line
    df2_v <- data_n[!duplicated(data_n$velocity_ratio), ]
    df2_v <-
        df2_v[which(df2_v$velocity_ratio < vel_threshold), ]
    if (nrow(df2_v) != 0) {
        p6 <-
            my_segment_NS(
                p6,
                data = df2_v,
                "V",
                y = 0,
                yend = 3,
                dis = 10,
                ytext = 3.4,
                color = vel_threshold_color,
                linetype = "dashed",
                fontface = fontface
            )
    }
    ############################Title and the plot#####################
    if (nrow(data_n) != 0) {
        Title <-
            paste0(
                "Term: termination (",
                forggtitle_ter,
                "), NS: new start (",
                forggtitle_NS,
                "), PS: pausing site (",
                length(which(data_n$pausing_site == "+")),
                "), iTSS_I: internal starting site (",
                length(which(data_n$iTSS_I == "+")),
                ")"
            )
    }
    p <- list(p6, p5, p4, Title)
    return(p)
}
