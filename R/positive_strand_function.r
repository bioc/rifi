##########################positive strand############################
#This part of the script runs only in case of the dataframe of
#positive strand is not empty
########################ggplot positive strand#######################
positive_strand_function <- function(data_p, data, tmp.c1, df1_1, frag, i,
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
                                     HL_threshold_color_2 = 
                                         HL_threshold_color_2,
                                     HL_threshold_color_1 = 
                                         HL_threshold_color_1,
                                     vel_threshold_color = vel_threshold_color,
                                     ps_color = ps_color,
                                     iTSS_I_color = iTSS_I_color) {
       #first plot for intensity segments
    p1 <-
        ggplot(data_p, aes(
            x = get('position'),
            y = get('intensity')
        )) +
        scale_x_continuous(limits = c(frag[i], frag[i + 1])) +
        scale_y_continuous(
            trans = 'log2',
            labels = label_log2_function,
            limits = c(NA, NA),
            sec.axis = sec_axis( ~ . * 1, name = "Coverage",
                                 labels = label_square_function)
        ) +
        labs(y = "Intensity [A.U]") +
        theme_bw() +
        background_grid(major = "xy", minor = "none") +
        theme(
            legend.title = element_blank(),
            axis.title.y = element_text(colour = 5, size = axis_title_y_size),
            axis.text.y = element_text(
                angle = 90,
                hjust = 1,
                size = axis_text_y_size
            ),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            legend.position = "none",
            plot.margin = margin(.1, .2, .1, .2, "cm"),
            panel.border = element_blank()
        )
    ###########################TI plot positive strand####################
    #plot Transcription interference, upon the p_value, if significant, an
    #asterisk is added.
    if (length(grep("TI", data_p$flag)) != 0) {
        #grep only flagged probes/bins with TI
        df1_ti_a <- data_p[grep("TI", data_p$flag),]
        for (l in seq_along(unique(df1_ti_a$TU))) {
            df1_ti <- df1_ti_a[which(df1_ti_a$TU == unique(df1_ti_a$TU)[l]),]
            #grep only TI without any termination as NA, T or O
            df1_ti <-
                df1_ti[grep(paste0("\\TI_\\d+", "$"),
                            df1_ti$TI_termination_fragment),]
            #grep only non-NA rows
            df1_ti <-
                df1_ti[!is.na(df1_ti$TI_mean_termination_factor),]
            # Eliminate TI factor above 1
            df1_ti <-
                df1_ti[which(df1_ti$TI_termination_factor < 1),]
            #TI is ignored if it contains only 3 probes/bins
            if (nrow(df1_ti) < 3 |
                nrow(df1_ti %>% group_by(get(
                    'TI_termination_fragment'
                ))) ==
                length(unique(df1_ti$TI_termination_fragment)) |
                all(df1_ti$TI_mean_termination_factor == 0)) {
                p1 <- p1
                #geom_step runs only with more observations for each group
                #exists, in opposite case only geom_line from intensity is
                #plotted.
            } else{
                # add a value corresponding to the TI factor relative to the
                # intensity values the value plotted dividing the intensity
                # mean by (1 - TI factor)
                df1_ti$value_d <-
                    df1_ti$intensity_mean_fragment /
                    (1 - df1_ti$TI_termination_factor)
                df1_ti$value_l <-
                    df1_ti$intensity_mean_fragment /
                    (1 - df1_ti$TI_mean_termination_factor)
                if (length(unique(df1_ti$TI_termination_fragment)) == 1) {
                    p1 <- p1 +
                        geom_line(
                            data = df1_ti,
                            aes(
                                y = get('value_l'),
                                group = get('TI_mean_termination_factor')
                            ),
                            size = .2
                        ) +
                        geom_point(
                            data = df1_ti,
                            aes(
                                y = get('value_d'),
                                group = get('TI_termination_fragment')
                            ),
                            size = .3,
                            shape = 1
                        )
                } else{
                    #plot the corresponding TI fragment, only an horizontal line
                    # is added if two TIs segments are present
                    df1_ti <-
                        indice_function(df1_ti,
                                        "TI_termination_fragment")
                    df1_ti <-
                        df1_ti %>% filter(get('indice') == 1)
                    p1 <- p1 +
                        geom_step(
                            data = df1_ti,
                            aes(
                                y = get('value_l'),
                                group = get('TI_mean_termination_factor')
                            ),
                            size = .2
                        ) +
                        geom_point(
                            data = df1_ti,
                            aes(
                                y = get('value_d'),
                                group = get('TI_termination_fragment')
                            ),
                            size = .3,
                            shape = 1
                        )
                }
                # To draw green lines, we need the first dataframe with all TI
                # fragments. TI factors equal to 0
                # at the beginning of the TI represent the start of the TI and
                # is plotted.
                # we select only the first position of each TI fragment
                df1.1_ti <-
                    df1_ti[!duplicated(df1_ti$TI_termination_fragment),]
                #plot vertical lines on the middle of each 2 fragments.
                if (nrow(df1.1_ti) > 1) {
                    TIs <- TI_frag_threshold(df1.1_ti, TI_threshold)
                    df1.1_ti_thr <-
                        df1.1_ti[which(
                          df1.1_ti$TI_termination_fragment %in% TIs),]
                    p1 <- p1 +
                        geom_vline(
                            df1.1_ti_thr,
                            xintercept = df1.1_ti_thr[, "position"],
                            colour = "green",
                            linetype = "dotted",
                            size = .2
                        )
                    # the condition below is in case TI t-test is NA for some
                    # reasons
                    if (length(na.omit(df1.1_ti$p_value_TI)) != 0) {
                        if (nrow(df1.1_ti %>%
                                 filter(
                                     get('p_value_TI') < p_value_TI
                                 )) != 0) {
                            p1 <- p1 +
                                geom_text(
                                    data = df1.1_ti,
                                    aes(
                                        x = get('position'),
                                        y = get(
                                            'intensity_mean_fragment'
                                        ),
                                        label = "Tinterf*"
                                    ),
                                    fontface = fontface,
                                    size = 1.3,
                                    check_overlap = TRUE,
                                    colour = "green"
                                )
                        } else if (nrow(df1.1_ti %>% filter(
                            get('p_value_TI') >
                            p_value_TI
                        )) != 0) {
                            p1 <- p1 +
                                geom_text(
                                    data = df1.1_ti,
                                    aes(
                                        x = get('position'),
                                        y = get(
                                            'intensity_mean_fragment'
                                        ),
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
    #######################segment plot positive strand###################
    #first plot for half-life segments
    #select segments without outliers
    data_p <- indice_function(data_p, "HL_fragment")
    # increase the limit to 20 in case 3 or more probes/bins have a HL
    # above 10
    Limit_h_df1 <-
        limit_function(data_p, "half_life", ind = 1)
    if (Limit_h_df1 == 20) {
        Breaks_h <- seq(0, Limit_h_df1, by = 4)
    } else{
        Breaks_h <- seq(0, Limit_h_df1, by = 2)
    }
    df1.h <- secondaryAxis(data_p, "half_life", ind = 1)
    #in case only one bin is available and the HL is above 20
    if (all(data_p$half_life > 20)) {
        data_p$half_life <- 20
    }
    p2 <-
        ggplot(data_p, aes(
            x = get('position'),
            y = get('half_life')
        )) +
        scale_x_continuous(limits = c(frag[i], frag[i + 1])) +
        scale_y_continuous(
            limits = c(0, Limit_h_df1),
            breaks = Breaks_h,
            sec.axis = sec_axis( ~ . * 1, name = "Half-life [min]", breaks =
                                     Breaks_h)
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
            axis.text.y = element_text(
                angle = 90,
                hjust = 1,
                size = axis_text_y_size
            ),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            plot.margin = margin(.1, .2, .1, .2, "cm"),
            panel.border = element_blank()
        )
    #add the second axis for half-life segments plot
    if (length(unique(data_p$half_life)) == 1) {
        if (is.na(unique(data_p$half_life))) {
            p2 <- p2 +
                geom_blank() +
                scale_y_continuous(
                    limits = c(0, Limit_h_df1),
                    breaks = Breaks_h,
                    sec.axis = sec_axis( ~ . * 1, name = "Half-life [min]",
                                         breaks = Breaks_h)
                )
        }
    }
    #select segments without outliers
    data_p <- indice_function(data_p, "delay_fragment")
    #increase the limit to 20 in case 3 or more probes/bins have a delay
    #above 10
    Limit_df1 <- limit_function(data_p, "delay", ind = 1)
    if (Limit_df1 == 20) {
        Breaks_d <- seq(0, Limit_df1, by = 4)
    } else{
        Breaks_d <- seq(0, Limit_df1, by = 2)
    }
    #first plot for delay segments
    df1.d <- secondaryAxis(data_p, "delay", ind = 1)
    #in case only one bin is available and the delay is above 20
    if (all(data_p$delay > 20)) {
        data_p$delay <- 20
    }
    p3 <-
        ggplot(data_p, aes(x = get('position'), y = get('delay'))) +
        scale_x_continuous(limits = c(frag[i], frag[i + 1])) +
        scale_y_continuous(
            limits = c(0, Limit_df1),
            breaks = Breaks_d,
            sec.axis = sec_axis( ~ . * 1, name = "Delay [min]", breaks =
                                     Breaks_d)
        ) +
        labs(y = "Delay [min]") +
        theme_bw() +
        background_grid(major = "xy", minor = "none") +
        theme(
            legend.title = element_blank(),
            axis.title.x = element_blank(),
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
            plot.margin = margin(.1, .2, .2, .2, "cm"),
            panel.border = element_blank()
        )
    #add the second axis for delay segments plot
    if (length(unique(data_p$delay)) == 1) {
        if (is.na(unique(data_p$delay))) {
            p3 <- p3 +
                scale_y_continuous(
                    limits = c(0, Limit_df1),
                    breaks = Breaks_d,
                    sec.axis = sec_axis( ~ . * 1, name = "Delay [min]",
                                         breaks = Breaks_d)
                )
        }
    }
    #add the segments to the plot base and check
    #intensity is plotted independently if delay/HL data are present or not.
    if (length(na.omit(data_p$delay)) == 0) {
        p1 <- p1 +
            geom_point(size = .5)
        p2 <- p2
        p3 <- p3
    } else if (length(na.omit(data_p$delay)) < 3) {
        p1 <- p1 +
            geom_point(size = .5)
        p2 <- p2 +
            geom_point(size = .5)
        p3 <- p3 +
            geom_point(size = .5)
    } else{
        ######################intensity plot positive strand################
        #add a reference to outliers probes or bins
        data_p <-
            indice_function(data_p, "intensity_fragment")
        #intensities/HL and delay with NA are not plotted, high intensities
        #are plotted with a different shape,
        #outliers probes/bins are plotted with a different color, each
        #segment has a color, a line is added for each segment
        #indicating the mean. "FC*" indicate significant t-test and each
        #segment is labeled.
        if (nrow(data_p %>% filter(get('indice') == 1)) != 0) {
            p1 <- p1 +
                geom_point(
                    data = data_p %>% filter(get('indice') == 1),
                    aes(col = get(
                        'intensity_fragment'
                    )),
                    size = .5
                )
            #eliminate outliers probes or bins
            df1_wo <- data_p[which(data_p$indice == 1),]
            #assign a mean position for the plot
            df1_wo <-
                meanPosition(df1_wo, "intensity_fragment")
            if (length(df1_wo$intensity_fragment) != 0) {
                p1 <- p1 +
                    geom_line(data = data_p %>% filter(get('indice') == 1),
                              aes(
                                  x = get('position'),
                                  y = get('intensity_mean_fragment'),
                                  col = get('intensity_fragment')
                              )) +
                    geom_text(
                        data = df1_wo,
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
        #plot intensity outliers
        if (nrow(data_p %>% filter(get('indice') == 2)) != 0) {
            p1 <- p1 +
                geom_point(
                    data = data_p %>% filter(get('indice') == 2),
                    col = col_outiler,
                    shape = shape_outlier,
                    size = .5
                )
        }
        #in case of coverage is available.
        if (coverage == 1) {
            p1 <- p1 +
                geom_line(data = tmp.c1,
                          aes(x = get('position'), y = coverage),
                          col = col_coverage)
        }
        #######################HL plot positive strand###################
        #plot bins and mean fragments
        data_p <- indice_function(data_p, "HL_fragment")
        if (nrow(data_p %>% filter(get('indice') == 1)) != 0) {
            p2 <- p2 +
                geom_point(
                    data = data_p %>% filter(get('indice') == 1),
                    aes(col = get('HL_fragment')),
                    size = .5
                )
            #eliminate outliers probes or bins
            df1_wo <- data_p[which(data_p$indice == 1),]
            #assign a mean position for the plot
            df1_wo <-
                meanPosition(df1_wo, "HL_fragment")
            if (length(df1_wo$HL_fragment) != 0) {
                p2 <- p2 +
                    geom_line(data = data_p %>% filter(get('indice') == 1),
                              aes(
                                  x = get('position'),
                                  y = get('HL_mean_fragment'),
                                  col = get('HL_fragment')
                              )) +
                    geom_text(
                        data = df1_wo,
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
        #plot HL outliers
        if (nrow(data_p %>% filter(get('indice') == 2)) != 0) {
            p2 <- p2 +
                geom_point(
                    data = data_p %>% filter(get('indice') == 2),
                    col = col_outiler,
                    shape = shape_outlier,
                    size = .5
                )
        }
        #add replaced values to 20
        if (nrow(df1.h %>%
                 filter(get('indice') == 1) %>%
                 filter(get('half_life') >= 20)) >= 1) {
            p2 <- p2 +
                geom_point(
                    data = df1.h %>%
                        filter(get('indice') == 1) %>%
                        filter(get('half_life') >= 20),
                    aes(y = Limit_h_df1),
                    col = col_above20,
                    shape = shape_above20,
                    size = .5
                )
        }
        ########################delay plot positive strand#################
        #assign the indice to df1
        data_p <- indice_function(data_p, "delay_fragment")
        #plot delay bins
        if (nrow(data_p %>%
                 filter(get('indice') == 1)) != 0) {
            p3 <- p3 +
                geom_point(
                    data = data_p %>%
                        filter(get('indice') == 1),
                    aes(col = get('delay_fragment')),
                    size = .5
                )
        }
        #plot delay outliers
        if (nrow(data_p %>%
                 filter(get('indice') == 2)) != 0) {
            p3 <- p3 +
                geom_point(
                    data = data_p %>%
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
        #in case of slope of delay fragments is 0, the intercept is plotted
        #otherwise the delay (delay.p) is predicted from the slope and
        #plotted.
        if (length(na.omit(data_p$delay)) > 2) {
            df.c <- regr(data_p, ind = 1, data = data)
            #delay_mean column is added to draw a line in case of velocity
            #is NA or equal to 60.
            #the mean of the delay is calculated excluding outliers
            data_p$delay_mean <-
                delay_mean(data_p, "delay", ind = 1)
            if (nrow(df.c) != 0) {
                p3 <- p3 +
                    geom_line(
                        data = df.c %>%
                            filter(get('indice') == 1),
                        aes(
                            y = get('predicted_delay'),
                            col = get('delay_fragment')
                        ),
                        size = .4
                    )
            }
            if (nrow(data_p %>%
                     filter(get('indice') == 1) %>%
                     filter(get('slope') == 0)) > 3) {
                p3 <- p3 +
                    geom_line(
                        data = data_p %>%
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
        #plot those bins over 20 minutes
        if (nrow(df1.d %>%
                 filter(get('delay') >= 20)) >= 1) {
            p3 <- p3 +
                geom_point(
                    data = df1.d %>%
                        filter(get('delay') >= 20),
                    aes(y = Limit_df1),
                    col = col_above20,
                    shape = shape_above20,
                    size = .5
                )
        }
        df1_wo <-
            data_p[grep(paste0("\\D_\\d+", "$"), data_p$delay_fragment),]
        if (nrow(df1_wo) != 0) {
            df1_wo <- meanPosition(df1_wo, "delay_fragment")
            df1_wo <-
                velocity_fun(df1_wo,
                             minVelocity = minVelocity,
                             medianVelocity = medianVelocity)
            if (length(unique(na.omit(df1_wo$velocity))) != 0) {
                p3 <- p3 +
                    geom_text(
                        data = df1_wo,
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
                        data = df1_wo,
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
    #########################Events plot positive strand#################
    #plot termination event from synthesis_ratio_event column
    forggtitle_ter <- 0
    forggtitle_NS <- 0
    df1_syR_T <- NA
    if (length(which(!is.na(data_p$synthesis_ratio))) != 0) {
        df1_syR <- data_p[which(!is.na(data_p$synthesis_ratio)),]
        #in case last position matches with an event which needs to be
        #on the next page of the plot.
        if (last(df1_syR$position) == frag[c(i + 1)]) {
            fc_seg <-
                df1_syR[which(df1_syR$position ==
                                  frag[c(i + 1)]), "FC_HL_intensity_fragment"]
            df1_syR[which(df1_syR$FC_HL_intensity_fragment == fc_seg),
                    c("synthesis_ratio_event",
                      "p_value_Manova")] <- NA
        }
        if (length(which(
            df1_syR$synthesis_ratio < termination_threshold &
            !is.na(df1_syR$FC_HL_intensity_fragment)
        )) != 0) {
            df1_syR_T <-
                df1_syR[which(
                    df1_syR$synthesis_ratio < termination_threshold &
                        !is.na(df1_syR$FC_HL_intensity_fragment)
                ),]
            df1_syR_T <-
                arrange_byGroup(df1_syR_T, "FC_HL_intensity_fragment")
            forggtitle_ter <- nrow(df1_syR_T)
        }
        df1_syR_T <- NA
        #plot New_start event from synthesis_ratio_event
        if (length(which(
            df1_syR$synthesis_ratio > iTSS_threshold &
            !is.na(df1_syR$FC_HL_intensity_fragment)
        )) != 0) {
            df1_syR_T <-
                df1_syR[which(
                    df1_syR$synthesis_ratio > iTSS_threshold &
                        !is.na(df1_syR$FC_HL_intensity_fragment)
                ),]
            df1_syR_T <-
                arrange_byGroup(df1_syR_T, "FC_HL_intensity_fragment")
            forggtitle_NS <- nrow(df1_syR_T)
        }
    }
    ####################pausing site positive strand##############
    #plot pausing sites events
    if (nrow(data_p %>%
             filter(get('pausing_site') == "+")) != 0) {
        #select pausing site
        df1_ps <- data_p %>%
            filter(get('pausing_site') == "+")
        #select pausing site duration
        df1_ps <-
            df1_ps %>%
            filter(na.omit(get('event_duration')) >= event_duration_ps)
        if (nrow(df1_ps %>%
                 filter(
                     get('event_ps_itss_p_value_Ttest')
                     < p_value_event
                 )) != 0) {
            df1_ps_s <-
                df1_ps %>%
                filter(get('event_ps_itss_p_value_Ttest')
                       < p_value_event)
            p3 <-
                my_segment_T(
                    p3,
                    data = df1_ps_s,
                    "PS*",
                    y = 0,
                    yend = 3,
                    dis = 50,
                    ytext = 3.8,
                    color = "orange",
                    linetype = "dashed",
                    df = "pausing",
                    fontface = fontface
                )
        }
        if (nrow(df1_ps %>%
                 filter(
                     get('event_ps_itss_p_value_Ttest')
                     > p_value_event
                 )) != 0) {
            df1_ps_b <-
                df1_ps %>%
                filter(get('event_ps_itss_p_value_Ttest')
                       > p_value_event)
            p3 <-
                my_segment_T(
                    p3,
                    data = df1_ps_b,
                    "PS",
                    y = 0,
                    yend = 3,
                    dis = 50,
                    ytext = 3.8,
                    color = "orange",
                    linetype = "dashed",
                    df = "pausing",
                    fontface = fontface
                )
        }
        if (nrow(df1_ps %>%
                 filter(is.na(
                     get('event_ps_itss_p_value_Ttest')
                 ))) != 0) {
            df1_ps_s <- df1_ps %>%
                filter(is.na('event_ps_itss_p_value_Ttest'))
            p3 <-
                my_segment_T(
                    p3,
                    data = df1_ps_s,
                    "PS",
                    y = 0,
                    yend = 3,
                    dis = 50,
                    ytext = 3.8,
                    color = "orange",
                    linetype = "dashed",
                    df = "pausing",
                    fontface = fontface
                )
        }
    }
    ####################iTSS_I positive strand###################
    #plot internal starting sites events
    if (nrow(data_p %>%
             filter(get('iTSS_I') == "+")) != 0) {
        #select iTSS
        df1_itss <- data_p %>%
            filter(get('iTSS_I') == "+")
        #select iTSS duration
        df1_itss <-
            df1_itss %>%
            filter(na.omit(get('event_duration'))
                   <= event_duration_itss)
        if (nrow(df1_itss %>%
                 filter(
                     get('event_ps_itss_p_value_Ttest')
                     < p_value_event
                 )) != 0) {
            df1_itss_s <-
                df1_itss %>%
                filter(get('event_ps_itss_p_value_Ttest')
                       < p_value_event)
            p3 <-
                my_segment_NS(
                    p3,
                    data = df1_itss_s,
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
        if (nrow(df1_itss %>%
                 filter(
                     get('event_ps_itss_p_value_Ttest')
                     > p_value_event
                 )) != 0) {
            df1_itss_b <-
                df1_itss %>%
                filter(get('event_ps_itss_p_value_Ttest')
                       > p_value_event)
            p3 <-
                my_segment_NS(
                    p3,
                    data = df1_itss_b,
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
        if (nrow(df1_itss %>%
                 filter(is.na(
                     'event_ps_itss_p_value_Ttest'
                 ))) != 0) {
            df1_itss_b <-
                df1_itss %>%
                filter(is.na('event_ps_itss_p_value_Ttest'))
            p3 <-
                my_segment_NS(
                    p3,
                    data = df1_itss_b,
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
    ####################FC positive strand###################
    #add FC for intensity ratio test if p_value is significant
    #select the last row for each segment and add 40 nucleotides in case
    #of negative strand for a nice plot
    data_p <- indice_function(data_p, "intensity_fragment")
    df1_wo <- data_p[which(data_p$indice == 1),]
    #select the last row for each segment and add 40 nucleotides in case
    #of negative strand for a nice plot
    df1_wo_pvalue <-
        arrange_byGroup(df1_wo, "p_value_intensity")
    if (nrow(df1_wo_pvalue) != 0) {
        df1_p_val_int <-
            df1_wo_pvalue[which(df1_wo_pvalue$p_value_intensity
                                < p_value_int),]
        if (nrow(df1_p_val_int) != 0) {
            p1 <- p1 +
                geom_text(
                    data = df1_p_val_int,
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
        df1_p_val_int <-
            df1_wo_pvalue[which(df1_wo_pvalue$p_value_intensity
                                > p_value_int),]
        if (nrow(df1_p_val_int) != 0) {
            p1 <- p1 +
                geom_text(
                    data = df1_p_val_int,
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
    #add FC for HL ratio lower than HL_threshold upon p_value significance
    data_p <- indice_function(data_p, "HL_fragment")
    df1_wo <- data_p[which(data_p$indice == 1),]
    df1_hl <- arrange_byGroup(data_p, "FC_HL")
    df1_p_val_hl <-
        df1_hl[which(df1_hl$FC_HL >= HL_threshold_1),]
    if (nrow(df1_p_val_hl) != 0) {
        df1_p_val_hl_p1 <-
            df1_p_val_hl[which(df1_p_val_hl$p_value_HL <= p_value_hl),]
        if (nrow(df1_p_val_hl_p1) != 0) {
            p2 <- my_segment_NS(
                p2,
                data = df1_p_val_hl_p1,
                "HL*",
                y = 0,
                yend = 3,
                dis = 10,
                ytext = 3.4,
                color = HL_threshold_color_1,
                linetype = "dashed",
                fontface = fontface
            )
        }
        df1_p_val_hl_p2 <-
            df1_p_val_hl[which(df1_p_val_hl$p_value_HL > p_value_hl),]
        if (nrow(df1_p_val_hl_p2) != 0) {
            p2 <- my_segment_NS(
                p2,
                data = df1_p_val_hl_p2,
                "HL",
                y = 0,
                yend = 3,
                dis = 10,
                ytext = 3.4,
                color = HL_threshold_color_1,
                linetype = "dashed",
                fontface = fontface
            )
        }
    }
    #add FC for HL ratio higher than HL_threshold upon p_value significance
    df1_p_val_hl <-
        df1_hl[which(df1_hl$FC_HL <= HL_threshold_2),]
    if (nrow(df1_p_val_hl) != 0) {
        df1_p_val_hl_p1 <-
            df1_p_val_hl[which(df1_p_val_hl$p_value_HL <= p_value_hl),]
        if (nrow(df1_p_val_hl_p1) != 0) {
            p2 <- my_segment_NS(
                p2,
                data = df1_p_val_hl_p1,
                "HL*",
                y = 0,
                yend = 3,
                dis = 10,
                ytext = 3.4,
                color = HL_threshold_color_1,
                linetype = "dashed",
                fontface = fontface
            )
        }
        df1_p_val_hl_p2 <-
            df1_p_val_hl[which(df1_p_val_hl$p_value_HL > p_value_hl),]
        if (nrow(df1_p_val_hl_p2) != 0) {
            p2 <- my_segment_NS(
                p2,
                data = df1_p_val_hl_p2,
                "HL",
                y = 0,
                yend = 3,
                dis = 10,
                ytext = 3.4,
                color = HL_threshold_color_1,
                linetype = "dashed",
                fontface = fontface
            )
        }
    }
    df1_p_val_hl <-
        df1_hl[which(df1_hl$FC_HL > HL_threshold_2 & 
                         df1_hl$FC_HL < HL_threshold_1),]
    if (nrow(df1_p_val_hl) != 0) {
        df1_p_val_hl_p <-
            df1_p_val_hl[which(df1_p_val_hl$p_value_HL <= p_value_hl),]
        if (nrow(df1_p_val_hl_p) != 0) {
            p2 <- my_segment_NS(
                p2,
                data = df1_p_val_hl_p,
                "HL*",
                y = 0,
                yend = 3,
                dis = 10,
                ytext = 3.4,
                color = HL_threshold_color_1,
                linetype = "dashed",
                fontface = fontface
            )
        }
    }
    
    #Select rows with velocity_ratio event and draw a line
    df1_v <- data_p[!duplicated(data_p$velocity_ratio), ]
    df1_v <-
        df1_v[which(df1_v$velocity_ratio < vel_threshold), ]
    if (nrow(df1_v) != 0) {
        p3 <-
            my_segment_NS(
                p3,
                data = df1_v,
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
    ##########################title positive strand#######################
    if (i == 1) {
        p1 <- p1 +
            ggtitle(
                paste0(
                    "ID: ",
                    df1_1$ID[1],
                    "-",
                    last(df1_1$ID),
                    "; ",
                    "FC*: ",
                    "significant t-test of two consecutive segments",
                    ";
                       Term: termination (",
                    forggtitle_ter,
                    "), NS: new start (",
                    forggtitle_NS,
                    "), PS: pausing site (",
                    length(which(data_p$pausing_site == "+")),
                    "), iTSS_I: internal starting site (",
                    length(which(data_p$iTSS_I == "+")),
                    "), (*): p_value below 0.05; ",
                    "TI: transcription interferance."
                )
            ) +
            theme(plot.title = element_text(size = 6, hjust = .5))
    } else{
        p1 <- p1 +
            ggtitle(
                paste0(
                    "ID: ",
                    df1_1$ID[1],
                    "-",
                    last(df1_1$ID),
                    "; Term: termination (",
                    forggtitle_ter,
                    "), NS: new start (",
                    forggtitle_NS,
                    "), PS: pausing site (",
                    length(which(data_p$pausing_site == "+")),
                    "), iTSS_I: internal starting site (",
                    length(which(data_p$iTSS_I == "+")),
                    ")"
                )
            ) +
            theme(plot.title = element_text(size = 6, hjust = .5))
    }
    p <- list(p1, p2, p3)
    return(p)
}
