#' rifi_visualization: plots all the data with fragments and events from both
#' strands.
#' rifi_visualization: plots the whole genome with genes, transcription units
#' (TUs), delay, half-life (HL), intensity fragments
#' features, events, velocity, annotation, coverage if available.
#' rifi_visualization uses several functions to plot the genes including
#' as-RNA and ncRNA and TUs as segments.
#' The function plots delay, HL and intensity fragments with statistical t-test
#' between the neighboring fragment,
#' significant t-test is assigned with '*'.
#' t-test and Manova statistical test are also depicted as '*'.
#' The functions used are:
#' strand_selection: check if data is stranded and arrange by position.
#' splitGenome_function: splits the genome into fragments
#' indice_function: assign a new column to the data to distinguish between
#' fragments, outliers from delay or HL or intensity.
#' TU_annotation: designs the segments border for the genes and TUs annotation
#' gene_annot_function: it requires gff3 file, returns a dataframe adjusting
#' each fragment according to its annotation. It allows as well the plot of
#' genes and TUs shared into two pages
#' label_log2_function: used to add log scale to intensity values.
#' label_square_function: used to add square scale to coverage values.
#' coverage_function: this function is used only in case of coverage is
#' available.
#' secondaryAxis: adjusts the half-life or delay to 20 in case of the dataframe
#' row numbers is equal to 1 and
#' the half-life or delay exceed the limit, they are plotted with different
#' shape and color.
#' add_genomeBorders: when the annotated genes are on the borders, they can
#' not be plotted, therefore the region was split
#' in 2 adding the row corresponding to the split part to the next annotation
#' (i + 1) except for the first page.
#' my_arrow: creates an arrow for the annotation.
#' arrange_byGroup: selects the last row for each segment and add 40 nucleotides
#'  in case of negative strand for a nice plot.
#' regr: plots the predicted delay from linear regression if the data is on
#' negative strand
#' meanPosition: assign a mean position for the plot.
#' delay_mean: adds a column in case of velocity is NA or equal to 60.
#' The mean of the delay is calculated outliers.
#' my_segment_T: plots terminals and pausing sites labels.
#' my_segment_NS: plots internal starting sites 'iTSS'.
#' min_value: returns minimum value for event plots in intensity plot.
#' velocity_fun: function for velocity plot
#' limit_function: for values above 10 or 20 in delay and hl. Limit of the axis
#' is set differently. y-axis limit is applied only if we have
#' more than 3 values above 10 and lower or equal to 20. An exception is added
#' in case a dataframe has less than 3 rows and 1
#' or more values are above 10, the rest of the values above 20 are adjusted to
#' 20 on "secondaryAxis" function.
#' empty_boxes: used only in case the dataframe from the positive strand is not
#' empty, the TU are annotated.
#' function_TU_arrow: used to avoid plotting arrows when a TU is split into two
#' pages.
#' terminal_plot_lm: draws a linear regression line when terminal outliers have
#' an intensity above a certain
#' threshold and are consecutive. Usually are smallRNA (ncRNA, asRNA).
#' slope_function: replaces slope lower than 0.0009 to 0.
#' velo_function: replaces infinite velocity with NA.
#' plot the coverage of RNA_seq in exponential phase growth
#'
#' @param data dataframe: the probe based dataframe.
#' @param genomeLength integer: genome length output of gff3_preprocess
#' function.
#' @param annot dataframe: the annotation file.
#' @param coverage integer: in case the coverage is available.
#' @param chr_fwd string object: coverage of the forward strand.
#' @param chr_rev string object: coverage of the reverse strand.
#' @param region dataframe: gff3 features of the genome.
#' @param color_region string vector: vector of colors.
#' @param color_TU string: TU colors
#' @param fontface integer: value assigning labels font
#' @param color_text.1 string: TU color text
#' @param color_text.2 string: genes color text
#' @param size_tu integer: TU size
#' @param size_locusTag integer: locus_tag size
#' @param Limit integer: value for y-axis limit.
#' @param shape integer: value for shape.
#' @param shape_outlier integer: value for outlier shape.
#' @param col_outiler string: outlier color.
#' @param color_TU string. TU color
#' @param limit_intensity integer: intensity limit if applicable.
#' @param face string: label font.
#' @param tick_length integer: value for ticks.
#' @param arrow.color string: arrows color.
#' @param minVelocity integer: threshold to fix the minimum of velocity.
#' @param medianVelocity integer: threshold to fix the maximum of velocity.
#' @param threshold_intensity integer: Threshold for intensity to plot terminals
#' @param col_above20 string: color for probes/bin above value 20.
#' @param fontface integer: font type
#' @param shape_above20 integer: shape for probes/bins above value 20.
#' @param axis_text_y_size integer: text size for y-axis.
#' @param axis_title_y_size integer: title size for y-axis.
#' @param Alpha integer: color transparency degree.
#' @param size_gene integer: font size for gene annotation.
#' @param col_coverage integer: color for coverage plot.
#' @param TI_threshold integer: threshold for TI between two fragments in case
#' the TI termination factor drops from the first segment to the second,
#' default 1.1.
#' @param p_value_TI integer: p_value of TI fragments selected to be plotted,
#' default 0.05.
#' @param p_value_manova integer: p_value of manova test fragments to plot,
#' default 0.05.
#' @param p_value_int integer: p_value of intensity fragments fold-change to
#' plot, default 0.05.
#' @param p_value_hl integer: p_value of half_life fragments fold-change to
#' plot, default 0.05.
#' @param p_value_event integer: p_value of t-test from pausing site and
#' iTSS_I events to plot, default 0.05.
#' @param termination_threshold integer: threshold for termination to plot,
#' default .8.
#' @param HL_threshold integer: threshold for HL fold change selected to plot,
#' default 20.
#' @param vel_threshold integer: threshold for velocity ratio selected to plot,
#' default 200.
#' @param iTSS_threshold integer: threshold for iTSS_II selected to plot,
#' default 1.2.
#' @param event_duration integer: threshold for pausing sites and iTSS_I
#' selected to plot, default plot all events selecting the maximum value.
#' @param HL_threshold_color string: color for HL fold change plot
#' @param vel_threshold_color string: color for velocity ratio plot
#' @param ps_color string: color for pausing site plot
#' @param iTSS_I_color string: color for iTSS_I plot
#'
#' @return The visualization.
#'
#' @examples
#' data(stats_minimal)
#' data(annot_g_minimal)
#' rifi_visualization(data = stats_minimal, genomeLength = annot_g_minimal[[2]],
#' annot = annot_g_minimal[[1]], coverage = 0, chr_fwd = NA, chr_rev = NA,
#' region = c("CDS","asRNA","5'UTR","ncRNA","3'UTR","tRNA"),
#' color_region = c("grey0", "red", "blue", "orange", "yellow", "green",
#' "white", "darkseagreen1", "grey50", "black"),
#' color_text.1 = "grey0", color_text.2 = "black", color_TU = "blue",
#' size_tu = 1.6, size_locusTag = 1.6, size_gene = 1.6, Limit = 10,
#' shape=22, col_outiler = "grey50", Alpha=0.5,
#' col_coverage = "grey", shape_outlier = 13, limit_intensity = NA,
#' face="bold", tick_length = .3, arrow.color = "darkseagreen1",
#' minVelocity = 3000, medianVelocity = 6000, threshold_intensity = 4000,
#' col_above20 = "#00FFFF", fontface = "plain", shape_above20 = 14,
#' axis_text_y_size = 3, axis_title_y_size = 6, TI_threshold = 1.1,
#' p_value_TI=0.05, p_value_manova = 0.05, termination_threshold = 1,
#' iTSS_threshold = 1.01, p_value_int = 0.05, p_value_event = 0.05,
#' p_value_hl = 0.05, event_duration = 40, HL_threshold=20,
#' vel_threshold = 200, HL_threshold_color="green4",
#' vel_threshold_color="grey52", ps_color="orange", iTSS_I_color="blue")
#'
#' @export

rifi_visualization <-
  function(data,
           genomeLength,
           annot,
           coverage = 0,
           chr_fwd = NA,
           chr_rev = NA,
           region = c("CDS", "asRNA", "5'UTR", "ncRNA", "3'UTR", "tRNA"),
           color_region = c(
             "grey0",
             "red",
             "blue",
             "orange",
             "yellow",
             "green",
             "white",
             "darkseagreen1",
             "grey50",
             "black"
           ),
           color_text.1 = "grey0",
           color_text.2 = "black",
           color_TU = "blue",
           size_tu = 1.6,
           size_locusTag = 1.6,
           size_gene = 1.6,
           Limit = 10,
           shape = 22,
           col_outiler = "grey50",
           Alpha = 0.5,
           col_coverage = "grey",
           shape_outlier = 13,
           limit_intensity = NA,
           face = "bold",
           tick_length = .3,
           arrow.color = "darkseagreen1",
           minVelocity = 3000,
           medianVelocity = 6000,
           threshold_intensity = 4000,
           col_above20 = "#00FFFF",
           fontface = "plain",
           shape_above20 = 14,
           axis_text_y_size = 3,
           axis_title_y_size = 6,
           TI_threshold = 1.1,
           p_value_TI = 0.05,
           p_value_manova = 0.05,
           termination_threshold = 1,
           iTSS_threshold = 1.01,
           p_value_int = 0.05,
           p_value_event = 0.05,
           p_value_hl = 0.05,
           event_duration = 40,
           HL_threshold = 20,
           vel_threshold = 200,
           HL_threshold_color = "green4",
           vel_threshold_color = "grey52",
           ps_color = "orange",
           iTSS_I_color = "blue") {
    ##########################data preparation##################################
    input <- event_dataframe(data, data_annotation = annot)
    #I. add coverage if its available from RNA-seq
    tmp <-
      coverage_function(coverage = coverage,
                        chr_fwd = chr_fwd,
                        chr_rev = chr_rev)
    if (!is.na(tmp)) {
      tmp.c1 <- strand_selection(tmp, "+")
      tmp.c2 <- strand_selection(tmp, "-")
    }
    #II. input for the main features split into 2 data frames according to
    #strand orientation
    tmp.1 <- strand_selection(data, "+")
    tmp.2 <- strand_selection(data, "-")
    #replace infinitive in velocity fragment with NA
    tmp.1 <- velo_function(tmp.1)
    tmp.2 <- velo_function(tmp.2)
    #replace slope lower than 0.0009 to 0
    tmp.1 <- slope_function(tmp.1)
    tmp.2 <- slope_function(tmp.2)
    #III. split the genome into fragments
    gLength <- seq_len(genomeLength)
    names(gLength) <- seq_along(gLength)
    fl <- floor(gLength / 10000)
    frag <- splitGenome_function(x = fl, gLength = gLength)
    #################################plot###############################
    #IV. the general plot
    pdf.options(
      reset = TRUE,
      onefile = TRUE,
      width = 8,
      height = 5.3
    )
    pdf("genome_fragments.pdf")
    an.newLine <- data.frame()
    suppressWarnings(for (i in seq_len(length(frag) - 1)) {
      p <- list()
      print(i)
      if (i == 1) {
        frag[i] <- 0
      } else if (i == (length(frag) - 1)) {
        #to have homogeneous annotation scaling, 10000 is added to the last
        #frag vector.
        frag[i + 1] <- frag[i] + 10000
      }
      ###########################data adjustment###########################
      #define the main dataframe with segments positive strand df1, negative
      #strand df2
      df1 <-
        tmp.1[between(tmp.1$position, frag[i], frag[c(i + 1)]),]
      df2 <-
        tmp.2[between(tmp.2$position, frag[i], frag[c(i + 1)]),]
      df1_1 <- df1[!is.na(df1$ID),]
      #avoid plot empty pages in case of small data
      if (nrow(df1) == 0 & nrow(df2) == 0 & nrow(data) < 10000) {
        next ()
      }
      #an is the annotation dataframe upon the position on the plot, its used
      # to loop into exactly the number of region contained in the gff3
      an <- annot[between(annot$start, frag[i], frag[c(i + 1)]),]
      an <- an[!duplicated(an),]
      #in case of no data nor annotation are available
      if (nrow(an) == 0 & nrow(df1) == 0 & nrow(df2) == 0) {
        next ()
      }
      ##########################annotation section#########################
      if (nrow(an.newLine) != 0) {
        an <- add_row(an.newLine, an)
      }
      if (nrow(an) != 0 & last(an$end) > frag[i + 1]) {
        an.newLine <- add_genomeBorders(data = an,
                                        frag = frag,
                                        i = i)
        firstValue <- an$start[1]
        lastValue <- last(an$end)
        dif <- lastValue - frag[i + 1]
        an[nrow(an), "end"] <- frag[i + 1]
      }
      an.p <- an %>% filter(an$strand == "+")
      an.p <- as.character(unique(an.p$region))
      an.m <- an %>% filter(an$strand == "-")
      an.m <- as.character(unique(an.m$region))
      breaks <- seq(frag[i], frag[i + 1], 1000)
      df.annt <- data.frame()
      segment_data_t.1 <- data.frame()
      segment_data_g.1 <- data.frame()
      segment_data_g.2 <- data.frame()
      segment_data_g.p1 <- data.frame()
      segment_data_g.p2 <- data.frame()
      #TU_annotation and gene_annot_function functions are used to annotate
      #TUs and genes on the positive strand.
      tryCatch({
        #only in case the dataframe from the positive strand is not empty,
        #the TU are annotated
        if (nrow(df1) != 0) {
          segment_data_t.1 <- TU_annotation(df1, "TU", 3, 6, color_TU)
          segment_data_t.1 <-
            segment_data_t.1[grep("_NA", segment_data_t.1$annotation,
                                  invert = TRUE),]
          if (nrow(segment_data_t.1) == 0) {
            segment_data_t.1 <- empty_boxes(
              ystart = 3,
              yend = 6,
              frag = frag,
              i = i
            )
          } else{
            df.annt <- rbind(df.annt, segment_data_t.1)
            # Function to find TUs split into 2 pages
            tu_border <-
              function_TU_arrow(df1, tmp.1, frag = frag, i = i)
            # Eliminate the TU split on the border to avoid arrow plot
            if (length(tu_border) != 0) {
              segment_data_t.1.tuLeft <-
                segment_data_t.1[-which(tu_border ==
                                          segment_data_t.1$annotation),]
            } else{
              segment_data_t.1.tuLeft <- segment_data_t.1
            }
          }
        }
        for (k in seq_along(an.p)) {
          segment_data_g.1 <-
            gene_annot_function(
              pos.1 = frag[i],
              pos.2 = frag[i + 1],
              yint = 11.8,
              yend = 15.8,
              annot = an,
              Strand = "+",
              reg = an.p[k],
              col = color_region[which(an.p[k] == region)]
            )
          segment_data_g.2 <-
            gene_annot_function(
              pos.1 = frag[i],
              pos.2 = frag[i + 1],
              yint = 7.4,
              yend = 11.4,
              annot = an,
              Strand = "+",
              reg = an.p[k],
              col = color_region[which(an.p[k] == region)]
            )
          segment_data_g.p1 <-
            rbind(segment_data_g.p1, segment_data_g.1)
          segment_data_g.p2 <-
            rbind(segment_data_g.p2, segment_data_g.2)
        }
      }, warning = function(war) {
        
      },
      error = function(err) {
        
      })
      segment_data_t.2 <- data.frame()
      segment_data_g.1 <- data.frame()
      segment_data_g.2 <- data.frame()
      segment_data_g.m1 <- data.frame()
      segment_data_g.m2 <- data.frame()
      #TU and genes annotation in the negative strand
      tryCatch({
        if (nrow(df2) != 0) {
          segment_data_t.2 <- TU_annotation(df2, "TU",-3,-6, color_TU)
          segment_data_t.2 <-
            segment_data_t.2[grep("_NA", segment_data_t.2$annotation,
                                  invert = TRUE),]
          if (nrow(segment_data_t.2) == 0) {
            segment_data_t.2 <- empty_boxes(-3,-6, frag = frag, i = i)
          } else{
            #function to find TUs split into 2 pages
            tu_border <-
              function_TU_arrow(df2, tmp.2, frag = frag, i = i)
            # Eliminate the TU split on the border to avoid arrow plot
            if (length(tu_border) != 0) {
              segment_data_t.2.tuLeft <-
                segment_data_t.2[-which(tu_border ==
                                          segment_data_t.2$annotation),]
            } else{
              segment_data_t.2.tuLeft <- segment_data_t.2
            }
          }
        }
        for (k in seq_along(an.m)) {
          segment_data_g.1 <-
            gene_annot_function(
              pos.1 = frag[i],
              pos.2 = frag[i + 1],
              yint = -11.8,
              yend = -15.8,
              annot = an,
              Strand = "-",
              reg = an.m[k],
              col = color_region[which(an.m[k] == region)]
            )
          segment_data_g.2 <-
            gene_annot_function(
              pos.1 = frag[i],
              pos.2 = frag[i + 1],
              yint = -7.4,
              yend = -11.4,
              annot = an,
              Strand = "-",
              reg = an.m[k],
              col = color_region[which(an.m[k] == region)]
            )
          segment_data_g.m1 <-
            rbind(segment_data_g.m1, segment_data_g.1)
          segment_data_g.m2 <-
            rbind(segment_data_g.m2, segment_data_g.2)
        }
      }, warning = function(war) {
        
      },
      error = function(err) {
        
      })
      df.annt <-
        rbind(
          segment_data_t.1,
          segment_data_t.2,
          segment_data_g.p1,
          segment_data_g.p2,
          segment_data_g.m1,
          segment_data_g.m2
        )
      if (nrow(df.annt) != 0) {
        p7 <- ggplot(data = df.annt) +
          scale_x_continuous(limits = c(frag[i], frag[i + 1]),
                             breaks = breaks) +
          scale_y_continuous(limits = c(-15.8, 15.8), expand = c(0, 0)) +
          geom_hline(yintercept = 0) +
          theme_minimal(base_size = 11) +
          theme(
            legend.title = element_blank(),
            axis.title.x = element_blank(),
            legend.position = "none",
            axis.ticks.y =  element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_blank(),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            plot.title = element_blank(),
            panel.border = element_blank()
          )
        if (nrow(segment_data_t.1) != 0) {
          p7 <- p7 +
            annotate(
              "segment",
              x = frag[i + 1],
              xend = frag[i],
              y = unique((
                segment_data_t.1$ystart + segment_data_t.1$yend
              ) / 2),
              yend = unique((
                segment_data_t.1$ystart + segment_data_t.1$yend
              ) / 2),
              size = .1,
              arrow = my_arrow(1, "open")
            )
        }
        if (nrow(segment_data_t.2) != 0) {
          p7 <- p7 +
            annotate(
              "segment",
              x = frag[i],
              xend = frag[i + 1],
              y = unique((
                segment_data_t.2$ystart + segment_data_t.2$yend
              ) / 2),
              yend = unique((
                segment_data_t.2$ystart + segment_data_t.2$yend
              ) / 2),
              size = .1,
              arrow = my_arrow(1, "open")
            )
        }
        if (nrow(df.annt) != 0) {
          p7 <- p7 +
            geom_rect(
              data = df.annt,
              aes(
                xmin = get('xstart'),
                ymin = get('ystart'),
                xmax = get('xend'),
                ymax = get('yend')
              ),
              alpha = Alpha,
              fill = df.annt$color
            ) +
            geom_text(
              data = df.annt[grep("TU", df.annt$annotation),],
              aes((get('xstart') + get('xend')) / 2, (get('ystart') +
                                                        get('yend')) / 2,
                  label = get('annotation')),
              size = size_tu,
              fontface = fontface,
              color = color_text.1,
              check_overlap = TRUE
            )
        }
        #avoid 5'UTR and 3'UTR label
        df.annt.1 <-
          df.annt[grep("UTR", df.annt$region, invert = TRUE),]
        if (nrow(df.annt.1) != 0) {
          p7 <- p7 +
            geom_text(
              data = df.annt.1[which(df.annt.1$ystart %in% c(7.4,-7.4)),],
              aes((get('xstart') + get('xend')) / 2, (get('ystart') +
                                                        get('yend')) / 2,
                  label = get('annotation')),
              size = size_locusTag,
              fontface = fontface,
              color = color_text.1,
              check_overlap = TRUE
            ) +
            geom_text(
              data = df.annt.1[which(df.annt.1$ystart %in% c(11.8,-11.8)),],
              aes((get('xstart') + get('xend')) / 2, (get('ystart') +
                                                        get('yend')) / 2,
                  label = get('gene')),
              size = size_gene,
              fontface = fontface,
              color = color_text.2,
              check_overlap = TRUE
            )
        }
        if (nrow(segment_data_t.1) != 0) {
          p7 <- p7 +
            annotate(
              "segment",
              x = unique(segment_data_t.1.tuLeft$xstart) - 10,
              xend = unique(segment_data_t.1.tuLeft$xstart) - 10,
              y = unique(segment_data_t.1.tuLeft$ystart),
              yend = unique(segment_data_t.1.tuLeft$yend) + .3,
              size = .2
            ) +
            annotate(
              "segment",
              xend = unique(segment_data_t.1.tuLeft$xstart) - 5,
              x = unique(segment_data_t.1.tuLeft$xstart) + 60,
              y = unique(segment_data_t.1.tuLeft$yend) + .3,
              yend = unique(segment_data_t.1.tuLeft$yend) + .3,
              size = .2,
              arrow = my_arrow(.4, "open")
            )
        }
        if (nrow(segment_data_t.2) != 0) {
          p7 <- p7 +
            annotate(
              "segment",
              x = unique(segment_data_t.2.tuLeft$xend) + 10,
              xend = unique(segment_data_t.2.tuLeft$xend) + 10,
              y = unique(segment_data_t.2.tuLeft$ystart),
              yend = unique(segment_data_t.2.tuLeft$yend) - .3,
              size = .2
            ) +
            annotate(
              "segment",
              xend = unique(segment_data_t.2.tuLeft$xend) + 5,
              x = unique(segment_data_t.2.tuLeft$xend) - 60,
              y = unique(segment_data_t.2.tuLeft$yend) - .3,
              yend = unique(segment_data_t.2.tuLeft$yend) - .3,
              size = .2,
              arrow = my_arrow(.4, "open")
            )
        }
      }
      #########################empty data positive strand#####################
      #in case the dataframe from the positive strand is empty, a dataframe is
      #created to be plotted for the genome annotation.
      if (nrow(df1) == 0) {
        print(paste0(i, ": no data on positive strand"))
        df1_f <- data.frame(matrix(NA, nrow = 1, ncol = ncol(df1)))
        colnames(df1_f) <- colnames(df1)
        df1_f$ID <- "ID_fake"
        df1_f$position <- frag[i]
        df1_f$strand <- "+"
        df1_f$delay <- .01
        df1_f$half_life <- .01
        df1_f$intensity <- 1000
        Title <- NA
        #add arrow and virtual boxes to adjust the scale in ggplot
        if (nrow(segment_data_t.1) == 0) {
          segment_data_t.1 <- empty_boxes(2, 4, frag = frag, i = i)
        } else{
          p7 <- p7 +
            geom_rect(
              data = segment_data_t.1,
              aes(
                xmin = get('xstart'),
                ymin = get('ystart'),
                xmax = get('xend'),
                ymax = get('yend')
              ),
              alpha = 0.2,
              fill = color_TU
            )
        }
        p1 <-
          ggplot(df1_f, aes(x = get('position'), y = get('intensity'))) +
          scale_x_continuous(limits = c(frag[i], frag[i + 1])) +
          labs(y = "Intensity [A.U]") +
          scale_y_continuous(
            trans = 'log2',
            labels = label_log2_function,
            sec.axis = sec_axis( ~ . * 1, name = "Coverage",
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
        if (nrow(df2) != 0) {
          p1 <- p1 +
            ggtitle(
              paste0(
                "ID: ",
                df2$ID[1],
                "-",
                last(df2$ID),
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
            sec.axis = sec_axis( ~ . * 1, name = "Delay [min]", breaks =
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
      }
      ##########################positive strand############################
      #This part of the script runs only in case of the dataframe from the
      #positive strand is not empty
      if (nrow(df1) != 0) {
        ########################ggplot positive strand#######################
        #first plot for intensity segments
        p1 <-
          ggplot(df1, aes(x = get('position'), y = get('intensity'))) +
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
        if (length(grep("TI", df1$flag)) != 0) {
          #grep only flagged probes/bins with TI
          df1_ti_a <- df1[grep("TI", df1$flag),]
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
                nrow(df1_ti %>% group_by(get('TI_termination_fragment'))) ==
                length(unique(df1_ti$TI_termination_fragment)) |
                all(df1_ti$TI_mean_termination_factor == 0)) {
              p1 <- p1
              #geom_step runs only with more observations for each group exists,
              #in opposite case only geom_line from intensity is plotted.
            } else{
              # add a value corresponding to the TI factor relative to the
              # intensity values the value plotted dividing the intensity
              # mean by (1 - TI factor)
              df1_ti$value_d <-
                df1_ti$intensity_mean_fragment /
                (1 - df1_ti$TI_termination_factor)
              df1_ti$value_l <- df1_ti$intensity_mean_fragment /
                (1 - df1_ti$TI_mean_termination_factor)
              if (length(unique(df1_ti$TI_termination_fragment)) == 1) {
                p1 <- p1 +
                  geom_line(data = df1_ti,
                            aes(
                              y = get('value_l'),
                              group = get('TI_mean_termination_factor')
                            ),
                            size = .2) +
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
                  indice_function(df1_ti, "TI_termination_fragment")
                df1_ti <- df1_ti %>% filter(get('indice') == 1)
                p1 <- p1 +
                  geom_step(data = df1_ti,
                            aes(
                              y = get('value_l'),
                              group = get('TI_mean_termination_factor')
                            ),
                            size = .2) +
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
                  df1.1_ti[which(df1.1_ti$TI_termination_fragment %in% TIs),]
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
                           filter(get('p_value_TI') < p_value_TI)) != 0) {
                    p1 <- p1 +
                      geom_text(
                        data = df1.1_ti,
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
                  } else if (nrow(df1.1_ti %>% filter(get('p_value_TI') >
                                                      p_value_TI)) != 0) {
                    p1 <- p1 +
                      geom_text(
                        data = df1.1_ti,
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
        #######################segment plot positive strand###################
        #first plot for half-life segments
        #select segments without outliers
        df1 <- indice_function(df1, "HL_fragment")
        # increase the limit to 20 in case 3 or more probes/bins have a HL
        # above 10
        Limit_h_df1 <- limit_function(df1, "half_life", ind = 1)
        if (Limit_h_df1 == 20) {
          Breaks_h <- seq(0, Limit_h_df1, by = 4)
        } else{
          Breaks_h <- seq(0, Limit_h_df1, by = 2)
        }
        df1.h <- secondaryAxis(df1, "half_life", ind = 1)
        #in case only one bin is available and the HL is above 20
        if (all(df1$half_life > 20)) {
          df1$half_life <- 20
        }
        p2 <-
          ggplot(df1, aes(x = get('position'), y = get('half_life'))) +
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
        if (length(unique(df1$half_life)) == 1) {
          if (is.na(unique(df1$half_life))) {
            p2 <- p2 +
              geom_blank() +
              scale_y_continuous(
                limits = c(0, Limit_h_df1),
                breaks = Breaks_h,
                sec.axis = sec_axis( ~ . * 1, name = "Half-life [min]", breaks =
                                       Breaks_h)
              )
          }
        }
        #select segments without outliers
        df1 <- indice_function(df1, "delay_fragment")
        #increase the limit to 20 in case 3 or more probes/bins have a delay
        #above 10
        Limit_df1 <- limit_function(df1, "delay", ind = 1)
        if (Limit_df1 == 20) {
          Breaks_d <- seq(0, Limit_df1, by = 4)
        } else{
          Breaks_d <- seq(0, Limit_df1, by = 2)
        }
        #first plot for delay segments
        df1.d <- secondaryAxis(df1, "delay", ind = 1)
        #in case only one bin is available and the delay is above 20
        if (all(df1$delay > 20)) {
          df1$delay <- 20
        }
        p3 <-
          ggplot(df1, aes(x = get('position'), y = get('delay'))) +
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
        if (length(unique(df1$delay)) == 1) {
          if (is.na(unique(df1$delay))) {
            p3 <- p3 +
              scale_y_continuous(
                limits = c(0, Limit_df1),
                breaks = Breaks_d,
                sec.axis = sec_axis( ~ . * 1, name = "Delay [min]", breaks =
                                       Breaks_d)
              )
          }
        }
        #add the segments to the plot base and check
        #intensity is plotted independently if delay/HL data are present or not.
        if (length(na.omit(df1$delay)) == 0) {
          p1 <- p1 +
            geom_point(size = .5)
          p2 <- p2
          p3 <- p3
        } else if (length(na.omit(df1$delay)) < 3) {
          p1 <- p1 +
            geom_point(size = .5)
          p2 <- p2 +
            geom_point(size = .5)
          p3 <- p3 +
            geom_point(size = .5)
        } else{
          ######################intensity plot positive strand################
          #add a reference to outliers probes or bins
          df1 <- indice_function(df1, "intensity_fragment")
          #intensities/HL and delay with NA are not plotted, high intensities
          #are plotted with a different shape,
          #outliers probes/bins are plotted with a different color, each
          #segment has a color, a line is added for each segment
          #indicating the mean. "FC*" indicate significant t-test and each
          #segment is labeled.
          if (nrow(df1 %>% filter(get('indice') == 1)) != 0) {
            p1 <- p1 +
              geom_point(data = df1 %>% filter(get('indice') == 1),
                         aes(col = get('intensity_fragment')),
                         size = .5)
            #eliminate outliers probes or bins
            df1_wo <- df1[which(df1$indice == 1),]
            #assign a mean position for the plot
            df1_wo <- meanPosition(df1_wo, "intensity_fragment")
            if (length(df1_wo$intensity_fragment) != 0) {
              p1 <- p1 +
                geom_line(data = df1 %>% filter(get('indice') == 1),
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
          if (nrow(df1 %>% filter(get('indice') == 2)) != 0) {
            p1 <- p1 +
              geom_point(
                data = df1 %>% filter(get('indice') == 2),
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
          df1 <- indice_function(df1, "HL_fragment")
          if (nrow(df1 %>% filter(get('indice') == 1)) != 0) {
            p2 <- p2 +
              geom_point(data = df1 %>% filter(get('indice') == 1),
                         aes(col = get('HL_fragment')),
                         size = .5)
            #eliminate outliers probes or bins
            df1_wo <- df1[which(df1$indice == 1),]
            #assign a mean position for the plot
            df1_wo <- meanPosition(df1_wo, "HL_fragment")
            if (length(df1_wo$HL_fragment) != 0) {
              p2 <- p2 +
                geom_line(data = df1 %>% filter(get('indice') == 1),
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
          if (nrow(df1 %>% filter(get('indice') == 2)) != 0) {
            p2 <- p2 +
              geom_point(
                data = df1 %>% filter(get('indice') == 2),
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
          df1 <- indice_function(df1, "delay_fragment")
          #plot delay bins
          if (nrow(df1 %>%
                   filter(get('indice') == 1)) != 0) {
            p3 <- p3 +
              geom_point(data = df1 %>%
                           filter(get('indice') == 1),
                         aes(col = get('delay_fragment')),
                         size = .5)
          }
          #plot delay outliers
          if (nrow(df1 %>%
                   filter(get('indice') == 2)) != 0) {
            p3 <- p3 +
              geom_point(
                data = df1 %>%
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
          if (length(na.omit(df1$delay)) > 2) {
            df.c <- regr(df1, ind = 1, data = data)
            #delay_mean column is added to draw a line in case of velocity
            #is NA or equal to 60.
            #the mean of the delay is calculated excluding outliers
            df1$delay_mean <- delay_mean(df1, "delay", ind = 1)
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
            if (nrow(df1 %>%
                     filter(get('indice') == 1) %>%
                     filter(get('slope') == 0)) > 2) {
              p3 <- p3 +
                geom_line(
                  data = df1 %>%
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
            df1[grep(paste0("\\D_\\d+", "$"), df1$delay_fragment),]
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
        if (length(which(!is.na(df1$synthesis_ratio))) != 0) {
          df1_syR <- df1[which(!is.na(df1$synthesis_ratio)),]
          #in case last position matches with an event which needs to be
          #on the next page of the plot.
          if (last(df1_syR$position) == frag[c(i + 1)]) {
            fc_seg <-
              df1_syR[which(df1_syR$position ==
                              frag[c(i + 1)]), "FC_HL_intensity_fragment"]
            df1_syR[which(df1_syR$FC_HL_intensity_fragment == fc_seg),
                    c("synthesis_ratio_event", "p_value_Manova")] <-
              NA
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
            if (length(which(df1_syR_T$p_value_Manova < p_value_manova)) != 0) {
              df1_syR_T.m <-
                df1_syR_T %>%
                filter(get('p_value_Manova') < p_value_manova)
              p2 <-
                my_segment_T(
                  p2,
                  data = df1_syR_T.m,
                  "Ter*",
                  y = 0,
                  yend = 2.5,
                  dis = 50,
                  ytext = 2.9,
                  color = 2,
                  linetype = "dotted",
                  df = "termination",
                  fontface = fontface
                )
            }
            if (length(which(df1_syR_T$p_value_Manova > p_value_manova)) != 0) {
              df1_syR_T.t <-
                df1_syR_T %>%
                filter(get('p_value_Manova') > p_value_manova)
              p2 <-
                my_segment_T(
                  p2,
                  data = df1_syR_T.t,
                  "Ter",
                  y = 0,
                  yend = 2.5,
                  dis = 50,
                  ytext = 2.9,
                  color = 2,
                  linetype = "dotted",
                  df = "termination",
                  fontface = fontface
                )
            }
            if (length(which(is.na(df1_syR_T$p_value_Manova))) != 0) {
              df1_syR_T.t <- df1_syR_T %>%
                filter(is.na(get('p_value_Manova')))
              p2 <-
                my_segment_T(
                  p2,
                  data = df1_syR_T.t,
                  "Ter",
                  y = 0,
                  yend = 2.5,
                  dis = 50,
                  ytext = 2.9,
                  color = 2,
                  linetype = "dotted",
                  df = "termination",
                  fontface = fontface
                )
            }
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
            if (length(which(df1_syR_T$p_value_Manova < p_value_manova)) != 0) {
              df1_syR_T.m <-
                df1_syR_T %>%
                filter(get('p_value_Manova') < p_value_manova)
              p2 <-
                my_segment_NS(
                  p2,
                  data = df1_syR_T.m,
                  "NS*",
                  y = 0,
                  yend = 2.5,
                  dis = 10,
                  ytext = 2.9,
                  color = "#00FFFF",
                  linetype = "dashed",
                  fontface = fontface
                )
            }
            if (length(which(df1_syR_T$p_value_Manova > p_value_manova)) != 0) {
              df1_syR_T.t <-
                df1_syR_T %>%
                filter(get('p_value_Manova') > p_value_manova)
              p2 <-
                my_segment_NS(
                  p2,
                  data = df1_syR_T.t,
                  "NS",
                  y = 0,
                  yend = 2.5,
                  dis = 10,
                  ytext = 2.9,
                  color = "#00FFFF",
                  linetype = "dashed",
                  fontface = fontface
                )
            }
            if (length(which(is.na(df1_syR_T$p_value_Manova))) != 0) {
              df1_syR_T.t <- df1_syR_T %>%
                filter(is.na(get('p_value_Manova')))
              p2 <-
                my_segment_NS(
                  p2,
                  data = df1_syR_T.t,
                  "NS",
                  y = 0,
                  yend = 2.5,
                  dis = 10,
                  ytext = 2.9,
                  color = "#00FFFF",
                  linetype = "dashed",
                  fontface = fontface
                )
            }
          }
        }
        ####################pausing site positive strand##############
        #plot pausing sites events
        if (nrow(df1 %>%
                 filter(get('pausing_site') == "+")) != 0) {
          #select pausing site
          df1_ps <- df1 %>%
            filter(get('pausing_site') == "+")
          #select pausing site duration
          df1_ps <-
            df1_ps %>%
            filter(na.omit(get('event_duration')) <= event_duration)
          if (nrow(df1_ps %>%
                   filter(get('event_ps_itss_p_value_Ttest')
                          < p_value_event)) != 0) {
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
                color = ps_color,
                linetype = "dashed",
                df = "pausing",
                fontface = fontface
              )
          }
          if (nrow(df1_ps %>%
                   filter(get('event_ps_itss_p_value_Ttest')
                          > p_value_event)) != 0) {
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
                color = ps_color,
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
                color = ps_color,
                linetype = "dashed",
                df = "pausing",
                fontface = fontface
              )
          }
        }
        ####################iTSS_I positive strand###################
        #plot internal starting sites events
        if (nrow(df1 %>%
                 filter(get('iTSS_I') == "+")) != 0) {
          #select iTSS
          df1_itss <- df1 %>%
            filter(get('iTSS_I') == "+")
          #select iTSS duration
          df1_itss <-
            df1_itss %>%
            filter(na.omit(get('event_duration'))
                   <= event_duration)
          if (nrow(df1_itss %>%
                   filter(get('event_ps_itss_p_value_Ttest')
                          < p_value_event)) != 0) {
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
                color = iTSS_I_color,
                linetype = "dotted",
                fontface = fontface
              )
          }
          if (nrow(df1_itss %>%
                   filter(get('event_ps_itss_p_value_Ttest')
                          > p_value_event)) != 0) {
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
                color = iTSS_I_color,
                linetype = "dotted",
                fontface = fontface
              )
          }
          if (nrow(df1_itss %>%
                   filter(is.na('event_ps_itss_p_value_Ttest'))) != 0) {
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
                color = iTSS_I_color,
                linetype = "dotted",
                fontface = fontface
              )
          }
        }
        ####################FC positive strand###################
        #add FC for intensity ratio test if p_value is significant
        #select the last row for each segment and add 40 nucleotides in case
        #of negative strand for a nice plot
        df1 <- indice_function(df1, "intensity_fragment")
        df1_wo <- df1[which(df1$indice == 1),]
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
        #add FC for HL ratio test if p_value is significant
        df1 <- indice_function(df1, "HL_fragment")
        df1_wo <- df1[which(df1$indice == 1),]
        df1_p_val_hl <-  arrange_byGroup(df1_wo, "p_value_HL")
        if (nrow(df1_p_val_hl) != 0) {
          df1_p_val_hl <-
            df1_wo_pvalue[which(df1_wo_pvalue$p_value_HL < p_value_hl),]
          if (nrow(df1_p_val_hl) != 0) {
            p2 <- p2 +
              geom_text(
                data = df1_p_val_hl,
                aes(
                  x = get('position'),
                  y = get('HL_mean_fragment')
                ),
                label = "FC*",
                size = 1,
                check_overlap = TRUE
              )
          }
          df1_p_val_hl <-
            df1_wo_pvalue[which(df1_wo_pvalue$p_value_HL > p_value_hl),]
          if (nrow(df1_p_val_hl) != 0) {
            p2 <- p2 +
              geom_text(
                data = df1_p_val_hl,
                aes(
                  x = get('position'),
                  y = get('HL_mean_fragment')
                ),
                label = "FC",
                size = 1,
                check_overlap = TRUE
              )
          }
        }
        #Select rows with FC_HL event and draw a line
        df1_hl <- arrange_byGroup(df1, "FC_HL")
        df1_hl <- df1_hl[which(df1_hl$FC_HL < HL_threshold),]
        if (nrow(df1_hl) != 0) {
          p2 <-
            my_segment_NS(
              p2,
              data = df1_hl,
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
        #Select rows with velocity_ratio event and draw a line
        df1_v <- arrange_byGroup(df1, "velocity_ratio")
        df1_v <-
          df1_v[which(df1_v$velocity_ratio < vel_threshold),]
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
                length(which(df1$pausing_site == "+")),
                "), iTSS_I: internal starting site (",
                length(which(df1$iTSS_I == "+")),
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
                length(which(df1$pausing_site == "+")),
                "), iTSS_I: internal starting site (",
                length(which(df1$iTSS_I == "+")),
                ")"
              )
            ) +
            theme(plot.title = element_text(size = 6, hjust = .5))
        }
      }
      #########################empty data reverse strand###################
      if (nrow(df2) == 0) {
        print(paste0(i, ": no data on negative strand"))
        df2_f <- data.frame(matrix(NA, nrow = 1, ncol = ncol(df2)))
        colnames(df2_f) <- colnames(df2)
        df2_f$ID <- "ID_fake"
        df2_f$position <- frag[i]
        df2_f$strand <- "-"
        df2_f$delay <- .01
        df2_f$half_life <- .01
        df2_f$intensity <- 1000
        Title <- NA
        # #add empty boxes
        if (nrow(segment_data_t.2) == 0) {
          segment_data_t.2 <- empty_boxes(-3.5,-5.5, frag = frag, i = i)
        } else {
          p7 <- p7 +
            geom_rect(
              data = segment_data_t.2,
              aes(
                xmin = get('xstart'),
                ymin = get('ystart'),
                xmax = get('xend'),
                ymax = get('yend')
              ),
              alpha = Alpha,
              fill = color_TU
            )
        }
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
      }
      #############################reverse strand#############################
      if (nrow(df2) != 0) {
        ##############################ggplot reverse strand#####################
        p4 <-
          ggplot(df2, aes(x = get('position'), y = get('intensity'))) +
          scale_x_continuous(limits = c(frag[i], frag[i + 1])) +
          scale_y_continuous(
            trans = 'log2',
            labels = label_log2_function,
            limits = c(NA, NA),
            sec.axis = sec_axis( ~ . * 1, name = "Coverage",
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
        df2 <- indice_function(df2, "HL_fragment")
        #increase the limit to 20 in case 3 or more probes/bins have a
        #delay above 10
        Limit_h_df2 <- limit_function(df2, "half_life", ind = 1)
        if (Limit_h_df2 == 20) {
          Breaks_h2 <- seq(0, Limit_h_df2, by = 4)
        } else{
          Breaks_h2 <- seq(0, Limit_h_df2, by = 2)
        }
        df2.h <- secondaryAxis(df2, "half_life", ind = 1)
        #in case only one bin is available and the HL is above 20
        if (all(df2$half_life > 20)) {
          df2$half_life <- 20
        }
        p5 <-
          ggplot(df2, aes(x = get('position'), y = get('half_life'))) +
          scale_x_continuous(limits = c(frag[i], frag[c(i + 1)])) +
          scale_y_continuous(
            limits = c(0, Limit_h_df2),
            breaks = Breaks_h2,
            sec.axis = sec_axis( ~ . * 1, name = "Half-life [min]", breaks =
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
        if (length(unique(df2$half_life)) == 1) {
          if (is.na(unique(df2$half_life))) {
            df2$half_life[c(1, 2)] <- c(.1, .2)
            p5 <-
              ggplot(df2, aes(
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
        df2 <- indice_function(df2, "delay_fragment")
        #increase the limit to 20 in case 3 or more probes/bins have a
        #delay above 10
        Limit_df2 <- limit_function(df2, "delay", ind = 1)
        if (Limit_df2 == 20) {
          Breaks_d2 <- seq(0, Limit_df2, by = 4)
        } else{
          Breaks_d2 <- seq(0, Limit_df2, by = 2)
        }
        df2.d <- secondaryAxis(df2, "delay", ind = 1)
        #in case only one bin is available and the delay is above 20
        if (all(df2$delay > 20)) {
          df2$delay <- 20
        }
        p6 <-
          ggplot(df2, aes(x = get('position'), y = get('delay'))) +
          scale_x_continuous(limits = c(frag[i], frag[i + 1])) +
          scale_y_continuous(
            limits = c(0, Limit_df2),
            breaks = Breaks_d2,
            sec.axis = sec_axis( ~ . * 1, name = "Delay [min]", breaks =
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
        if (length(unique(df2$delay)) == 1) {
          if (is.na(unique(df2$delay))) {
            p6 <- p6 +
              scale_y_continuous(
                limits = c(0, Limit_df2),
                breaks = Breaks_d2,
                sec.axis = sec_axis( ~ . * 1, name = "Delay [min]", breaks =
                                       Breaks_d2)
              ) +
              coord_trans(y = 'reverse')
          }
        }
        #############################TI plot reverse strand###################
        #plot transcription interference on the background
        if (length(grep("TI", df2$flag)) != 0) {
          df2_ti_a <- df2[grep("TI", df2$flag),]
          #looping into TUs
          for (l in seq_along(unique(df2_ti_a$TU))) {
            df2_ti <- df2_ti_a[which(df2_ti_a$TU == unique(df2_ti_a$TU)[l]),]
            df2_ti <-
              df2_ti[grep(paste0("\\TI_\\d+", "$"),
                          df2_ti$TI_termination_fragment),]
            df2_ti <-
              df2_ti[!is.na(df2_ti$TI_mean_termination_factor),]
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
                  indice_function(df2_ti, "TI_termination_fragment")
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
                  df2.1_ti[which(df2.1_ti$TI_termination_fragment %in% TIs),]
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
        if (length(na.omit(df2$delay)) == 0) {
          p4 <- p4 +
            geom_point(size = .5)
          p5 <- p5
          p6 <- p6
        } else if (length(na.omit(df2$delay)) < 3) {
          p4 <- p4 +
            geom_point(size = .5)
          p5 <- p5 +
            geom_point(size = .5)
          p6 <- p6 +
            geom_point(size = .5)
        } else{
          ######################intensity plot reverse strand################
          df2 <- indice_function(df2, "intensity_fragment")
          if (nrow(df2 %>%
                   filter(get('indice') == 1)) != 0) {
            p4 <- p4 +
              geom_point(data = df2 %>%
                           filter(get('indice') == 1),
                         aes(col = get('intensity_fragment')),
                         size = .5)
            df2_wo <- df2[which(df2$indice == 1),]
            df2_wo <- meanPosition(df2_wo, "intensity_fragment")
            if (length(df2_wo$intensity_fragment) != 0) {
              p4 <- p4 +
                geom_line(data = df2 %>%
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
          if (nrow(df2 %>%
                   filter(get('indice') == 2)) != 0) {
            p4 <- p4 +
              geom_point(
                data = df2 %>%
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
          df2 <- indice_function(df2, "HL_fragment")
          if (nrow(df2 %>%
                   filter(get('indice') == 1)) != 0) {
            p5 <- p5 +
              geom_point(data = df2 %>%
                           filter(get('indice') == 1),
                         aes(col = get('HL_fragment')),
                         size = .5)
            df2_wo <- df2[which(df2$indice == 1),]
            df2_wo <- meanPosition(df2_wo, "HL_fragment")
            if (length(df2_wo$HL_fragment) != 0) {
              p5 <- p5 +
                geom_line(data = df2 %>%
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
          if (nrow(df2 %>%
                   filter(get('indice') == 2)) != 0) {
            p5 <- p5 +
              geom_point(
                data = df2 %>%
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
          #######################delay plot reverse strand##################
          #plot delay segment and outliers
          df2 <- indice_function(df2, "delay_fragment")
          if (nrow(df2 %>%
                   filter(get('indice') == 1)) != 0) {
            p6 <- p6 +
              geom_point(data = df2 %>%
                           filter(get('indice') == 1),
                         aes(col = get('delay_fragment')),
                         size = .5)
          }
          if (nrow(df2 %>%
                   filter(get('indice') == 2)) != 0) {
            p6 <- p6 +
              geom_point(
                data = df2 %>%
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
          if (length(na.omit(df2$delay)) > 2) {
            df.c <- regr(df2, ind = 1, data = data)
            #delay_mean column is added to draw a line in case of velocity
            #is NA or equal to 60.
            #the mean of the delay is calculated excluding outliers
            df2$delay_mean <- delay_mean(df2, "delay", ind = 1)
            #in case a fragment is split on two pages, each break get a
            #different velocity
            #upon the corresponding delay therefore a segment is checked
            #first if its flaking the borders
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
            if (nrow(df2 %>%
                     filter(get('indice') == 1) %>%
                     filter(get('slope') == 0)) >= 3) {
              p6 <- p6 +
                geom_line(
                  data = df2 %>%
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
            df2[grep(paste0("\\D_\\d+", "$"), df2$delay_fragment),]
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
        if (length(which(!is.na(df2$synthesis_ratio))) != 0) {
          df2_syR <- df2[which(!is.na(df2$synthesis_ratio)),]
          if (length(which(
            df2_syR$synthesis_ratio < termination_threshold &
            !is.na(df2_syR$FC_HL_intensity_fragment)
          )) != 0) {
            df2_syR_T <-
              df2_syR[which(
                df2_syR$synthesis_ratio < termination_threshold &
                  !is.na(df2_syR$FC_HL_intensity_fragment)
              ),]
            df2_syR_T <-
              df2_syR_T[!duplicated(df2_syR_T$FC_fragment_intensity),]
            forggtitle_ter <- nrow(df2_syR_T)
            if (length(which(df2_syR_T$p_value_Manova < p_value_manova)) != 0) {
              df2_syR_T.m <-
                df2_syR_T %>%
                filter(get('p_value_Manova') < p_value_manova)
              p5 <-
                my_segment_T(
                  p5,
                  data = df2_syR_T.m,
                  "Ter*",
                  y = 0,
                  yend = 2.5,
                  dis = 50,
                  ytext = 2.7,
                  color = 2,
                  linetype = "dotted",
                  df = "termination",
                  fontface = fontface
                )
            }
            if (length(which(df2_syR_T$p_value_Manova > p_value_manova)) != 0) {
              df2_syR_T.t <-
                df2_syR_T %>%
                filter(get('p_value_Manova') > p_value_manova)
              p5 <-
                my_segment_T(
                  p5,
                  data = df2_syR_T.t,
                  "Ter",
                  y = 0,
                  yend = 2.5,
                  dis = 50,
                  ytext = 2.7,
                  color = 2,
                  linetype = "dotted",
                  df = "termination",
                  fontface = fontface
                )
            }
            if (length(which(is.na(df2_syR_T$p_value_Manova))) != 0) {
              df2_syR_T.t <- df2_syR_T %>%
                filter(is.na(get('p_value_Manova')))
              p5 <-
                my_segment_T(
                  p5,
                  data = df2_syR_T.t,
                  "Ter",
                  y = 0,
                  yend = 2.5,
                  dis = 50,
                  ytext = 2.7,
                  color = 2,
                  linetype = "dotted",
                  df = "termination",
                  fontface = fontface
                )
            }
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
              ),]
            df2_syR_T <-
              df2_syR_T[!duplicated(df2_syR_T$FC_fragment_intensity),]
            forggtitle_NS <- nrow(df2_syR_T)
            if (length(which(df2_syR_T$p_value_Manova < p_value_manova)) != 0) {
              df2_syR_T.m <-
                df2_syR_T %>%
                filter(get('p_value_Manova') < p_value_manova)
              p6 <-
                my_segment_NS(
                  p6,
                  data = df2_syR_T.m,
                  "NS*",
                  y = 0,
                  yend = 2.5,
                  dis = 10,
                  ytext = 2.7,
                  color = "#00FFFF",
                  linetype = "dashed",
                  fontface = fontface
                )
            }
            if (length(which(df2_syR_T$p_value_Manova > p_value_manova)) != 0) {
              df2_syR_T.t <-
                df2_syR_T %>%
                filter(get('p_value_Manova') > p_value_manova)
              p6 <-
                my_segment_NS(
                  p6,
                  data = df2_syR_T.t,
                  "NS",
                  y = 0,
                  yend = 2.5,
                  dis = 10,
                  ytext = 2.7,
                  color = "#00FFFF",
                  linetype = "dashed",
                  fontface = fontface
                )
            }
            if (length(which(is.na(df2_syR_T$p_value_Manova))) != 0) {
              df2_syR_T.t <- df2_syR_T %>%
                filter(is.na(get('p_value_Manova')))
              p6 <-
                my_segment_NS(
                  p6,
                  data = df2_syR_T.t,
                  "NS",
                  y = 0,
                  yend = 2.5,
                  dis = 10,
                  ytext = 2.7,
                  color = "#00FFFF",
                  linetype = "dashed",
                  fontface = fontface
                )
            }
          }
        }
        ####################pausing site reverse strand##############
        #Select rows with pausing site event and draw a line
        if (nrow(df2 %>%
                 filter(get('pausing_site') == "+")) != 0) {
          df2_ps <- df2 %>%
            filter(get('pausing_site') == "+")
          #select pausing site duration
          df2_ps <-
            df2_ps %>%
            filter(na.omit(get('event_duration')) <= event_duration)
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
                color = ps_color,
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
                color = ps_color,
                linetype = "dotted",
                df = "pausing",
                fontface = fontface
              )
          }
          if (nrow(df2_ps %>%
                   filter(is.na(
                     get('event_ps_itss_p_value_Ttest')
                   ))) != 0) {
            df2_ps_b <-
              df2_ps %>% filter(is.na(get(
                'event_ps_itss_p_value_Ttest'
              )))
            p6 <-
              my_segment_T(
                p6,
                data = df2_ps_b,
                "PS",
                y = 0,
                yend = 3,
                dis = 50,
                ytext = 3.8,
                color = ps_color,
                linetype = "dotted",
                df = "pausing",
                fontface = fontface
              )
          }
        }
        ####################iTSS_I reverse strand###################
        #Select rows with internal starting site event and draw a line
        if (nrow(df2 %>%
                 filter(get('iTSS_I') == "+")) != 0) {
          #select iTSS
          df2_itss <- df2 %>%
            filter(get('iTSS_I') == "+")
          #select iTSS duration
          df2_itss <-
            df2_itss %>%
            filter(na.omit(get('event_duration')) <= event_duration)
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
                color = iTSS_I_color,
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
                color = iTSS_I_color,
                linetype = "dotted",
                fontface = fontface
              )
          }
          if (nrow(df2_itss %>%
                   filter(is.na(
                     get('event_ps_itss_p_value_Ttest')
                   ))) != 0) {
            df2_itss_b <-
              df2_itss %>%
              filter(is.na(get(
                'event_ps_itss_p_value_Ttest'
              )))
            p6 <-
              my_segment_NS(
                p6,
                data = df2_itss_b,
                "iTSS",
                y = 0,
                yend = 3.2,
                dis = 10,
                ytext = 3.6,
                color = iTSS_I_color,
                linetype = "dotted",
                fontface = fontface
              )
          }
        }
        ####################FC reverse strand###################
        #add FC for intensity ratio test if p_value is significant
        df2 <- indice_function(df2, "intensity_fragment")
        df2_wo <- df2[which(df2$indice == 1),]
        #select the last row for each segment and add 40 nucleotides in case
        #of negative strand for a nice plot
        df2_wo_pvalue <-
          df2_wo[!duplicated(df2_wo$p_value_intensity),]
        if (nrow(df2_wo_pvalue) != 0) {
          df2_p_val_int <-
            df2_wo_pvalue[which(df2_wo_pvalue$p_value_intensity
                                < p_value_int),]
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
                                > p_value_int),]
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
        #add FC for HL ratio test if p_value is significant
        df2 <- indice_function(df2, "HL_fragment")
        df2_wo <- df2[which(df2$indice == 1),]
        df2_wo_pvalue <-  df2_wo[!duplicated(df2_wo$p_value_HL),]
        if (nrow(df2_wo_pvalue) != 0) {
          df2_p_val_hl <-
            df2_wo_pvalue[which(df2_wo_pvalue$p_value_HL < p_value_hl),]
          if (nrow(df2_p_val_hl) != 0) {
            p5 <- p5 +
              geom_text(
                data = df2_p_val_hl,
                aes(
                  x = get('position'),
                  y = get('HL_mean_fragment')
                ),
                label = "FC*",
                size = 1,
                check_overlap = TRUE
              )
          }
          df2_p_val_hl <-
            df2_wo_pvalue[which(df2_wo_pvalue$p_value_HL > p_value_hl),]
          if (nrow(df2_p_val_hl) != 0) {
            p5 <- p5 +
              geom_text(
                data = df2_p_val_hl,
                aes(
                  x = get('position'),
                  y = get('HL_mean_fragment')
                ),
                label = "FC",
                size = 1,
                check_overlap = TRUE
              )
          }
        }
        #Select rows with FC_HL event and draw a line
        df2_hl <- df2[!duplicated(df2$FC_HL),]
        df2_hl <- df2_hl[which(df2_hl$FC_HL < HL_threshold),]
        if (nrow(df2_hl) != 0) {
          p5 <-
            my_segment_NS(
              p5,
              data = df2_hl,
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
        #Select rows with velocity_ratio event and draw a line
        df2_v <- df2[!duplicated(df2$velocity_ratio),]
        df2_v <-
          df2_v[which(df2_v$velocity_ratio < vel_threshold),]
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
        if (nrow(df2) != 0) {
          Title <-
            paste0(
              "Term: termination (",
              forggtitle_ter,
              "), NS: new start (",
              forggtitle_NS,
              "), PS: pausing site (",
              length(which(df2$pausing_site == "+")),
              "), iTSS_I: internal starting site (",
              length(which(df2$iTSS_I == "+")),
              ")"
            )
        } else{
          Title <- NA
        }
      }
      p <- list(p1, p2, p3, p7, p6, p5, p4)
      egg::ggarrange(
        plots = p,
        ncol = 1,
        nrow = 7,
        heights = c(4.5, 4.5, 4.5, 6, 4.5, 4.5, 4.5),
        bottom = textGrob(Title, gp = gpar(fontsize = 6))
      )
    })
    dev.off()
  }
