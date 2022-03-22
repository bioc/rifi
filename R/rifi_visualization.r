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
#' annotation_plot: plots the corresponding annotation 
#' positive_strand_function: plots delay, HL, intensity and events of positive 
#' strand
#' negative_strand_function: plots delay, HL, intensity and events of negative
#' strand
#' empty_data_positive: plots empty boxes in case no data is available for 
#' positive strand
#' empty_data_negative: plots empty boxes in case no data is available for
#' negative strand 
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
#' @param data SummarizedExperiment: the input data frame with correct format.
#' @param genomeLength integer: genome length output of gff3_preprocess
#' function and element of metadata of SummarizedExperiment.
#' @param annot dataframe: the annotation file, output of gff3_preprocess
#' function and element of metadata of SummarizedExperiment.
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
#' @param HL_threshold_1 integer: threshold for log2FC(HL) selected to plot, 
#' default log2(1.5). log2FC(HL) >= log2(1.5) are indicated by black color.
#' If p_value <= p_value_hl (default 0.05), log2FC(HL) is indicated by HL* 
#' otherwise HL.
#' @param HL_threshold_2 integer: threshold for log2FC(HL) selected to plot, 
#' default -log2(1.5). log2FC(HL) <= -log2(1.5) are indicated by green color.
#' If p_value <= p_value_hl (default 0.05), log2FC(HL) is indicated by HL* 
#' otherwise HL. 
#' In case of p_value is significant and the log2FC(HL) is between -log2FC(1.5) 
#' and log2FC(1.5), FC is assigned by green color and HL*. 
#' @param HL_threshold_color_1 string: color for HL fold change plot.
#' @param vel_threshold integer: threshold for velocity ratio selected to plot,
#' default 200.
#' @param termination_threshold integer: threshold for termination to plot,
#' default .8.
#' @param vel_threshold_color string: color for velocity ratio plot.
#' @param iTSS_threshold integer: threshold for iTSS_II selected to plot,
#' default 1.2.
#' @param event_duration_ps integer: threshold for pausing sites selected to 
#' plot, default -2.
#' @param event_duration_itss integer: threshold for iTSS_I selected to 
#' plot, default 2.
#' @param ps_color string: color for pausing site plot.
#' @param iTSS_I_color string: color for iTSS_I plot.
#'
#' @return The visualization.
#'
#' @examples
#' data(stats_minimal)
#' rifi_visualization(data = as.data.frame(rowRanges(stats_minimal)), 
#' genomeLength = metadata(stats_minimal)$annot[[2]],
#' annot = metadata(stats_minimal)$annot[[1]], coverage = 0, 
#' chr_fwd = NA, chr_rev = NA,
#' region = c("CDS","asRNA","5'UTR","ncRNA","3'UTR","tRNA"),
#' color_region = c("grey0", "red", "blue", "orange", "yellow", "green",
#' "white", "darkseagreen1", "grey50", "black"),
#' color_text.1 = "grey0", color_text.2 = "black", color_TU = "blue",
#' size_tu = 1.6, size_locusTag = 1.6, size_gene = 1.6, Limit = 10,
#' shape=22, col_outiler = "grey50", Alpha=0.5,
#' col_coverage = "grey", shape_outlier = 13, limit_intensity = NA,
#' face="bold", tick_length = .3, arrow.color = "darkseagreen1",
#' minVelocity = 3000, medianVelocity = 6000, 
#' col_above20 = "#00FFFF", fontface = "plain", shape_above20 = 14,
#' axis_text_y_size = 3, axis_title_y_size = 6, TI_threshold = 1.1,
#' p_value_TI=0.05, p_value_manova = 0.05, termination_threshold = 1,
#' iTSS_threshold = 1.01, p_value_int = 0.05, p_value_event = 0.05,
#' p_value_hl = 0.05, event_duration_ps = -2, event_duration_itss = 2,
#' HL_threshold=20, vel_threshold = 200, HL_threshold_color_1="black",
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
           Alpha = 0.5,
           size_tu = 1.6,
           size_locusTag = 1.6,
           size_gene = 1.6,
           Limit = 10,
           shape = 22,
           col_outiler = "grey50",
           col_coverage = "grey",
           shape_outlier = 13,
           limit_intensity = NA,
           face = "bold",
           tick_length = .3,
           arrow.color = "darkseagreen1",
           minVelocity = 3000,
           medianVelocity = 6000,
           col_above20 = "#00FFFF",
           fontface = "plain",
           shape_above20 = 14,
           axis_text_y_size = 3,
           axis_title_y_size = 6,
           TI_threshold = 1.1,
           termination_threshold = 0.8,
           iTSS_threshold = 1.2,
           p_value_int = 0.05,
           p_value_event = 0.05,
           p_value_hl = 0.05,
           p_value_TI = 0.05,
           p_value_manova = 0.05,
           event_duration_ps = 1,
           event_duration_itss = -1,
           HL_threshold_1 = log2(1.5),
           HL_threshold_2 = -log2(1.5),
           vel_threshold = 200,
           HL_threshold_color_1 = "black",
           HL_threshold_color_2 = "green",
           vel_threshold_color = "grey52",
           ps_color = "orange",
           iTSS_I_color = "blue") {
    ##########################data preparation##################################
    #I. add coverage if its available from RNA-seq
    data <- as.data.frame(rowRanges(data))
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
      #adjust position for genes split on two pages
      pos.1 <- frag[i] - 2000
      pos.2 <- frag[c(i + 1)] + 2000
      ###########################data adjustment###########################
      #define the main dataframe with segments positive strand df1, negative
      #strand df2
      df1 <-
        tmp.1[between(tmp.1$position, frag[i], frag[c(i + 1)]), ]
      df2 <-
        tmp.2[between(tmp.2$position, frag[i], frag[c(i + 1)]), ]
      df1_1 <- df1[!is.na(df1$ID), ]
      #avoid plot empty pages in case of small data
      if (nrow(df1) == 0 & nrow(df2) == 0 & nrow(data) < 10000) {
        next ()
      }
      ##########################annotation section#########################
      #an is the annotation dataframe upon the position on the plot, its used
      # to loop into exactly the number of region contained in the gff3
      an.1 <- annot[between(annot$start, frag[i], frag[c(i + 1)]),]
      an <- annot[between(annot$start, pos.1, pos.2),]
      an <- an[!duplicated(an),]
      an <- as.data.frame(an)
      #in case of no data nor annotation are available
      if (nrow(an.1) == 0 & nrow(df1) == 0 & nrow(df2) == 0) {
        next ()
      }
      p7 <-
        annotation_plot(
          data_p = df1,
          data_n = df2,
          annot = annot,
          tmp.1 = tmp.1,
          tmp.2 = tmp.2,
          frag = frag,
          i = i,
          an = an,
          region = region,
          color_region = color_region,
          fontface = fontface,
          color_text.1 = color_text.1,
          color_text.2 = color_text.2,
          color_TU = color_TU,
          Alpha = Alpha,
          size_tu = size_tu,
          size_locusTag = size_locusTag,
          termination_threshold = 
            termination_threshold,
          iTSS_threshold =
            iTSS_threshold,
          p_value_manova = 
            p_value_manova,
          size_gene = size_gene,
          pos.1 = pos.1,
          pos.2 = pos.2)
      #########################empty data positive strand###################
      if (nrow(df1) == 0) {
        p_positive <- empty_data_positive(data_p = df1, data_n = df2,
                                          frag = frag, i = i,  
                                          axis_title_y_size = axis_title_y_size, 
                                          axis_text_y_size = axis_text_y_size,
                                          Limit = Limit)
        p1 <- p_positive[[1]]
        p2 <- p_positive[[2]]
        p3 <- p_positive[[3]]
      }
      ############################positive_strand_plot#######################
      if (nrow(df1) != 0) {
        p_positive <- positive_strand_function(data_p = df1, data = data, 
                                               tmp.c1 = tmp.c1, df1_1 = df1_1,
                                               frag = frag, i = i,
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
                                               axis_text_y_size = 
                                                 axis_text_y_size,
                                               axis_title_y_size = 
                                                 axis_title_y_size,
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
                                               vel_threshold_color = 
                                                 vel_threshold_color,
                                               ps_color = ps_color,
                                               iTSS_I_color = iTSS_I_color)
        p1 <- p_positive[[1]]
        p2 <- p_positive[[2]]
        p3 <- p_positive[[3]]
      }
      #########################empty data reverse strand###################
      if (nrow(df2) == 0) {
        p_negative <- empty_data_negative(data_n = df2, frag = frag, i = i,
                                          axis_title_y_size = axis_title_y_size, 
                                          axis_text_y_size = axis_text_y_size,
                                          Limit = Limit)
        p6 <- p_negative[[1]]
        p5 <- p_negative[[2]]
        p4 <- p_negative[[3]]
      }
      #############################reverse_strand_plot#######################
      if (nrow(df2) != 0) {
        p_negative <- negative_strand_function(data_n = df2, data = data, 
                                               frag = frag, i = i,
                                               Limit = Limit,
                                               shape = shape,
                                               col_outiler = col_outiler,
                                               col_coverage = col_coverage,
                                               shape_outlier = shape_outlier,
                                               limit_intensity = 
                                                 limit_intensity,
                                               face = face,
                                               tick_length = tick_length,
                                               arrow.color = arrow.color,
                                               minVelocity = minVelocity,
                                               medianVelocity = medianVelocity,
                                               shape_above20 = shape_above20,
                                               col_above20 = col_above20,
                                               fontface = fontface,
                                               coverage = coverage,
                                               axis_text_y_size = 
                                                 axis_text_y_size,
                                               axis_title_y_size = 
                                                 axis_title_y_size,
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
                                               vel_threshold_color = 
                                                 vel_threshold_color,
                                               ps_color = ps_color,
                                               iTSS_I_color = iTSS_I_color)
        p6 <- p_negative[[1]]
        p5 <- p_negative[[2]]
        p4 <- p_negative[[3]]
      }
      ############################Title and plot#############################
      p <- list(p1, p2, p3, p7, p6, p5, p4)
      egg::ggarrange(
        plots = p,
        ncol = 1,
        nrow = 7,
        heights = c(4.5, 4.5, 4.5, 6, 4.5, 4.5, 4.5),
        bottom = textGrob(p_negative[[4]], gp = gpar(fontsize = 6))
      )
    })
    dev.off()
  }
