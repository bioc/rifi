##########################annotation section#########################
annotation_plot <-
  function(data_p,
           data_n,
           annot,
           tmp.1,
           tmp.2,
           frag,
           i,
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
           size_gene = size_gene) {
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
      if(nrow(data_p) != 0){
        segment_data_t.1 <- TU_annotation(data_p, "TU", 3, 6, color_TU)
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
            function_TU_arrow(data_p, tmp.1, frag = frag, i = i)
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
      if(length(an.p) != 0){
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
      if (nrow(data_n) != 0) {
        segment_data_t.2 <- TU_annotation(data_n, "TU", -3, -6, color_TU)
        segment_data_t.2 <-
          segment_data_t.2[grep("_NA", segment_data_t.2$annotation,
                                invert = TRUE),]
        if (nrow(segment_data_t.2) == 0) {
          segment_data_t.2 <- empty_boxes(-3, -6, frag = frag, i = i)
        } else{
          #function to find TUs split into 2 pages
          tu_border <-
            function_TU_arrow(data_n, tmp.2, frag = frag, i = i)
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
      if(length(an.m) != 0){
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
      if(nrow(data_p) == 0){
        segment_data_t.1 <- empty_boxes(3, 6, frag = frag, i = i)
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
      
      if(nrow(data_n) == 0){
      segment_data_t.2 <- empty_boxes(-3, -6, frag = frag, i = i)
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
    
return(p7)
}
