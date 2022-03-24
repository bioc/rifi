#' dataframe_summary_TI creates one table with all TI fragments, p_value and
#' the coordinates.
#' The dataframe_summary creates one table with the following columns: event,
#' TI_fragment, TI_factor, TI_fragments_TU, p_value, feature_type,
#' gene, locus_tag, strand, TU, features, event_position, position_1 and
#' position_2.
#' The columns are:
#' 1. event: event type, transcription interference.
#' 2. TI_fragment: Transcription interference fragment.
#' 3. TI_factor: Transcription interference factor.
#' 4. TI_fragments_TU: Transcription interference fragments included on the
#' TU.
#' 5. p_value: TI p_value between two successive fragments is assigned.
#' 6. feature_type: indicated on the output data frame as region, are the
#' feature type covering the TI.
#' 7. gene: the genes covering the TI.
#' 8. locus_tag: the locus_tags covering the TI.
#' 9. strand: +/- indicated in case of stranded data.
#' 10. TU: TU covering the TI.
#' 11. features: number of segment TI involved on a TU.
#' 12. event_position : position between two TI fragments.
#' 13. position_1 : the first position of TI fragment, if 2 fragments, first
#' position is from the first fragment.
#' 14. position_2 : the last position of TI fragment, if 2 fragments, last
#' position is from the second fragment.
#'
#' @param data SummarizedExperiment: the input data frame with correct format.
#' @param input dataframe: dataframe from event_dataframe function.
#' 
#' @return WIP
#' 
#' @examples
#' data(stats_minimal)
#' data(res_minimal)
#' dataframe_summary_TI(data = stats_minimal, input = res_minimal)
#' 
#' @export

dataframe_summary_TI <- function(data, input) {
  tmp <-
    as.data.frame(
    rowRanges(data)[, c(
      "ID",
      "position",
      "position_segment",
      "flag",
      "TU",
      "delay_fragment",
      "HL_fragment",
      "intensity_fragment",
      "velocity_fragment",
      "event_duration",
      "delay_frg_slope",
      "p_value_slope",
      "TI_termination_fragment",
      "TI_mean_termination_factor",
      "p_value_TI",
      "TI_fragments_p_value"
    )]
    )
  tmp <- tmp[,-c(1:4)] 
  tmp_event <-
    input[, c(
      "region",
      "gene",
      "locus_tag",
      "FC_fragment_HL",
      "FC_HL",
      "p_value_HL",
      "FC_fragment_intensity",
      "FC_intensity",
      "p_value_intensity",
      "FC_HL_intensity",
      "FC_HL_intensity_fragment",
      "synthesis_ratio",
      "synthesis_ratio_event",
      "p_value_Manova",
      "pausing_site",
      "iTSS_I",
      "event_ps_itss_p_value_Ttest",
      "ps_ts_fragment",
      "event_position"
    )]
  tmp_merged <- cbind(tmp, tmp_event)
  tmp_merged <-
    tmp_merged[grep("\\TU_\\d+$", tmp_merged$TU), ]
  tmp <- tmp_merged[grep("TI", tmp_merged$flag), ]
  uniqTU <- unique(na.omit(tmp$TU))
  uniqTU <- uniqTU[grep("_T|_O|_NA", uniqTU, invert = TRUE)]
  df <- data.frame()
  event <- c()
  TI_fragment <- c()
  TI_termination_factor <- c()
  p_value <- c()
  feature_type <- c()
  gene <- c()
  locus_tag <- c()
  strand <- c()
  TU <- c()
  features <- c()
  event_position <- c()
  position_1 <- c()
  position_2 <- c()
  for (i in seq_along(uniqTU)) {
    d <- tmp[which(tmp$TU %in% uniqTU[i]), ]
    d <-
      d[grep("_T|_O|_NA", d$TI_termination_fragment, invert = TRUE), ]
    d[which(d$velocity_fragment == Inf), "velocity_fragment"] <- NA
    ev_fragments <- unique(na.omit(d$TI_termination_fragment))
    if (is_empty(ev_fragments)) {
      next ()
    }
    TI_frg <- unique(na.omit(d$TI_fragments_p_value))
    if (!is_empty(TI_frg)) {
      for (j in seq_along(TI_frg)) {
        ti_frg <- unique(unlist(strsplit(TI_frg[j], split = ":")))
        event <- c(event, "TI")
        TI_fragment <- c(TI_fragment, paste(ti_frg, collapse = ":"))
        TI_termination_factor <-
          c(TI_termination_factor, paste(c(
            round(unique(tmp[which(
              tmp$TI_termination_fragment == ti_frg[1]),
              "TI_mean_termination_factor"]), digits = 2),
            round(unique(tmp[which(tmp$TI_termination_fragment == ti_frg[2]),
              "TI_mean_termination_factor"]), digits = 2)
          ),
          collapse = "|"))
        p_value <-
          c(p_value, formatC(unique(tmp[which(
            tmp$TI_termination_fragment == ti_frg[1]), "p_value_TI"]),
                             format = "e", digits = 2))
        feature_type <-
          c(feature_type, paste(c(unique(tmp[which(
            tmp$TI_termination_fragment == ti_frg[1]), "region"]),
            (unique(tmp[which(tmp$TI_termination_fragment == ti_frg[2]),
                        "region"]))), collapse = "|"))
        gene <-
          c(gene, paste(c(unique(tmp[which(
            tmp$TI_termination_fragment == ti_frg[1]), "gene"]),
                          (unique(tmp[which(
                            tmp$TI_termination_fragment ==
                              ti_frg[2]), "gene"]))), collapse = "|"))
        locus_tag <-
          c(locus_tag, paste(c(unique(tmp[which(
            tmp$TI_termination_fragment == ti_frg[1]), "locus_tag"]),
            (unique(tmp[which(tmp$TI_termination_fragment ==
                ti_frg[2]), "locus_tag"]))), collapse = "|"))
        strand <- c(strand, as.character(unique(d$strand)))
        TU <- c(TU, unique(d$TU))
        features <- c(features, length(ev_fragments))
        event_position <-
          c(event_position, (tmp[last(which(
            tmp$TI_termination_fragment == ti_frg[1])), "position"] +
                               tmp[which(tmp$TI_termination_fragment ==
                                           ti_frg[2]), "position"][1]) / 2)
        position_1 <-
          c(position_1, tmp[which(tmp$TI_termination_fragment ==
                                    ti_frg[1]), "position"][1])
        position_2 <-
          c(position_2, last(tmp[which(tmp$TI_termination_fragment ==
                                         ti_frg[2]), "position"]))
      }
    } else{
      event <- c(event, "TI")
      p_value <- c(p_value, NA)
      TU <- c(TU, unique(d$TU))
      strand <- c(strand, as.character(unique(d$strand)))
      if (length(ev_fragments) == 2) {
        TI_fragment <- c(TI_fragment, paste(ev_fragments, collapse = ":"))
        TI_termination_factor <-
          c(TI_termination_factor, paste(c(
            round(unique(tmp[which(
              tmp$TI_termination_fragment ==
                ev_fragments[1]), "TI_mean_termination_factor"]), digits = 2),
            round(unique(tmp[which(
              tmp$TI_termination_fragment ==
                ev_fragments[2]), "TI_mean_termination_factor"]), digits = 2)
          ),
          collapse = "|"))
        feature_type <-
          c(feature_type, paste(c(unique(tmp[which(
            tmp$TI_termination_fragment == ev_fragments[1]), "region"]),
            (unique(tmp[which(tmp$TI_termination_fragment ==
                ev_fragments[2]), "region"]))), collapse = "|"))
        gene <-
          c(gene, paste(c(unique(tmp[which(
            tmp$TI_termination_fragment == ev_fragments[1]), "gene"]),
            (unique(tmp[which(tmp$TI_termination_fragment ==
                                ev_fragments[2]), "gene"]))), collapse = "|"))
        locus_tag <-
          c(locus_tag, paste(c(unique(tmp[which(
            tmp$TI_termination_fragment == ev_fragments[1]), "locus_tag"]),
            (unique(tmp[which(tmp$TI_termination_fragment ==
                                ev_fragments[2]), "locus_tag"]))),
            collapse = "|"))
        position_1 <-
          c(position_1, tmp[which(tmp$TI_termination_fragment ==
                                    ev_fragments[2]), "position"][1])
        position_2 <-
          c(position_2, last(tmp[which(tmp$TI_termination_fragment ==
                                         ev_fragments[2]), "position"]))
        features <- c(features, 2)
        event_position <-
          c(event_position, (tmp[last(which(tmp$TI_termination_fragment ==
                                              ev_fragments[1])), "position"] +
                               tmp[which(tmp$TI_termination_fragment ==
                                           ev_fragments[2]), "position"][1]) /
              2)
      } else{
        TI_fragment <-
          c(TI_fragment, unique(na.omit(d$TI_termination_fragment)))
        TI_termination_factor <-
          c(TI_termination_factor, round(unique(
            na.omit(d$TI_mean_termination_factor)
          ), digits = 2))
        feature_type <-
          c(feature_type, paste(unique(unlist(
            strsplit(d$region, split = ";")
          )), collapse = "|"))
        gene <-
          c(gene, paste(unique(unlist(
            strsplit(d$gene, split = ";")
          )), collapse = "|"))
        locus_tag <-
          c(locus_tag, paste(unique(unlist(
            strsplit(d$locus_tag, split = ";")
          )), collapse = "|"))
        position_1 <-
          c(position_1, tmp[which(tmp$TI_termination_fragment ==
                                    ev_fragments[1]), "position"][1])
        position_2 <-
          c(position_2, last(tmp[which(tmp$TI_termination_fragment ==
                                         ev_fragments[1]), "position"]))
        features <- c(features, 1)
        event_position <- c(event_position, NA)
      }
    }
  }
  df <-
    cbind.data.frame(
      event,
      TI_fragment,
      TI_termination_factor,
      p_value,
      feature_type,
      gene,
      locus_tag,
      strand,
      TU,
      features,
      event_position,
      position_1,
      position_2
    )
  if(nrow(df) != 0){
  p_adjusted <-
      p.adjust(as.numeric(as.character(df$p_value)), method = "fdr")
  df <-
    tibble::add_column(df, formatC(p_adjusted, format = "e", digits = 2),
                       .after = 4)
  colnames(df)[5] <- "p_adjusted"
  }
  return(df)
}
