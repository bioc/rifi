#' dataframe_summary_events creates one table with all events between the
#' segments.
#' The dataframe_summary_events creates one table with the following columns:
#' event, features, p_value, event_position, event_duration, position, region,
#' gene, locus_tag, strand, TU, segment_1, segment_2, length, velocity_ratio,
#' FC_HL, FC_intensity, FC_HL/FC_intensity.
#' The columns are:
#' 1. event: event type, pausing site, iTSS_I, iTSS_II, Termination, HL_event,
#' Int_event, HL_Int_event and velocity_change.
#' 2. FC_HL: fold change between 2 half-life fragments.
#' 3. FC_intensity: fold change between 2 intensity fragments.
#' 4. FC_HL/FC_intensity: ratio of fold change between 2 half-life fragments
#' and fold change between 2 intensity fragments.
#' 5. velocity_ratio: ratio between any two fragment where the event happen.
#' 6. p_value: depending on the event, t-test, manova test p_value is assigned.
#' 7. feature_type: indicated on the output data frame as region, are the
#' feature type covering the event.
#' 8. gene: gene covering the event.
#' 9. locus_tag: locus_tag covering the event.
#' 10. strand: +/- indicated in case of stranded data.
#' 11. TU: TU covering the event.
#' 12. segment_1: the first segment of the event, includes the segment, TU,
#' delay fragment in case of ps or iTSS_I. The rest of the events include HL
#' fragment and intensity fragment.
#' 13. segment_2: same description as segment_1 but is the second fragment
#' of the event.
#' 14. event_position: the position of event, calculated dividing the last
#' position of the first fragment and the first position of the next fragment
#' on 2.
#' 15. event_duration: the difference (min) between 2 delay fragment when
#' ps or iTSS_I happen.
#' 16. gap_fragments: length in position (nt), calculated by the difference
#' between the last position of the first fragment and the first position
#' of the second fragment.
#' 17. features: number of segment involved on the event.
#' @param data SummarizedExperiment: the input data frame with correct format.
#' @param data_annotation dataframe: dataframe from processed gff3 file.
#' 
#' @return 
#' \describe{
#'     \item{event:}{}
#'     \item{FC_HL:}{}
#'     \item{FC_intensity:}{}
#'     \item{FC_HL_FC_intensity:}{}
#'     \item{p_adjusted:}{}
#'     \item{velocity_ratio:}{}
#'     \item{p_value:}{}
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{segment_1:}{}
#'     \item{segment_2:}{}
#'     \item{event_position:}{}
#'     \item{event_duration:}{}
#'     \item{gap_fragments:}{}
#'     \item{features:}{}
#'     }
#'  
#' @examples
#' if(!require(SummarizedExperiment)){
#' suppressPackageStartupMessages(library(SummarizedExperiment))
#' }
#' data(stats_minimal)
#' dataframe_summary_events(data = stats_minimal, 
#' data_annotation = metadata(stats_minimal)$annot[[1]])
#'  
#' @export

dataframe_summary_events <- function(data, data_annotation) {
  tmp_merged <-
    as.data.frame(
      rowRanges(data)[, c(
      "ID",
      "position",
      "position_segment",
      "TU",
      "delay_fragment",
      "HL_fragment",
      "intensity_fragment",
      "velocity_fragment",
      "FC_fragment_HL",
      "FC_HL",
      "p_value_HL",
      "FC_fragment_intensity",
      "FC_intensity",
      "p_value_intensity",
      "FC_HL_adapted",
      "FC_HL_intensity_fragment",
      "synthesis_ratio",
      "synthesis_ratio_event",
      "p_value_Manova",
      "pausing_site",
      "iTSS_I",
      "event_ps_itss_p_value_Ttest",
      "ps_ts_fragment",
      "event_duration",
      "delay_frg_slope",
      "p_value_slope"
    )]
    )
  if(nrow(tmp_merged) != 0){
  tmp_merged <- tmp_merged[,-c(1:4)]
  tmp_merged <-
    tmp_merged[grep("\\TU_\\d+$", tmp_merged$TU), ]
  df <- data.frame()
  event <- c()
  FC_HL <- c()
  FC_intensity <- c()
  FC_HL_adapted <- c()
  FC_HL_FC_intensity <- c()
  velocity_ratio <- c()
  p_value <- c()
  feature_type <- c()
  gene <- c()
  locus_tag <- c()
  strand <- c()
  TU <- c()
  segment_1 <- c()
  segment_2 <- c()
  event_position <- c()
  event_duration <- c()
  gap_fragments <- c()
  features <- c()
  #ps and iTSSI
  ps_frg <- which(tmp_merged$pausing_site == "+")
  itss1_frg <- which(tmp_merged$iTSS_I == "+")
  ps_its <- c(ps_frg, itss1_frg)
  if(length(ps_its) != 0){
    for (i in seq_along(ps_its)) {
    d <- tmp_merged[ps_its[i], ]
    d[which(d$velocity_fragment == Inf), "velocity_fragment"] <- NA
    if (d$pausing_site == "+") {
      event <- c(event, "ps")
    } else if (d$iTSS_I == "+") {
      event <- c(event, "iTSS_I")
    }
    ev_fragments <- unlist(strsplit(d$ps_ts_fragment, split = ":"))
    if (unique(as.character(d$strand) == "-")) {
      ev_fragments <- rev(ev_fragments)
    }
    FC_HL <- c(FC_HL, NA)
    FC_intensity <- c(FC_intensity, NA)
    FC_HL_adapted <- c(FC_HL_adapted, NA)
    FC_HL_FC_intensity <- c(FC_HL_FC_intensity, NA)
    event_position <- c(event_position,
                        (tmp_merged[last(which(tmp_merged$delay_fragment ==
                                                 ev_fragments[1])), "position"]
                         + tmp_merged[which(tmp_merged$delay_fragment ==
                                              ev_fragments[2]), "position"][1])
                        / 2)
    velocity_ratio <- c(velocity_ratio, NA)
    p_value <- c(p_value, as.numeric(d$event_ps_itss_p_value_Ttest))
    feature_type <-
      c(
        feature_type,
        annotation_function_df(
          feature = "region",
          pos = last(event_position),
          strand = d$strand,
          data_annotation = data_annotation
        )
      )
    gene <-
      c(
        gene,
        annotation_function_df(
          feature = "gene",
          pos = last(event_position),
          strand = d$strand,
          data_annotation = data_annotation
        )
      )
    locus_tag <-
      c(
        locus_tag,
        annotation_function_df(
          feature = "locus_tag",
          pos = last(event_position),
          strand = d$strand,
          data_annotation = data_annotation
        )
      )
    strand <- c(strand, as.character(d$strand))
    TU <- c(TU, d$TU)
    segment_1 <-
      c(segment_1, paste(c(
        d$position_segment, d$TU, ev_fragments[1]
      ), collapse = "|"))
    segment_2 <-
      c(segment_2, paste0(c(
        d$position_segment,  d$TU, ev_fragments[2]
      ), collapse = "|"))
    event_duration <- c(event_duration, d$event_duration)
    gap_fragments <- c(gap_fragments,
                       abs(tmp_merged[last(which(tmp_merged$delay_fragment ==
                                                   ev_fragments[1])),
                                      "position"] -
                             tmp_merged[which(tmp_merged$delay_fragment ==
                                                ev_fragments[2]),
                                        "position"][1]))
    features <- c(features, length(unique(ev_fragments)))
  }
  }
  #termination and iTSSII
  tmp <-
    tmp_merged[!duplicated(tmp_merged$FC_HL_intensity_fragment), ]
  ter_frg <- which(tmp$synthesis_ratio_event == "Termination")
  itss2_frg <- which(tmp$synthesis_ratio_event == "iTSS_II")
  ter_its <- c(ter_frg, itss2_frg)
  if(length(ter_its) != 0){
    for (i in seq_along(ter_its)) {
    d <- tmp[ter_its[i], ]
    d[which(d$velocity_fragment == Inf), "velocity_fragment"] <- NA
    event <- c(event, d$synthesis_ratio_event)
    ev_fragments <-
      unlist(strsplit(d$FC_HL_intensity_fragment, split = ";"))
    ev_fragments <- unlist(strsplit(ev_fragments, split = ":"))
    FC_HL <- c(FC_HL, NA)
    FC_intensity <- c(FC_intensity, NA)
    FC_HL_adapted <- c(FC_HL_adapted, d$FC_HL_adapted)
    FC_HL_FC_intensity <- c(FC_HL_FC_intensity, d$synthesis_ratio)
    event_position <-
      c(event_position,
        (tmp_merged[last(which(
          tmp_merged$intensity_fragment == ev_fragments[3])), "position"] +
           tmp_merged[which(
             tmp_merged$intensity_fragment == ev_fragments[4]),
             "position"][1]) / 2)
    velocity_ratio <- c(velocity_ratio, NA)
    p_value <- c(p_value, as.numeric(d$p_value_Manova))
    feature_type <-
      c(
        feature_type,
        annotation_function_df(
          feature = "region",
          pos = last(event_position),
          strand = d$strand,
          data_annotation = data_annotation
        )
      )
    gene <-
      c(
        gene,
        annotation_function_df(
          feature = "gene",
          pos = last(event_position),
          strand = d$strand,
          data_annotation = data_annotation
        )
      )
    locus_tag <-
      c(
        locus_tag,
        annotation_function_df(
          feature = "locus_tag",
          pos = last(event_position),
          strand = d$strand,
          data_annotation = data_annotation
        )
      )
    strand <- c(strand, as.character(d$strand))
    TU <- c(TU, d$TU)
    segment_1 <-
      c(segment_1, paste(
        c(
          d$position_segment,
          d$TU,
          d$delay_fragment,
          ev_fragments[1],
          ev_fragments[3]
        ),
        collapse = "|"
      ))
    segment_2 <-
      c(segment_2, paste0(
        c(
          d$position_segment,
          d$TU,
          d$delay_fragment,
          ev_fragments[2],
          ev_fragments[4]
        ),
        collapse = "|"
      ))
    event_duration <- c(event_duration, NA)
    gap_fragments <-
      c(gap_fragments, abs(tmp_merged[last(which(
        tmp_merged$intensity_fragment == ev_fragments[3])), "position"] -
                             tmp_merged[which(tmp_merged$intensity_fragment ==
                                 ev_fragments[4]), "position"][1]))
    features <- c(features, length(unique(ev_fragments)))
  }
  }
  #FC intensity fragments
  tmp <- tmp_merged[!duplicated(tmp_merged$FC_fragment_intensity), ]
  tmp <- tmp[!is.na(tmp$FC_fragment_intensity), ]
  if(nrow(tmp) != 0){
    for (i in seq_len(nrow(tmp))) {
    d <- tmp[i, ]
    d[which(d$velocity_fragment == Inf), "velocity_fragment"] <- NA
    ev_fragments <-
      unlist(strsplit(d$FC_fragment_intensity, split = ":"))
    event <- c(event, "Int_event")
    FC_HL <- c(FC_HL, NA)
    FC_intensity <- c(FC_intensity, d$FC_intensity)
    FC_HL_adapted <- c(FC_HL_adapted, d$FC_HL_adapted)
    FC_HL_FC_intensity <- c(FC_HL_FC_intensity, NA)
    event_position <-
      c(event_position,
        (tmp_merged[last(which(tmp_merged$intensity_fragment ==
                                 ev_fragments[1])), "position"] +
           tmp_merged[which(tmp_merged$intensity_fragment ==
                              ev_fragments[2]), "position"][1]) / 2)
    velocity_ratio <- c(velocity_ratio, NA)
    p_value <- c(p_value, as.numeric(d$p_value_intensity))
    feature_type <-
      c(
        feature_type,
        annotation_function_df(
          feature = "region",
          pos = last(event_position),
          strand = d$strand,
          data_annotation = data_annotation
        )
      )
    gene <-
      c(
        gene,
        annotation_function_df(
          feature = "gene",
          pos = last(event_position),
          strand = d$strand,
          data_annotation = data_annotation
        )
      )
    locus_tag <-
      c(
        locus_tag,
        annotation_function_df(
          feature = "locus_tag",
          pos = last(event_position),
          strand = d$strand,
          data_annotation = data_annotation
        )
      )
    strand <- c(strand, as.character(d$strand))
    TU <- c(TU, d$TU)
    segment_1 <-
      c(segment_1, paste(
        c(
          d$position_segment,
          d$TU,
          d$delay_fragment,
          d$HL_fragment,
          ev_fragments[1]
        ),
        collapse = "|"
      ))
    segment_2 <-
      c(segment_2, paste0(
        c(
          d$position_segment,
          d$TU,
          d$delay_fragment,
          d$HL_fragment,
          ev_fragments[2]
        ),
        collapse = "|"
      ))
    event_duration <- c(event_duration, NA)
    gap_fragments <-
      c(gap_fragments, abs(tmp_merged[last(
        which(tmp_merged$intensity_fragment ==
                ev_fragments[1])), "position"] -
          tmp_merged[which(tmp_merged$intensity_fragment ==
                             ev_fragments[2]), "position"][1]))
    features <- c(features, length(unique(ev_fragments)))
  }
  }
  #FC HL fragments
  tmp <- tmp_merged[!duplicated(tmp_merged$FC_fragment_HL), ]
  tmp <- tmp[!is.na(tmp$FC_fragment_HL), ]
  if(nrow(tmp) != 0){
    for (i in seq_len(nrow(tmp))) {
    d <- tmp[i, ]
    ev_fragments <- unlist(strsplit(d$FC_fragment_HL, split = ":"))
    d[which(d$velocity_fragment == Inf), "velocity_fragment"] <- NA
    event <- c(event, "HL_event")
    FC_HL <- c(FC_HL, d$FC_HL)
    FC_intensity <- c(FC_intensity, NA)
    FC_HL_adapted <- c(FC_HL_adapted, NA)
    FC_HL_FC_intensity <- c(FC_HL_FC_intensity, NA)
    event_position <-
      c(event_position, (tmp_merged[last(which(
        tmp_merged$HL_fragment == ev_fragments[1])), "position"] +
          tmp_merged[which(tmp_merged$HL_fragment ==
                             ev_fragments[2]), "position"][1]) / 2)
    velocity_ratio <- c(velocity_ratio, NA)
    p_value <- c(p_value, as.numeric(d$p_value_HL))
    feature_type <-
      c(
        feature_type,
        annotation_function_df(
          feature = "region",
          pos = last(event_position),
          strand = d$strand,
          data_annotation = data_annotation
        )
      )
    gene <- c(
      gene,
      annotation_function_df(
        feature = "gene",
        pos = last(event_position),
        strand = d$strand,
        data_annotation = data_annotation
      )
    )
    locus_tag <-
      c(
        locus_tag,
        annotation_function_df(
          feature = "locus_tag",
          pos = last(event_position),
          strand = d$strand,
          data_annotation = data_annotation
        )
      )
    strand <- c(strand, as.character(d$strand))
    TU <- c(TU, d$TU)
    segment_1 <-
      c(segment_1, paste(
        c(d$position_segment, d$TU, d$delay_fragment, ev_fragments[1]),
        collapse = "|"
      ))
    segment_2 <-
      c(segment_2, paste0(
        c(d$position_segment, d$TU, d$delay_fragment, ev_fragments[2]),
        collapse = "|"
      ))
    event_duration <- c(event_duration, NA)
    gap_fragments <-
      c(gap_fragments, abs(tmp_merged[last(which(tmp_merged$HL_fragment ==
                                             ev_fragments[1])), "position"] -
                             tmp_merged[which(tmp_merged$HL_fragment ==
                                             ev_fragments[2]), "position"][1]))
    features <- c(features, length(unique(ev_fragments)))
  }
  }
  #ratio velocity
  uniqDelay <- unique(na.omit(tmp_merged$delay_frg_slope))
  tmp <- tmp_merged[!duplicated(tmp_merged$delay_frg_slope), ]
  if(nrow(tmp) != 0){
    for (i in seq_along(uniqDelay)) {
    ev_fragments <- unlist(strsplit(uniqDelay[i], split = ":"))
    d <- tmp[which(tmp$delay_frg_slope == uniqDelay[i]), ]
    d[which(d$velocity_fragment == Inf), "velocity_fragment"] <- NA
    event <- c(event, "velocity")
    FC_HL <- c(FC_HL, NA)
    FC_intensity <- c(FC_intensity, NA)
    FC_HL_adapted <- c(FC_HL_adapted, NA)
    FC_HL_FC_intensity <- c(FC_HL_FC_intensity, NA)
    event_position <-
      c(event_position, (tmp_merged[last(which(
        tmp_merged$delay_fragment == ev_fragments[1])), "position"] +
        tmp_merged[which(tmp_merged$delay_fragment ==
                           ev_fragments[2]), "position"][1]) / 2)
    velocity_ratio <-
      c(velocity_ratio, tmp_merged[which(
        tmp_merged$delay_fragment == ev_fragments[2])[1], "velocity_fragment"] /
          tmp_merged[which(tmp_merged$delay_fragment ==
                             ev_fragments[1])[1], "velocity_fragment"])
    p_value <- c(p_value, as.numeric(d$p_value_slope))
    feature_type <-
      c(
        feature_type,
        annotation_function_df(
          feature = "region",
          pos = last(event_position),
          strand = d$strand,
          data_annotation = data_annotation
        )
      )
    gene <-
      c(
        gene,
        annotation_function_df(
          feature = "gene",
          pos = last(event_position),
          strand = d$strand,
          data_annotation = data_annotation
        )
      )
    locus_tag <-
      c(
        locus_tag,
        annotation_function_df(
          feature = "locus_tag",
          pos = last(event_position),
          strand = d$strand,
          data_annotation = data_annotation
        )
      )
    strand <- c(strand, as.character(d$strand))
    TU <- c(TU, d$TU)
    segment_1 <-
      c(segment_1, paste(c(
        d$position_segment, d$TU, ev_fragments[1]
      ), collapse = "|"))
    segment_2 <-
      c(segment_2, paste0(c(
        d$position_segment, d$TU, ev_fragments[2]
      ), collapse = "|"))
    event_duration <- c(event_duration, NA)
    gap_fragments <-
      c(gap_fragments, abs(tmp_merged[last(which(
        tmp_merged$delay_fragment == ev_fragments[1])), "position"] -
          tmp_merged[which(tmp_merged$delay_fragment ==
                             ev_fragments[2]), "position"][1]))
    features <- c(features, length(unique(ev_fragments)))
    }
  }
  df <-
    cbind.data.frame(
      event,
      p_value,
      FC_HL,
      FC_intensity,
      FC_HL_adapted,
      FC_HL_FC_intensity,
      event_position,
      velocity_ratio,
      feature_type,
      gene,
      locus_tag,
      strand,
      TU,
      segment_1,
      segment_2,
      event_duration,
      gap_fragments,
      features
    )
  if(nrow(df) != 0){
  df$p_value <- formatC(as.numeric(as.character(df$p_value)), format = "e", 
                        digits = 2)
  df <-
    as.data.frame(df %>% mutate_if(is.numeric, round, digits = 2))
  p_adjusted <- as.numeric(as.character(p.adjust(df$p_value, method = "fdr")))
  df <-
    tibble::add_column(df, formatC(p_adjusted, format = "e", digits = 2),
                       .after = 2)
  colnames(df)[3] <- "p_adjusted"
    }
  }
  return(df)
}
