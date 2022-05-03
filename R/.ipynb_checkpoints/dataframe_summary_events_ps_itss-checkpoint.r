# =========================================================================
# dataframe_summary_events_ps_itss  Creates one table with all events 
#                                    between the segments.
# -------------------------------------------------------------------------
#'
#' 
#' The dataframe_summary_events_ps_itss creates one table with the following
#' columns: event, features, p_value, event_position, event_duration, position,
#' region, gene, locus_tag, strand, TU, segment_1, segment_2, length,
#' velocity_ratio.

#' The columns are:

#' 1. event: event type, pausing site, iTSS_I, iTSS_II, Termination, HL_event,
#' Int_event, HL_Int_event and velocity_change.

#' 2. velocity_ratio: velocity ratio between two delay fragments.

#' 3. p_value: depending on the event, t-test, manova test p_value is assigned.

#' 4. feature_type: indicated on the output data frame as region, are the
#' feature type covering the event.

#' 5. gene: gene covering the event.

#' 6. locus_tag: locus_tag covering the event.

#' 7. strand: +/- indicated in case of stranded data.

#' 8. TU: TU covering the event.

#' 9. segment_1: the first segment of the event, includes the segment, TU,
#' delay fragment in case of ps or iTSS_I.

#' 10. segment_2: same description as segment_1 but is the second fragment
#' of the event.

#' 11. event_position: the position of event, calculated dividing the last
#' position of the first fragment and the first position of the next fragment
#' on 2.

#' 12. event_duration: the difference (min) between 2 delay fragment when ps or
#' iTSS_I happen.

#' 13. gap_fragments: length in position (nt), calculated by the difference
#' between the last position of the first fragment and the first position of
#' the second fragment.

#' 14. features: number of segment involved on the event.
#'
#' @param data SummarizedExperiment: the input data frame with correct format.
#' @param data_annotation dataframe: dataframe from processed gff3 file.
#' 
#' @return
#'   \describe{
#'     \item{event:}{String, event type}
#'     \item{p_value:}{Integer, p_value of the event}
#'     \item{p_adjusted:}{Integer, p_value adjusted}
#'     \item{event_position:}{Integer, the position middle between 2 fragments with an event}
#'     \item{velocity_ratio:}{Integer, the ratio value of velocity from 2 delay fragments}
#'     \item{feature_type:}{String, region annotation covering the fragments}
#'     \item{gene:}{String, gene annotation covering the fragments}
#'     \item{locus_tag:}{String, locus_tag annotation covering the fragments}
#'     \item{strand:}{Boolean. The bin/probe specific strand (+/-)}
#'     \item{TU:}{String, The overarching transcription unit}
#'     \item{segment_1:}{String, the first fragment of the two of fragments subjected to analysis}
#'     \item{segment_2:}{String, the second fragment of the two of fragments subjected to analysis}
#'     \item{event_duration:}{Integer, the duration between two delay fragments}
#'     \item{gap_fragments:}{Integer, the distance between two delay fragments}
#'     \item{features:}{Integer, number of fragements involved on the event}
#'     }
#' 
#'  
#' @examples 
#' data(stats_minimal)
#' if(!require(SummarizedExperiment)){
#' suppressPackageStartupMessages(library(SummarizedExperiment))
#' }
#' dataframe_summary_events_ps_itss(data = stats_minimal,
#' data_annotation = metadata(stats_minimal)$annot[[1]])
#' 
#' @export

dataframe_summary_events_ps_itss <-
  function(data, data_annotation) {
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
    tmp_merged <- tmp_merged[,-c(1:4)]
    tmp_merged <-
      tmp_merged[grep("\\TU_\\d+$", tmp_merged$TU), ]
    df <- data.frame()
    event <- c()
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
    for (i in seq_along(ps_its)) {
      d <- tmp_merged[ps_its[i], ]
      d[which(d$velocity_fragment == Inf), "velocity_fragment"] <- NA
      if (d$pausing_site == "+") {
        event <- c(event, "ps")
      } else{
        event <- c(event, "iTSS_I")
      }
      ev_fragments <- unlist(strsplit(d$ps_ts_fragment, split = ":"))
      if (unique(as.character(d$strand) == "-")) {
        ev_fragments <- rev(ev_fragments)
      }
      event_position <-
        c(event_position, (tmp_merged[last(
          which(tmp_merged$delay_fragment ==
                  ev_fragments[1])), "position"] + tmp_merged[
                    which(tmp_merged$delay_fragment ==
                            ev_fragments[2]), "position"][1]) / 2)
      velocity_ratio <-
        c(velocity_ratio, tmp_merged[
          which(tmp_merged$delay_fragment ==
                  ev_fragments[2])[1], "velocity_fragment"] /
            tmp_merged[which(tmp_merged$delay_fragment ==
                               ev_fragments[1])[1], "velocity_fragment"])
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
      gap_fragments <-
        c(gap_fragments, abs(tmp_merged[last(which(
          tmp_merged$delay_fragment == ev_fragments[1])), "position"] -
            tmp_merged[which(tmp_merged$delay_fragment ==
                               ev_fragments[2]), "position"][1]))
      features <- c(features, length(unique(ev_fragments)))
    }
    df <-
      cbind.data.frame(
        event,
        p_value,
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
      df$p_value <- formatC(as.numeric(df$p_value), format = "e", digits = 2)
      df <-
        as.data.frame(df %>% mutate_if(is.numeric, round, digits = 2))
      p_adjusted <-
        as.numeric(p.adjust(as.numeric(df$p_value), method = "fdr"))
      df <-
        tibble::add_column(df, formatC(p_adjusted, format = "e", digits = 2),
                           .after = 2)
      colnames(df)[3] <- "p_adjusted"
    }
    return(df)
  }
