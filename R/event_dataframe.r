#' event_dataframe: creates a dataframe only with events for segments and genes.
#' The function used are:
#' position_function: adds the specific position of ps or iTSS event
#' annotation_function_event: adds the events to the annotated genes.
#' gff3 file has to be supplied. Strand is indicated in case of stranded data
#' The event_dataframe selects columns with statistical features. ID, position,
#' strand and TU columns are required.
#' Two major dataframe are generated, df gathers t-test and Manova test and df1
#' gathers ps and ITSS with the corresponding features. df selects only
#' unique intensity fragments since they are the lowest on the hierarchy.
#' One new column is added to df, "synthesis_ratio_event", it corresponds
#' to the FC-HL/FC-intensity and assignment of an event to the synthesis ratio
#' respectively. df adds a new column to indicate the position of the ps or
#' ITSS event.
#' @param data dataframe: the probe based data frame.
#' @param data_annotation dataframe: the coordinates are extracted from gff3
#' 
#' @return WIP
#' 
#' @examples
#' data(stats_minimal)
#' event_dataframe(data = stats_minimal,
#' data_annotation = metadata(stats_minimal)$annot[[1]])
#' 
#' @export
#' 
event_dataframe <- function(data, data_annotation) {
  data <- as.data.frame(rowRanges(data))
  data <- data[,-c(1:4)]
  events_df <-
    data[, c(
      "ID",
      "position",
      "strand",
      "delay",
      "delay_fragment",
      "half_life",
      "intensity",
      "TU",
      "pausing_site",
      "iTSS_I",
      "ps_ts_fragment",
      "event_ps_itss_p_value_Ttest",
      "p_value_HL",
      "FC_fragment_HL",
      "p_value_intensity",
      "FC_fragment_intensity",
      "FC_HL",
      "FC_intensity",
      "FC_HL_intensity",
      "p_value_Manova",
      "FC_HL_intensity_fragment",
      "synthesis_ratio",
      "synthesis_ratio_event",
      "FC_HL_adapted",
      "event_duration",
      "delay_frg_slope",
      "p_value_slope"
    )]
  # generate the first dataframe
  events_df1 <-
    events_df[, c(
      "ID",
      "position",
      "strand",
      "TU",
      "delay_fragment",
      "pausing_site",
      "iTSS_I",
      "event_ps_itss_p_value_Ttest",
      "ps_ts_fragment",
      "event_duration",
      "delay",
      "half_life",
      "intensity"
    )]
  # add the specific position to each event
  events_df1[, "event_position"] <- NA
  events_df1 <-
    position_function("pausing_site", events_df1, "position", "event_position")
  events_df1 <-
    position_function("iTSS_I", events_df1, "position", "event_position")
  # add the specific position to each event
  events_df <- full_join(events_df, events_df1, all = TRUE)
  # add gene annotation to the dataframe events
  events_df$region <- NA
  events_df$gene <- NA
  events_df$locus_tag <- NA
  events_df[which(events_df[, "strand"] == "+"), "region"] <-
    annotation_function_event("region", events_df, "+", data_annotation)
  events_df[which(events_df[, "strand"] == "+"), "gene"] <-
    annotation_function_event("gene", events_df, "+", data_annotation)
  events_df[which(events_df[, "strand"] == "+"), "locus_tag"] <-
    annotation_function_event("locus_tag", events_df, "+", data_annotation)
  events_df[which(events_df[, "strand"] == "-"), "region"] <-
    annotation_function_event("region", events_df, "-", data_annotation)
  events_df[which(events_df[, "strand"] == "-"), "gene"] <-
    annotation_function_event("gene", events_df, "-", data_annotation)
  events_df[which(events_df[, "strand"] == "-"), "locus_tag"] <-
    annotation_function_event("locus_tag", events_df, "-", data_annotation)
  events_df <-
    events_df[, c(
      "region",
      "gene",
      "locus_tag",
      "strand",
      "TU",
      "position",
      "FC_fragment_intensity",
      "FC_intensity",
      "p_value_intensity",
      "FC_fragment_HL",
      "FC_HL",
      "p_value_HL",
      "FC_HL_intensity_fragment",
      "FC_HL_intensity",
      "FC_HL_adapted",
      "p_value_Manova",
      "synthesis_ratio",
      "synthesis_ratio_event",
      "pausing_site",
      "iTSS_I",
      "event_ps_itss_p_value_Ttest",
      "ps_ts_fragment",
      "event_position",
      "event_duration",
      "delay_frg_slope",
      "p_value_slope",
      "delay",
      "half_life",
      "intensity"
    )]
  return(events_df)
}
