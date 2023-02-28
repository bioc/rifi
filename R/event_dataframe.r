#' =========================================================================
#' event_dataframe                                    
#' -------------------------------------------------------------------------
#' event_dataframe creates a dataframe only with events 
#'
#' event_dataframe creates a dataframe connecting segments, events and 
#' the annotation.
#' 
#' The functions used are:
#' 
#' position_function: adds the specific position of ps or iTSS event.
#' 
#' annotation_function_event: adds the events to the annotated genes.
#' 
#' annotation file needs to be supplied. Strand is indicated in case of 
#' stranded data
#' The event_dataframe selects columns with statistical features. 
#' ID, position, strand and TU columns are required.
#' 
#' @param data dataframe: the probe based data frame.
#' @param data_annotation dataframe: the coordinates are extracted from gff3
#' 
#' @return A dataframe with unique intensity fragments 
#'   \describe{
#'     \item{feature_type:}{String, region annotation covering the fragments}
#'     \item{gene:}{String, gene annotation covering the fragments}
#'     \item{locus_tag:}{String, locus_tag annotation covering the fragments}
#'     \item{strand:}{Boolean. The bin/probe specific strand (+/-)}
#'     \item{TU:}{String, The overarching transcription unit}
#'     \item{position:}{Integer, position of the bin/probe on the genome}
#'     \item{segment:}{String, the bin/probe segment on the genome}
#'     \item{FC_fragment_intensity:}{String, the fragments subjected to
#'           fold change}
#'     \item{FC_intensity:}{Integer, the fold change value of 2 intensity 
#'     fragments}
#'     \item{p_value_intensity:}{Integer, p_value of the FC_intensity}
#'     \item{FC_fragment_HL:}{String, the fragments subjected to fold change}
#'     \item{FC_HL:}{Integer, the fold change value of 2 HL fragments}
#'     \item{p_value_HL:}{Integer, p_value of the FC_HL}
#'     \item{FC_HL_FC_intensity_fragment:}{String, fragments subjected to 
#'          FC_HL/FC_intensity}
#'     \item{FC_HL_FC_intensity:}{Integer, the value of FC_HL/FC_intensity}
#'     \item{FC_HL_adapted:}{Integer, the fold change of half-life/ fold change
#'      of intensity,
#'     position of the half-life fragment is adapted to intensity fragment}
#'     \item{p_value_manova:}{Integer, p_value of the event FC_HL/FC_intensity}
#'     \item{synthesis_ratio:}{Integer, the value correspomding to synthesis 
#'     rate}
#'     \item{synthesis_ratio_event:}{String, the event assigned by synthesis 
#'     rate either 
#'       Termination or iTSS}
#'     \item{pausing_site:}{Boolean, presence or absence of pausing_site event 
#'     (ps)}
#'     \item{iTSS_I:}{Boolean, presence or absence of internal starting site 
#'     event (iTSS_I)}
#'     \item{ps_ts_fragment:}{String, fragments involved on the event ps or 
#'     iTSS_I}
#'     \item{event_position:}{Integer, the position middle between 2 fragments 
#'     with an event}
#'     \item{event_duration:}{Integer, the duration between two delay fragments}
#'     \item{delay:}{Integer, the delay value of the bin/probe}
#'     \item{half_life:}{The half-life of the bin/probe}
#'     \item{intensity:}{The relative intensity at time point 0}
#'     \item{delay_frg_slope:}{Integer, the slope value of the fit through the
#'      respective delay fragment}
#'     \item{p_value_slope:}{Integer, the p_value added to the inp}
#'     }
#' 
#' @examples
#' if(!require(SummarizedExperiment)){
#' suppressPackageStartupMessages(library(SummarizedExperiment))
#' }
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
  events_df <- full_join(events_df, events_df1)
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
