#' dataframe_summary: creates two tables relating gene annotation to fragments.
#' dataframe_summary creates two tables summary of segments and their
#' half-lives. The first output is bin/probe features and the second one is
#' intensity fragment based.
#' The dataframe_summary creates one table with feature_type, gene, locus_tag,
#' position, strand, TU, delay_fragment, HL_fragment,
#' half_life, intensity_fragment, intensity and velocity. The second table is
#' similar to the first one but in compact form.
#' It contains the same columns, the only difference is on position where a
#' start and end position are indicated separately.
#' Strand is indicated in case of stranded data to select the corresponding
#' positions.
#' 
#' @param inp SummarizedExperiment: the input data frame with correct format.
#' @param input dataframe: the probe based data frame with events and
#' gene annotation.
#' 
#' @return WIP
#' 
#' @examples
#' data(stats_minimal)
#' data(annot_g_minimal)
#' input <- event_dataframe(data = as.data.frame(rowRanges(stats_minimal)),
#' data_annotation = annot_g_minimal[[1]])
#' dataframe_summary(inp = stats_minimal, input = input)
#' @export
#' 
dataframe_summary <- function(inp, input) {
  tmp <-
    as.data.frame(
      rowRanges(inp)[, c(
      "ID",
      "position",
      "TU",
      "position_segment",
      "delay",
      "delay_fragment",
      "HL_fragment",
      "half_life",
      "HL_mean_fragment",
      "intensity",
      "intensity_mean_fragment",
      "intensity_fragment",
      "flag",
      "TI_termination_factor",
      "velocity_fragment"
    )])
  tmp <- tmp[,-c(1:4)]
  tmp_event <- input[, c("region", "gene", "locus_tag")]
  tmp_merged <- cbind(tmp, tmp_event)
  tmp_merged <-
    as.data.frame(tmp_merged %>% mutate_if(is.numeric, round, digits = 2))
  tmp_merged <-
    tmp_merged[grep("\\TU_\\d+$", tmp_merged$TU), ]
  tmp_merged <-
    tmp_merged[, c(
      "ID",
      "region",
      "gene",
      "locus_tag",
      "position",
      "strand",
      "TU",
      "position_segment",
      "delay",
      "delay_fragment",
      "HL_fragment",
      "half_life",
      "HL_mean_fragment",
      "intensity",
      "intensity_mean_fragment",
      "intensity_fragment",
      "flag",
      "TI_termination_factor",
      "velocity_fragment"
    )]
  tmp_merged[grep("_O$", tmp_merged$HL_fragment), "HL_mean_fragment"] <-
    tmp_merged[grep("_O$", tmp_merged$HL_fragment), "half_life"]
  tmp_merged[grep("_O$", tmp_merged$intensity_fragment),
             "intensity_mean_fragment"] <-
    tmp_merged[grep("_O$", tmp_merged$intensity_fragment), "intensity"]
  tmp_merged <-
    tmp_merged[, c(
      "ID",
      "region",
      "gene",
      "locus_tag",
      "position",
      "strand",
      "position_segment",
      "TU",
      "delay_fragment",
      "delay",
      "HL_fragment",
      "half_life",
      "HL_mean_fragment",
      "intensity",
      "intensity_mean_fragment",
      "intensity_fragment",
      "flag",
      "TI_termination_factor",
      "velocity_fragment"
    )]
  colnames(tmp_merged) <-
    c(
      "ID",
      "feature_type",
      "gene",
      "locus_tag",
      "position",
      "strand",
      "segment",
      "TU",
      "delay_fragment",
      "delay",
      "HL_fragment",
      "half_life",
      "HL_mean_fragment",
      "intensity",
      "intensity_mean_fragment",
      "intensity_fragment",
      "flag",
      "TI_termination_factor",
      "velocity_fragment"
    )
  tmp_df <-
    tmp_merged[, c(
      "ID",
      "feature_type",
      "gene",
      "locus_tag",
      "position",
      "strand",
      "segment",
      "TU",
      "delay_fragment",
      "delay",
      "HL_fragment",
      "half_life",
      "intensity_fragment",
      "intensity",
      "flag",
      "TI_termination_factor"
    )]
  int_frg <- unique(tmp_merged$intensity_fragment)
  int_frg <- int_frg[grep(paste0("\\I_", "\\d+", "$"), int_frg)]
  df <- data.frame()
  for (i in seq_along(int_frg)) {
    d <-
      tmp_merged[which(tmp_merged$intensity_fragment == int_frg[i]),
                 ]
    df[i, "feature_type"] <-
      paste(unique(unlist(strsplit(
        d$feature_type, split = ";"
      ))), collapse = "|")
    df[i, "gene"] <-
      paste(unique(unlist(strsplit(d$gene, split = ";"))), collapse = "|")
    df[i, "locus_tag"] <-
      paste(unique(unlist(strsplit(
        d$locus_tag, split = ";"
      ))), collapse = "|")
    df[i, "first_position_frg"] <-
      tmp[which(tmp$intensity_fragment == int_frg[i]), "position"][1]
    df[i, "last_position_frg"] <-
      last(tmp[which(tmp$intensity_fragment == int_frg[i]), "position"])
    df[i, "strand"] <- unique(d$strand)
    df[i, "TU"] <- unique(d$TU)
    df[i, "segment"] <- unique(d$segment)
    delay_frg <-
      d[grep("\\D_\\d+$", d$delay_fragment)[1], "delay_fragment"]
    if (is.na(delay_frg)) {
      df[i, "delay_fragment"] <- NA
    } else{
      df[i, "delay_fragment"] <-
        unique(tmp_merged[which(tmp_merged$delay_fragment %in% delay_frg),
                          "delay_fragment"])
    }
    hl_frg <- d[grep("\\Dc_\\d+$", d$HL_fragment)[1], "HL_fragment"]
    if (is.na(hl_frg)) {
      df[i, "HL_fragment"] <- NA
      df[i, "half_life"] <- NA
      df[i, "HL_SD"] <- NA
      df[i, "HL_SE"] <- NA
    } else{
      df[i, "HL_fragment"] <-
        unique(tmp_merged[which(tmp_merged$HL_fragment %in% hl_frg),
                          "HL_fragment"])
      df[i, "half_life"] <-
        unique(tmp_merged[which(tmp_merged$HL_fragment %in% hl_frg),
                          "HL_mean_fragment"])
      df[i, "HL_SD"] <-
        sd(tmp_merged[which(tmp_merged$HL_fragment %in% hl_frg), "half_life"])
      df[i, "HL_SE"] <-
        se(tmp_merged[which(tmp_merged$HL_fragment %in% hl_frg), "half_life"])
    }
    df[i, "intensity_fragment"] <- int_frg[i]
    df[i, "intensity"] <- unique(d$intensity_mean_fragment)
    df[i, "intensity_SD"] <-
      sd(tmp_merged[which(tmp_merged$intensity_fragment %in% int_frg[i]),
                    "intensity"])
    df[i, "intensity_SE"] <-
      se(tmp_merged[which(tmp_merged$intensity_fragment %in% int_frg[i]),
                    "intensity"])
    df[i, "velocity"] <- unique(d$velocity)
  }
  df <-
    as.data.frame(df %>% mutate_if(is.numeric, round, digits = 2))
  tables <- list(tmp_df, df)
  return(tables)
}
