#' =========================================================================
#' dataframe_summary                                       
#' -------------------------------------------------------------------------
#' dataframe_summary creates two tables relating gene annotation to fragments  
#'
#' dataframe_summary creates two tables summary of segments and their
#' half-lives. The first output is bin/probe features and the second one is
#' intensity fragment based.The dataframe_summary creates one table with
#' feature_type, gene, locus_tag, position, strand, TU, delay_fragment,
#' HL_fragment, half_life, intensity_fragment, intensity and velocity. The 
#' second table is similar to the first one but in compact form. It contains the
#' same columns, the only difference is on position where a start and end 
#' position are indicated separately.
#' 
#' 
#' @param data SummarizedExperiment: the input data frame with correct format.
#' @param input dataframe: dataframe from event_dataframe function.
#' 
#' @return
#'   \item{bin_df:}{all information regarding bins:
#'   \describe{
#'     \item{position:}{Integer, position of the bin/probe on the genome}
#'     \item{feature_type:}{String, region annotation covering the fragments}
#'     \item{gene:}{String, gene annotation covering the fragments}
#'     \item{locus_tag:}{String, locus_tag annotation covering the fragments}
#'     \item{strand:}{Boolean. The bin/probe specific strand (+/-)}
#'     \item{segment:}{String, the bin/probe segment on the genome}
#'     \item{TU:}{String, The overarching transcription unit}
#'     \item{delay_fragment:}{The delay fragment the bin belongs to}
#'     \item{delay:}{Integer, the delay value of the bin/probe}
#'     \item{HL_fragment:}{The half-life fragment the bin belongs to}
#'     \item{half_life:}{The half-life of the bin/probe}
#'     \item{intensity_fragment:}{The intensity fragment the bin belongs to}
#'     \item{intensity:}{The relative intensity at time point 0}
#'     \item{flag:}{String, the flag of the bin/probe, contains information 
#'                  or the distribution for the #'different fitting models}
#'     \item{TI_termination_factor:}{String, the TI termination factor determined by TI}
#'     }
#'   }

#'   \item{frag_df:}{all information regarding fragments:
#'   \describe{
#'     \item{feature_type:}{String, region annotation covering the fragments}
#'     \item{gene:}{String, gene annotation covering the fragments}
#'     \item{locus_tag:}{String, locus_tag annotation covering the fragments}
#'     \item{first_position_frg:}{Integer, the bin/probe specific first position}
#'     \item{last_position_frg:}{Integer, the bin/probe specific last position}
#'     \item{strand:}{Boolean. The bin/probe specific strand (+/-)}
#'     \item{TU:}{String, The overarching transcription unit}
#'     \item{segment:}{String, the bin/probe segment on the genome}
#'     \item{delay_fragment:}{String, the delay fragment the bin belongs to}
#'     \item{HL_fragment:}{Integer, the half_life fragment of the bin/probe belongs to}
#'     \item{half_life:}{Integer, the half-life of the bin/probe}
#'     \item{HL_SD:}{Integer, the half-life standard deviation of the HL fragment, bin/probe based}
#'     \item{HL_SE:}{Integer, the half-life standard error of the HL fragment, bin/probe based}
#'     \item{intensity_fragment:}{Integer, the intensity fragment the bin belongs to}
#'     \item{intensity:}{Integer, the relative intensity of bin/probe at time point 0}
#'     \item{intensity_SD:}{Integer, the intensity standard deviation of the intensity fragment, bin/probe                  based}
#'     \item{intensity_SE:}{Integer, the intensity standard error of the intensity fragment, bin/probe based}
#'     \item{velocity:}{The velocity value of the respective delay fragment}
#'     }
#'   }
#' 
#' @examples
#' data(stats_minimal)
#' data(res_minimal)
#' dataframe_summary(data = stats_minimal, input = res_minimal)
#' @export
#' 
dataframe_summary <- function(data, input) {
  tmp <-
    as.data.frame(
      rowRanges(data)[, c(
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
  if(length(int_frg) != 0){
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
    }
  df <-
    as.data.frame(df %>% mutate_if(is.numeric, round, digits = 2))
  tables <- list(tmp_df, df)
  return(tables)
}
