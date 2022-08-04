#'=========================================================================
#' apply_manova
#'-------------------------------------------------------------------------
#' apply_manova checks if the ratio of hl ratio and intensity ratio is
#' statistically significant. 
#' 
#' apply_manova compares the variance between two fold-changes HL and intensity
#' within the same TU (half-life frgA/half-life frgB/intensity
#' frgA/intensity frgB). HL fragment could cover two intensity fragments 
#' therefore this function sets first fragments borders and uses manova_function.
#' Manova checks the variance between 2 segments (independent variables) and two
#' dependents variables (HL and intensity).

#' @param inp SummarizedExperiment: the input data frame with correct format.
#'
#' @return The probe data frame with the columns regarding statistics:
#' \describe{
#'   \item{ID:}{The bin/probe specific ID.}
#'   \item{position:}{The bin/probe specific position.}
#'   \item{strand:}{The bin/probe specific strand.}
#'   \item{intensity:}{The relative intensity at time point 0.}
#'   \item{probe_TI:}{An internal value to determine which fitting model
#'   is applied.}
#'   \item{flag:}{Information on which fitting model is applied.}
#'   \item{position_segment:}{The position based segment.}
#'   \item{delay:}{The delay value of the bin/probe.}
#'   \item{half_life:}{The half-life of the bin/probe.}
#'   \item{TI_termination_factor:}{String, the factor of TI fragment.}
#'   \item{delay_fragment:}{The delay fragment the bin belongs to.}
#'   \item{velocity_fragment:}{The velocity value of the respective delay
#'   fragment.}
#'   \item{intercept:}{The vintercept of fit through the respective delay
#'   fragment.}
#'   \item{slope:}{The slope of the fit through the respective delay fragment.}
#'   \item{HL_fragment:}{The half-life fragment the bin belongs to.}
#'   \item{HL_mean_fragment:}{The mean half-life value of the respective
#'   half-life fragment.}
#'   \item{intensity_fragment:}{The intensity fragment the bin belongs to.}
#'   \item{intensity_mean_fragment:}{The mean intensity value of the respective
#'   intensity fragment.}
#'   \item{TU:}{The overarching transcription unit.}
#'   \item{TI_termination_fragment:}{The TI fragment the bin belongs to.}
#'   \item{TI_mean_termination_factor:}{The mean termination factor of the
#'   respective TI fragment.}
#'   \item{seg_ID:}{The combined ID of the fragment.}
#'   \item{pausing_site:}{presence of pausing site indicated by +/-.}
#'   \item{iTSS_I:}{presence of iTSS_I indicated by +/-.}
#'   \item{ps_ts_fragment:}{The fragments involved in pausing site or iTSS_I.}
#'   \item{event_duration:}{Integer, the duration between two delay fragments.}
#'   \item{event_ps_itss_p_value_Ttest:}{p_value of pausing site or iTSS_I.}
#'   \item{p_value_slope:}{Integer, the p_value added to the inp.}
#'   \item{delay_frg_slope:}{Integer, the slope value of the fit through the 
#'   respective delay fragment.}
#'   \item{velocity_ratio:}{Integer, the ratio value of velocity from 2 delay 
#'   fragments.}
#'   \item{event_position:}{Integer, position of the event added to the input.}
#'   \item{FC_HL:}{Integer, the fold change value of 2 HL fragments.}
#'   \item{FC_fragment_HL:}{String, the fragments corresponding to HL fold 
#'   change.}
#'   \item{p_value_HL:}{Integer, the p_value added to the input of 2 HL 
#'   fragments.}
#'   \item{FC_intensity:}{Integer, the fold change value of 2 intensity 
#'   fragments.}
#'   \item{FC_fragment_intensity:}{String, the fragments corresponding to 
#'   intensity fold change.}
#'   \item{p_value_intensity:}{Integer, the p_value added to the input of 2 
#'   intensity fragments.}
#'   \item{synthesis_ratio:}{Integer, the value correspomding to synthesis rate.}
#'   \item{synthesis_ratio_event:}{String, the event assigned by synthesis 
#'   rate either Termination or iTSS.}
#'   \item{FC_HL_intensity:}{Integer, the value corresponding to HL and
#'   intensity fold change.}
#'   \item{FC_HL_intensity_fragment:}{String, the fragments corresponding to 
#'   intensity and HL fold change.}
#'   \item{FC_HL_adapted:}{Integer, the fold change of half-life/ fold change 
#'   of intensity,position of the half-life fragment is adapted to intensity 
#'   fragment.}
#'   \item{p_value_Manova:}{Integer, the p_value added to the input.}
#' }
#'
#' @examples
#' data(stats_minimal)
#' apply_manova(inp = stats_minimal)
#'
#' @export

# II. Manova statistical test
apply_manova <- function(inp) {
  rowRanges(inp)$p_value_Manova <- NA
  # select unique patterns subjected to ratio of fold change half-life...
  # ...and fold change intensity
  unique_Int_HL <-
    unique(na.omit(rowRanges(inp)$FC_HL_intensity_fragment))
  for (i in seq_along(unique_Int_HL)) {
    intHL <-
      rowRanges(inp)[which(rowRanges(inp)$FC_HL_intensity_fragment %in%
                             unique_Int_HL[i]),]
    # select half-life and intensity fragments
    frag_hl <- unique(intHL$FC_HL_intensity_fragment)
    frag_hl <- unlist(strsplit(frag_hl, ";"))
    frag_hl <- unlist(strsplit(frag_hl, ":"))
    # select the corresponding intensity fragments to the unique...
    # ...half-life fragment
    I_1 <-
      rowRanges(inp)[which(rowRanges(inp)$intensity_fragment %in%
                             frag_hl[length(frag_hl) - 1]),
                     c("intensity_fragment", "intensity", "position")]
    hl_1 <- rowRanges(inp)[which(
      rowRanges(inp)$position %in%
        I_1$position &
        rowRanges(inp)$intensity_fragment %in%
        I_1$intensity_fragment
    ), c("HL_fragment", "half_life", "position")]
    
    hl_1 <-
      hl_1[grep(paste0("\\Dc_", "\\d+", "$"), hl_1$HL_fragment),]
    
    I_1 <- I_1[which(hl_1$position %in% I_1$position),]
    
    I_2 <-
      rowRanges(inp)[which(rowRanges(inp)$intensity_fragment %in%
                             frag_hl[length(frag_hl)]),
                     c("intensity_fragment", "intensity", "position")]
    
    hl_2 <- rowRanges(inp)[which(
      rowRanges(inp)$position %in%
        I_2$position &
        rowRanges(inp)$intensity_fragment %in%
        I_2$intensity_fragment
    ), c("HL_fragment", "half_life", "position")]
    
    hl_2 <-
      hl_2[grep(paste0("\\Dc_", "\\d+", "$"), hl_2$HL_fragment),]
    
    I_2 <- I_2[which(hl_2$position %in% I_2$position),]
    
    I_1$segment <- "S1"
    I_2$segment <- "S2"
    hl_1$segment <- "S1"
    hl_2$segment <- "S2"
    if (length(I_1$segment) < 2 |
        length(I_2$segment) < 2 | 
        length(hl_1$segment) < 2 |
        length(hl_2$segment) < 2)
    {
      next ()
    }
    df <- cbind.data.frame(c(hl_1$half_life, hl_2$half_life),
                           c(I_1[, c("intensity", "segment")],
                             I_2[, c("intensity", "segment")]))
    colnames(df)[1] <- "half_life"
    tryCatch({
      p_value_Manova <-
        manova_function(
          x = df$half_life,
          y = df$intensity,
          z = df$segment,
          data = df
        )
      rowRanges(inp)$p_value_Manova[which(
        rowRanges(inp)$FC_HL_intensity_fragment %in% unique_Int_HL[i])] <-
        p_value_Manova
    }, error = function(e) {
    })
  }
  return(inp)
}