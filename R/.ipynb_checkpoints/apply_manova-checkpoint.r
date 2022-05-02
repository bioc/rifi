# =========================================================================
# apply_manova   Checks if the ratio of hl ratio and intensity ratio 
#'               is statistically significant
# -------------------------------------------------------------------------
#'
#'
#' apply_manova compares the variance between two fold-changes, HL and intensity
#' within the same TU (half-life frgA/half-life frgB/
#' intensity frgA/intensity frgB).
#' HL fragment could cover two intensity fragments therefore this function sets
#' first fragments borders and uses manova_function.
#' Manova checks the variance between 2 segments (independent variables) and two
#' dependents variables (HL and intensity).
#'
#' @param inp SummarizedExperiment: the input data frame with correct format.
#'
#' @return the probe data frame with the columns regarding statistics:
#' \describe{
#'   \item{position:}{Integer, the bin/probe specific position}
#'   \item{intensity:}{Integer, the relative intensity of bin/probe at time point 0}
#'   \item{half_life:}{Integer, the half-life of the bin/probe}
#'   \item{HL_fragment:}{Integer, the half-life fragment the bin/probe belongs to}
#'   \item{HL_mean_fragment:}{Integer, the mean half-life value of the respective
#'   half-life fragment}
#'   \item{HL_fragment:}{Integer, the half_life fragment of the bin/probe belongs to}
#'   \item{intensity_fragment:}{Integer, the intensity fragment of the bin/probe belongs to}
#'   \item{TU:}{String, the overarching transcription unit}
#'   \item{FC_HL_intensity:}{Integer, The fold change of half-life/ fold change of
#'     intensity}
#'   \item{FC_HL_intensity_fragment:}{String, the fragments corresponding to intensity 
#'    and HL fold change}
#'   \item{p_value_Manova:}{Integer, the p_value added to the input}
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
