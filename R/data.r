#' An example SummarizedExperiment from E. coli 
#' An example SummarizedExperiment from RNA-seq containing information about the
#' intensities at all time points (assay). Seqnames, IRanges and strand columns
#' (rowRanges)and colData with time point series and replicates. 
#'
#' @format A assay:
#' \describe{
#'   \item{0:}{relative intensities at 0 min}
#'   \item{1:}{relative intensities at 1 min}
#'   \item{10:}{relative intensities at 10 min}
#'   \item{15:}{relative intensities at 15 min}
#'   \item{2:}{relative intensities at 2 min}
#'   \item{20:}{relative intensities at 20 min}
#'   \item{3:}{relative intensities at 3 min}
#'   \item{4:}{relative intensities at 4 min}
#'   \item{5:}{relative intensities at 5 min}
#'   \item{6:}{relative intensities at 6 min}
#'   \item{8:}{relative intensities at 8 min}
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(example_input_e_coli)
#'
"example_input_e_coli"

#' An artificial example SummarizedExperiment
#' An example SummarizedExperiment containing information about the intensities
#' at all time points (assay). Seqnames, IRanges and strand columns (rowRanges)
#' and colData with time point series and replicates.
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(example_input_minimal)
#'
"example_input_minimal"

#' An example input data frame from Synechocystis PCC 6803
#' A SummarizedExperiment from microarrays data containing information about the
#' intensities at all time points (assay), Seqnames, IRanges and strand columns
#' (rowRanges) and colData with time point series and averaged replicates.
#' 
#' @format Assay with 3000 rows and 10 variables:
#' \describe{
#'   \item{0:}{relative intensities at 0 min}
#'   \item{2:}{relative intensities at 2 min}
#'   \item{4:}{relative intensities at 4 min}
#'   \item{8:}{relative intensities at 8 min}
#'   \item{16:}{relative intensities at 16 min}
#'   \item{32:}{relative intensities at 32 min}
#'   \item{64:}{relative intensities at 64 min}
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(example_input_synechocystis_6803)
#'
"example_input_synechocystis_6803"

#' The result of rifi_fit for E.coli example data
#' A SummarizedExperiment containing the output from rifi_fit as an extension of
#' rowRanges and metadata.
#' 
#' @format Three data frames with 290 rows and 10 variables, 155 rows
#' and 5 variables, and 135 rows and 9 variables are generated. The columns of
#' the first data frame are added to the rowRanges and the rest are added as
#' metadata.
#' \describe{
#'   \item{inp:}{The SummarizedExperiment:
#'   \describe{
#'     \item{ID:}{The bin/probe specific ID}
#'     \item{position:}{The bin/probe specific position}
#'     \item{intensity:}{The relative intensity at time point 0}
#'     \item{probe_TI:}{An internal value to determine which fitting model
#'     is applied}
#'     \item{flag:}{Information on which fitting model is applied}
#'     \item{postion_segment:}{The position based segment}
#'     \item{delay:}{The delay value of the bin/probe}
#'     \item{half_life:}{The half-life of the bin/probe}
#'     \item{TI_termination_factor:}{The termination factor of the bin/probe}
#'     }
#'   }
#'   \item{fit_obj_STD:}{the fit object for the standard fit:
#'   \describe{
#'     \item{ID:}{The bin/probe specific ID}
#'     \item{delay:}{The delay value of the bin/probe}
#'     \item{half_life:}{The half-life of the bin/probe}
#'     \item{inty_S0:}{The relative intensity at time point 0}
#'     \item{intyf:}{The background value of the fit}
#'     }
#'   }
#'   \item{fit_obj_TI:}{the fit object for the TI fit:
#'   \describe{
#'     \item{delay:}{The delay value of the bin/probe}
#'     \item{ti_delay:}{The ti-delay value of the bin/probe}
#'     \item{half_life:}{The half-life of the bin/probe}
#'     \item{ti_value:}{The ti-value of the bin/probe}
#'     \item{TI_termination_factor:}{The termination factor of the bin/probe}
#'     \item{synthesis_rate:}{The synthesis rate of the bin/probe}
#'     \item{TI_background:}{The background value of the fit}
#'     \item{position:}{The bin/probe specific position}
#'     \item{ID:}{The bin/probe specific ID}
#'     }
#'   }
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(fit_e_coli)
#'
"fit_e_coli"

#' The artificial result of rifi_fit for artificial example data
#' A SummarizedExperiment containing the output from rifi_fit.
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(fit_minimal)
#'
"fit_minimal"

#' The result of rifi_fit for Synechocystis 6803 example data
#' A SummarizedExperiment containing the output from rifi_fit as an extension of
#' rowRanges and metadata.
#' 
#' @format Three data frames with 3000 rows and 10 variables, 2811 rows
#' and 5 variables, and 189 rows and 9 variable are generated. The columns of
#' the first data frame are added to the rowRanges and the rest are added as
#' metadata.
#' \describe{
#'   \item{inp:}{the SummarizedExperiment:
#'   \describe{
#'     \item{ID:}{The bin/probe specific ID}
#'     \item{position:}{The bin/probe specific position}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{intensity:}{The relative intensity at time point 0}
#'     \item{probe_TI:}{An internal value to determine which fitting model
#'     is applied}
#'     \item{flag:}{Information on which fitting model is applied}
#'     \item{postion_segment:}{The position based segment}
#'     \item{delay:}{The delay value of the bin/probe}
#'     \item{half_life:}{The half-life of the bin/probe}
#'     \item{TI_termination_factor:}{The termination factor of the bin/probe}
#'     }
#'   }
#'   \item{fit_obj_STD:}{the fit object for the standard fit:
#'   \describe{
#'     \item{ID:}{The bin/probe specific ID}
#'     \item{delay:}{The delay value of the bin/probe}
#'     \item{half_life:}{The half-life of the bin/probe}
#'     \item{inty_S0:}{The relative intensity at time point 0}
#'     \item{intyf:}{The background value of the fit}
#'     }
#'   }
#'   \item{fit_obj_TI:}{the fit object for the TI fit:
#'   \describe{
#'     \item{delay:}{The delay value of the bin/probe}
#'     \item{ti_delay:}{The ti-delay value of the bin/probe}
#'     \item{half_life:}{The half-life of the bin/probe}
#'     \item{ti_value:}{The ti-value of the bin/probe}
#'     \item{TI_termination_factor:}{The termination factor of the bin/probe}
#'     \item{synthesis_rate:}{The synthesis rate of the bin/probe}
#'     \item{TI_background:}{The background value of the fit}
#'     \item{position:}{The bin/probe specific position}
#'     \item{ID:}{The bin/probe specific ID}
#'     }
#'   }
#' }
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(fit_synechocystis_6803)
#'
"fit_synechocystis_6803"

#' The result of rifi_fragmentation for E.coli example data
#' A SummarizedExperiment containing the output from rifi_fragmentation as an
#' extension of rowRanges
#' 
#' @format rowRanges of the SummarizedExperiment with 290 rows and 22 variables:
#' \describe{
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{position:}{The bin/probe specific position}
#'   \item{intensity:}{The relative intensity at time point 0}
#'   \item{probe_TI:}{An internal value to determine which fitting model is
#'    applied}
#'   \item{flag:}{Information on which fitting model is applied}
#'   \item{position_segment:}{The position based segment}
#'   \item{delay:}{The delay value of the bin/probe}
#'   \item{half_life:}{The half-life of the bin/probe}
#'   \item{TI_termination_factor:}{The termination factor of the bin/probe}
#'   \item{delay_fragment:}{The delay fragment the bin belongs to}
#'   \item{velocity_fragment:}{The velocity value of the respective delay
#'   fragment}
#'   \item{intercept:}{The vintercept of fit through the respective delay
#'   fragment}
#'   \item{slope:}{The slope of the fit through the respective delay fragment}
#'   \item{HL_fragment:}{The half-life fragment the bin belongs to}
#'   \item{HL_mean_fragment:}{The mean half-life value of the respective
#'   half-life fragment}
#'   \item{intensity_fragment:}{The intensity fragment the bin belongs to}
#'   \item{intensity_mean_fragment:}{The mean intensity value of the respective
#'   intensity fragment}
#'   \item{TU:}{The overarching transcription unit}
#'   \item{TI_termination_fragment:}{The TI fragment the bin belongs to}
#'   \item{TI_mean_termination_factor:}{The mean termination factor of the
#'   respective TI fragment}
#'   \item{seg_ID:}{The combined ID of the fragment}
#' }
#' 
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(fragmentation_e_coli)
#'
"fragmentation_e_coli"

#' The result of rifi_fragmentation for artificial example data
#' A SummarizedExperiment containing the output from rifi_fragmentation as an
#' extension of rowRanges and metadata.
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(fragmentation_minimal)
#'
"fragmentation_minimal"

#' The result of rifi_fragmentation for Synechocystis 6803 example data
#' A SummarizedExperiment containing the output from rifi_fragmentation as an
#' extension fo rowRanges
#' 
#' @format rowRanges of the SummarizedExperiment:
#' \describe{
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{position:}{The bin/probe specific position}
#'   \item{intensity:}{The relative intensity at time point 0}
#'   \item{probe_TI:}{An internal value to determine which fitting model is
#'   applied}
#'   \item{flag:}{Information on which fitting model is applied}
#'   \item{position_segment:}{The position based segment}
#'   \item{delay:}{The delay value of the bin/probe}
#'   \item{half_life:}{The half-life of the bin/probe}
#'   \item{TI_termination_factor:}{The termination factor of the bin/probe}
#'   \item{delay_fragment:}{The delay fragment the bin belongs to}
#'   \item{velocity_fragment:}{The velocity value of the respective delay
#'   fragment}
#'   \item{intercept:}{The vintercept of fit through the respective delay
#'   fragment}
#'   \item{slope:}{The slope of the fit through the respective delay fragment}
#'   \item{HL_fragment:}{The half-life fragment the bin belongs to}
#'   \item{HL_mean_fragment:}{The mean half-life value of the respective
#'   half-life fragment}
#'   \item{intensity_fragment:}{The intensity fragment the bin belongs to}
#'   \item{intensity_mean_fragment:}{The mean intensity value of the
#'   respective intensity fragment}
#'   \item{TU:}{The overarching transcription unit}
#'   \item{TI_termination_fragment:}{The TI fragment the bin belongs to}
#'   \item{TI_mean_termination_factor:}{The mean termination factor of the
#'   respective TI fragment}
#'   \item{seg_ID:}{The combined ID of the fragment}
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(fragmentation_synechocystis_6803)
#'
"fragmentation_synechocystis_6803"

#' The result of rifi_penalties for E.coli example data.
#' A SummarizedExperiment containing the output from rifi_penalties including
#' the logbook and the four penalty objects as metadata.
#' 
#' @format A list with 5 items:
#' \describe{
#'   \item{logbook:}{The logbook vector containing all penalty information}
#'   \item{pen_obj_delay:}{A list with 4 items:
#'     \describe{
#'       \item{logbook:}{The logbook vector containing all penalty information}
#'       \item{delay_penalties:}{a vetor with the delay penalty and delay
#'       outlier penalty}
#'       \item{correct:}{a matrix of the correct splits}
#'       \item{wrong:}{a matrix of the incorrect splits}
#'     }
#'   }
#'   \item{pen_obj_HL:}{A list with 4 items:
#'     \describe{
#'       \item{logbook:}{The logbook vector containing all penalty information}
#'       \item{HL_penalties:}{a vetor with the half-life penalty and half-life
#'       outlier penalty}
#'       \item{correct:}{a matrix of the correct splits}
#'       \item{wrong:}{a matrix of the incorrect splits}
#'     }
#'   }
#'   \item{pen_obj_inty:}{A list with 4 items:
#'     \describe{
#'       \item{logbook:}{The logbook vector containing all penalty information}
#'       \item{inty_penalties:}{a vetor with the intensity penalty and intensity
#'       outlier penalty}
#'       \item{correct:}{a matrix of the correct splits}
#'       \item{wrong:}{a matrix of the incorrect splits}
#'     }
#'   }
#'   \item{pen_obj_TI:}{A list with 4 items:
#'     \describe{
#'       \item{logbook:}{The logbook vector containing all penalty information}
#'       \item{TI_penalties:}{a vetor with the TI penalty and TI outlier
#'       penalty}
#'       \item{correct:}{a matrix of the correct splits}
#'       \item{wrong:}{a matrix of the incorrect splits}
#'     }
#'   }
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(penalties_e_coli)
#'
"penalties_e_coli"

#' The result of rifi_penalties for artificial example data
#' A SummarizedExperiment containing the output from rifi_penalties including
#' the logbook and the four penalty objects as metadata.
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(penalties_minimal)
#'
"penalties_minimal"

#' The result of rifi_penalties for Synechocystis 6803 example data.
#' A SummarizedExperiment containing the output from rifi_penalties including
#' the logbook and the four penalty objects as metadata.
#' 
#' @format A list with 5 items:
#' \describe{
#'   \item{logbook:}{The logbook vector containing all penalty information}
#'   \item{pen_obj_delay:}{A list with 4 items:
#'     \describe{
#'       \item{logbook:}{The logbook vector containing all penalty information}
#'       \item{delay_penalties:}{a vetor with the delay penalty and delay
#'       outlier penalty}
#'       \item{correct:}{a matrix of the correct splits}
#'       \item{wrong:}{a matrix of the incorrect splits}
#'     }
#'   }
#'   \item{pen_obj_HL:}{A list with 4 items:
#'     \describe{
#'       \item{logbook:}{The logbook vector containing all penalty information}
#'       \item{HL_penalties:}{a vetor with the half-life penalty and half-life
#'       outlier penalty}
#'       \item{correct:}{a matrix of the correct splits}
#'       \item{wrong:}{a matrix of the incorrect splits}
#'     }
#'   }
#'   \item{pen_obj_inty:}{A list with 4 items:
#'     \describe{
#'       \item{logbook:}{The logbook vector containing all penalty information}
#'       \item{inty_penalties:}{a vetor with the intensity penalty and intensity
#'       outlier penalty}
#'       \item{correct:}{a matrix of the correct splits}
#'       \item{wrong:}{a matrix of the incorrect splits}
#'     }
#'   }
#'   \item{pen_obj_TI:}{A list with 4 items:
#'     \describe{
#'       \item{logbook:}{The logbook vector containing all penalty information}
#'       \item{TI_penalties:}{a vetor with the TI penalty and TI outlier
#'       penalty}
#'       \item{correct:}{a matrix of the correct splits}
#'       \item{wrong:}{a matrix of the incorrect splits}
#'     }
#'   }
#' }
#' 
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(penalties_synechocystis_6803)
#'
"penalties_synechocystis_6803"

#' The result of rifi_preprocess for E.coli example data
#' A SummarizedExperiment containing the output from rifi_penalties including
#' the logbook and the four penalty objects as metadata.
#' A list containing the output from rifi_preprocess, including the inp
#' and the modified input_df.
#' 
#' @format A SummarizedExperiment:
#' \describe{
#'   \item{inp:}{the SummarizedExperiment:
#'   \describe{
#'     \item{ID:}{The bin/probe specific ID}
#'     \item{position:}{The bin/probe specific position}
#'     \item{intensity:}{The relative intensity at time point 0}
#'     \item{probe_TI:}{An internal value to determine which fitting model is
#'     applied}
#'     \item{flag:}{Information on which fitting model is applied}
#'     \item{postion_segment:}{The position based segment}
#'       }
#'     }
#'   \item{fit_obj_TI:}{the fit object for the TI fit:
#'     \describe{
#'     \item{0:}{relative intensities at 0 min}
#'     \item{1:}{relative intensities at 1 min}
#'     \item{10:}{relative intensities at 10 min}
#'     \item{15:}{relative intensities at 15 min}
#'     \item{2:}{relative intensities at 2 min}
#'     \item{20:}{relative intensities at 20 min}
#'     \item{3:}{relative intensities at 3 min}
#'     \item{4:}{relative intensities at 4 min}
#'     \item{5:}{relative intensities at 5 min}
#'     \item{6:}{relative intensities at 6 min}
#'     \item{8:}{relative intensities at 8 min}
#'     \item{ID:}{unique IDs}
#'     \item{position:}{genome positions}
#'     \item{filtration:}{indicator wether the replicate is filtered or not}
#'     }
#'   }
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(preprocess_e_coli)
#'
"preprocess_e_coli"

#' The result of rifi_preprocess for artificial example data
#' A SummarizedExperiment containing the output from rifi_preprocess
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(preprocess_minimal)
#'
"preprocess_minimal"

#' The result of rifi_preprocess for Synechocystis 6803 example data is a 
#' A SummarizedExperiment containing the output of rifi_preprocess as an
#' extention to rowRanges
#' 
#' @format A SummarizedExperiment:
#' \describe{
#'   \item{inp:}{the SummarizedExperiment:
#'   \describe{
#'     \item{ID:}{The bin/probe specific ID}
#'     \item{position:}{The bin/probe specific position}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{intensity:}{The relative intensity at time point 0}
#'     \item{probe_TI:}{An internal value to determine which fitting model
#'     is applied}
#'     \item{flag:}{Information on which fitting model is applied}
#'     \item{postion_segment:}{The position based segment}
#'       }
#'     }
#'   \item{fit_obj_TI:}{the fit object for the TI fit:
#'     \describe{
#'     \item{0:}{relative intensities at 0 min}
#'     \item{2:}{relative intensities at 2 min}
#'     \item{4:}{relative intensities at 4 min}
#'     \item{8:}{relative intensities at 8 min}
#'     \item{16:}{relative intensities at 16 min}
#'     \item{32:}{relative intensities at 32 min}
#'     \item{64:}{relative intensities at 64 min}
#'     \item{ID:}{unique IDs}
#'     \item{position:}{genome positions}
#'     \item{filtration:}{indicator wether the replicate is filtered or not}
#'     }
#'   }
#' }
#' 
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(preprocess_synechocystis_6803)
#'
"preprocess_synechocystis_6803"

#' The result of rifi_stats for E.coli example data
#' A SummarizedExperiment containing the output from rifi_stats
#' 
#' @format A SummarizedExperiment:
#' \describe{
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{position:}{The bin/probe specific position}
#'   \item{strand:}{The bin/probe specific strand}
#'   \item{intensity:}{The relative intensity at time point 0}
#'   \item{probe_TI:}{An internal value to determine which fitting model
#'   is applied}
#'   \item{flag:}{Information on which fitting model is applied}
#'   \item{position_segment:}{The position based segment}
#'   \item{delay:}{The delay value of the bin/probe}
#'   \item{half_life:}{The half-life of the bin/probe}
#'   \item{TI_termination_factor:}{The termination factor of the bin/probe}
#'   \item{delay_fragment:}{The delay fragment the bin belongs to}
#'   \item{velocity_fragment:}{The velocity value of the respective delay
#'   fragment}
#'   \item{intercept:}{The vintercept of fit through the respective delay
#'   fragment}
#'   \item{slope:}{The slope of the fit through the respective delay fragment}
#'   \item{HL_fragment:}{The half-life fragment the bin belongs to}
#'   \item{HL_mean_fragment:}{The mean half-life value of the respective
#'   half-life fragment}
#'   \item{intensity_fragment:}{The intensity fragment the bin belongs to}
#'   \item{intensity_mean_fragment:}{The mean intensity value of the respective
#'   intensity fragment}
#'   \item{TU:}{The overarching transcription unit}
#'   \item{TI_termination_fragment:}{The TI fragment the bin belongs to}
#'   \item{TI_mean_termination_factor:}{The mean termination factor of the
#'   respective TI fragment}
#'   \item{seg_ID:}{The combined ID of the fragment}
#'   \item{pausing_site:}{}
#'   \item{iTSS_I:}{}
#'   \item{ps_ts_fragment:}{}
#'   \item{event_ps_itss_p_value_Ttest:}{}
#'   \item{p_value_slope:}{}
#'   \item{delay_frg_slope:}{}
#'   \item{velocity_ratio:}{}
#'   \item{event_duration:}{}
#'   \item{event_position:}{}
#'   \item{FC_HL:}{}
#'   \item{FC_fragment_HL:}{}
#'   \item{p_value_HL:}{}
#'   \item{FC_intensity:}{}
#'   \item{FC_fragment_intensity:}{}
#'   \item{p_value_intensity:}{}
#'   \item{FC_HL_intensity:}{}
#'   \item{FC_HL_intensity_fragment:}{}
#'   \item{FC_HL_adapted:}{}
#'   \item{synthesis_ratio:}{}
#'   \item{synthesis_ratio_event:}{}
#'   \item{p_value_Manova:}{}
#'   \item{p_value_TI:}{}
#'   \item{TI_fragments_p_value:}{}
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(stats_e_coli)
#'
"stats_e_coli"

#' The result of rifi_stats for artificial example data
#' A SummarizedExperiment containing the output of rifi_stats as an
#' extention to rowRanges and metadata (gff file processed, see gff file
#' documentation)
#' 
#' @format A rowRanges of SummarizedExperiment with 24 rows and 45 variables:
#' \describe{
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{position:}{The bin/probe specific position}
#'   \item{intensity:}{The relative intensity at time point 0}
#'   \item{probe_TI:}{An internal value to determine which fitting model is
#'   applied}
#'   \item{flag:}{Information on which fitting model is applied}
#'   \item{position_segment:}{The position based segment}
#'   \item{delay:}{The delay value of the bin/probe}
#'   \item{half_life:}{The half-life of the bin/probe}
#'   \item{TI_termination_factor:}{The termination factor of the bin/probe}
#'   \item{delay_fragment:}{The delay fragment the bin belongs to}
#'   \item{velocity_fragment:}{The velocity value of the respective delay
#'   fragment}
#'   \item{intercept:}{The vintercept of fit through the respective delay
#'   fragment}
#'   \item{slope:}{The slope of the fit through the respective delay fragment}
#'   \item{HL_fragment:}{The half-life fragment the bin belongs to}
#'   \item{HL_mean_fragment:}{The mean half-life value of the respective
#'   half-life fragment}
#'   \item{intensity_fragment:}{The intensity fragment the bin belongs to}
#'   \item{intensity_mean_fragment:}{The mean intensity value of the respective
#'   intensity fragment}
#'   \item{TU:}{The overarching transcription unit}
#'   \item{TI_termination_fragment:}{The TI fragment the bin belongs to}
#'   \item{TI_mean_termination_factor:}{The mean termination factor of the
#'   respective TI fragment}
#'   \item{seg_ID:}{The combined ID of the fragment}
#'   \item{pausing_site:}{}
#'   \item{iTSS_I:}{}
#'   \item{ps_ts_fragment:}{}
#'   \item{event_ps_itss_p_value_Ttest:}{}
#'   \item{p_value_slope:}{}
#'   \item{delay_frg_slope:}{}
#'   \item{velocity_ratio:}{}
#'   \item{event_duration:}{}
#'   \item{event_position:}{}
#'   \item{FC_HL:}{}
#'   \item{FC_fragment_HL:}{}
#'   \item{p_value_HL:}{}
#'   \item{FC_intensity:}{}
#'   \item{FC_fragment_intensity:}{}
#'   \item{p_value_intensity:}{}
#'   \item{FC_HL_intensity:}{}
#'   \item{FC_HL_intensity_fragment:}{}
#'   \item{FC_HL_adapted:}{}
#'   \item{synthesis_ratio:}{}
#'   \item{synthesis_ratio_event:}{}
#'   \item{p_value_Manova:}{}
#'   \item{p_value_TI:}{}
#'   \item{TI_fragments_p_value:}{}
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(stats_minimal)
#'
"stats_minimal"

#' The result of rifi_stats for Synechocystis 6803 example data
#' A SummarizedExperiment containing the output of rifi_stats as an
#' extention to rowRanges
#' 
#' @format The rowRanges of SummarizedExperiment:
#' \describe{
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{position:}{The bin/probe specific position}
#'   \item{intensity:}{The relative intensity at time point 0}
#'   \item{probe_TI:}{An internal value to determine which fitting model is
#'   applied}
#'   \item{flag:}{Information on which fitting model is applied}
#'   \item{position_segment:}{The position based segment}
#'   \item{delay:}{The delay value of the bin/probe}
#'   \item{half_life:}{The half-life of the bin/probe}
#'   \item{TI_termination_factor:}{The termination factor of the bin/probe}
#'   \item{delay_fragment:}{The delay fragment the bin belongs to}
#'   \item{velocity_fragment:}{The velocity value of the respective delay
#'   fragment}
#'   \item{intercept:}{The vintercept of fit through the respective delay
#'   fragment}
#'   \item{slope:}{The slope of the fit through the respective delay fragment}
#'   \item{HL_fragment:}{The half-life fragment the bin belongs to}
#'   \item{HL_mean_fragment:}{The mean half-life value of the respective
#'   half-life fragment}
#'   \item{intensity_fragment:}{The intensity fragment the bin belongs to}
#'   \item{intensity_mean_fragment:}{The mean intensity value of the respective
#'   intensity fragment}
#'   \item{TU:}{The overarching transcription unit}
#'   \item{TI_termination_fragment:}{The TI fragment the bin belongs to}
#'   \item{TI_mean_termination_factor:}{The mean termination factor of the
#'   respective TI fragment}
#'   \item{seg_ID:}{The combined ID of the fragment}
#'   \item{pausing_site:}{}
#'   \item{iTSS_I:}{}
#'   \item{ps_ts_fragment:}{}
#'   \item{event_ps_itss_p_value_Ttest:}{}
#'   \item{p_value_slope:}{}
#'   \item{delay_frg_slope:}{}
#'   \item{velocity_ratio:}{}
#'   \item{event_duration:}{}
#'   \item{event_position:}{}
#'   \item{FC_HL:}{}
#'   \item{FC_fragment_HL:}{}
#'   \item{p_value_HL:}{}
#'   \item{FC_intensity:}{}
#'   \item{FC_fragment_intensity:}{}
#'   \item{p_value_intensity:}{}
#'   \item{FC_HL_intensity:}{}
#'   \item{FC_HL_intensity_fragment:}{}
#'   \item{FC_HL_adapted:}{}
#'   \item{synthesis_ratio:}{}
#'   \item{synthesis_ratio_event:}{}
#'   \item{p_value_Manova:}{}
#'   \item{p_value_TI:}{}
#'   \item{TI_fragments_p_value:}{}
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(stats_synechocystis_6803)
#'
"stats_synechocystis_6803"

#' The result of rifi_summary for E.coli example data
#' A SummarizedExperiment containing the output of rifi_stats as an
#' extention to rowRanges
#' 
#' @format The rowRanges of SummarizedExperiment:
#' \describe{
#'   \item{bin_df:}{all information regarding bins:
#'   \describe{
#'     \item{ID:}{}
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{position:}{}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{segment:}{The segment the bin/probe belongs to}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{delay_fragment:}{The delay fragment the bin/probe belongs to}
#'     \item{delay:}{The delay of the bin/probe}
#'     \item{HL_fragment:}{The half-life fragment the bin/probe belongs to}
#'     \item{half_life:}{The half-life of the bin/probe}
#'     \item{intensity_fragment:}{The intensity fragment the bin/probe belongs
#'      to}
#'     \item{intensity:}{The relative intensity at time point 0}
#'     \item{flag:}{The flag of the bin/probe(TI, PDD)}
#'     \item{TI_termination_factor:}{The TI_termination_factor of the bin/probe 
#'     (in case TI is detected)}
#'     }
#'   }
#'   \item{frag_df:}{all information regarding fragments:
#'   \describe{
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{first_position_frg:}{The first position of the fragment on the
#'      genome}
#'     \item{last_position_frg:}{The last position of the fragment on the
#'      genome}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{segment:}{The segment the fragment belongs to}
#'     \item{delay_fragment:}{The delay fragment of the fragment}
#'     \item{HL_fragment:}{The half-life fragment of the fragment}
#'     \item{half_life:}{The half-life mean of the fragment}
#'     \item{HL_SD:}{The half-life standard deviation of the fragment}
#'     \item{HL_SE:}{The half-life standard error of the fragment}
#'     \item{intensity_fragment:}{The intensity_fragment of the fragment}
#'     \item{intensity:}{The relative intensity at time point 0}
#'     \item{intensity_SD:}{The intensity standard deviation of the fragment}
#'     \item{intensity_SE:}{The intensity standard error of the fragment}
#'     \item{velocity:}{The velocity value of the respective delay fragment}
#'     }
#'   }
#'   \item{event_df:}{all information regarding events:
#'   \describe{
#'     \item{event:}{}
#'     \item{p_value:}{}
#'     \item{p_adjusted:}{}
#'     \item{FC_HL:}{Fold change of half-life}
#'     \item{FC_intensity:}{Fold change of intensity}
#'     \item{FC_HL_adapted:}{Fold change of half-life/ fold change of intensity,
#'     position of the half-life fragment is adapted to intensity fragment}
#'     \item{FC_HL_FC_intensity:}{Fold change of half-life/ fold change of
#'     intensity}
#'     \item{event_position:}{}
#'     \item{velocity_ratio:}{}
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{segment_1:}{}
#'     \item{segment_2:}{}
#'     \item{event_duration:}{}
#'     \item{gap_fragments:}{}
#'     \item{features:}{}
#'     }
#'   }
#'   \item{events_HL_int_df:}{all information regarding events related to 
#'   half-life and intensity:
#'   \describe{
#'     \item{event:}{}
#'     \item{p_value:}{}
#'     \item{p_adjusted:}{}
#'     \item{FC_HL:}{}
#'     \item{FC_intensity:}{}
#'     \item{FC_HL_adapted:}{Fold change of half-life/ fold change of intensity,
#'     position of the half-life fragment is adapted to intensity fragment}
#'     \item{FC_HL_FC_intensity:}{Fold change of half-life/ fold change of
#'     intensity}
#'     \item{event_position:}{}
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{segment_1:}{}
#'     \item{segment_2:}{}
#'     \item{event_duration:}{}
#'     \item{gap_fragments:}{}
#'     \item{features:}{}
#'     }
#'   }
#'   \item{events_ps_itss_df:}{all information regarding events related to 
#'   pausing sites and iTSS_I:
#'   \describe{
#'     \item{event:}{}
#'     \item{p_value:}{}
#'     \item{p_adjusted:}{}
#'     \item{event_position:}{}
#'     \item{velocity_ratio:}{}
#'     \item{FC_HL_adapted:}{}
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{segment_1:}{}
#'     \item{segment_2:}{}
#'     \item{event_duration:}{}
#'     \item{gap_fragments:}{}
#'     \item{features:}{}
#'     }
#'   }
#'   \item{events_velocity_df:}{all information regarding events related to 
#'   velocity:
#'   \describe{
#'     \item{event:}{}
#'     \item{p_value:}{}
#'     \item{p_adjusted:}{}
#'     \item{event_position:}{}
#'     \item{velocity_ratio:}{}
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{segment_1:}{}
#'     \item{segment_2:}{}
#'     \item{event_duration:}{}
#'     \item{gap_fragments:}{}
#'     \item{features:}{}
#'     }
#'   }
#'   \item{TI_df:}{all information regarding TI:
#'   \describe{
#'     \item{event:}{}
#'     \item{TI_fragment:}{}
#'     \item{TI_termination_factor:}{}
#'     \item{p_value:}{}
#'     \item{p_adjusted:}{}
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{features:}{}
#'     \item{event_position:}{}
#'     \item{position_1:}{}
#'     \item{position_2:}{}
#'     }
#'   }
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(summary_e_coli)
#'
"summary_e_coli"
#'
#' The result of rifi_summary for artificial example data
#' A SummarizedExperiment with the output from rifi_summary as metadata
#'
#' @format A list of 7 data frames with 290 rows and 11 variables, 36 rows
#' and 11 variables, 57 rows and 18 variables, and 8 rows and 14 variables:
#' \describe{
#'   \item{bin_df:}{all information regarding bins:
#'   \describe{
#'     \item{ID:}{}
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{position:}{}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{segment:}{The segment the bin/probe belongs to}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{delay_fragment:}{The delay fragment the bin/probe belongs to}
#'     \item{delay:}{The delay of the bin/probe}
#'     \item{HL_fragment:}{The half-life fragment the bin/probe belongs to}
#'     \item{half_life:}{The half-life of the bin/probe}
#'     \item{intensity_fragment:}{The intensity fragment the bin/probe belongs
#'      to}
#'     \item{intensity:}{The relative intensity at time point 0}
#'     \item{flag:}{The flag of the bin/probe(TI, PDD)}
#'     \item{TI_termination_factor:}{The TI_termination_factor of the bin/probe 
#'     (in case TI is detected)}
#'     }
#'   }
#'   \item{frag_df:}{all information regarding fragments:
#'   \describe{
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{first_position_frg:}{The first position of the fragment on the
#'      genome}
#'     \item{last_position_frg:}{The last position of the fragment on the
#'      genome}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{segment:}{The segment the fragment belongs to}
#'     \item{delay_fragment:}{The delay fragment of the fragment}
#'     \item{HL_fragment:}{The half-life fragment of the fragment}
#'     \item{half_life:}{The half-life mean of the fragment}
#'     \item{HL_SD:}{The half-life standard deviation of the fragment}
#'     \item{HL_SE:}{The half-life standard error of the fragment}
#'     \item{intensity_fragment:}{The intensity_fragment of the fragment}
#'     \item{intensity:}{The relative intensity at time point 0}
#'     \item{intensity_SD:}{The intensity standard deviation of the fragment}
#'     \item{intensity_SE:}{The intensity standard error of the fragment}
#'     \item{velocity:}{The velocity value of the respective delay fragment}
#'     }
#'   }
#'   \item{event_df:}{all information regarding events:
#'   \describe{
#'     \item{event:}{}
#'     \item{p_value:}{}
#'     \item{p_adjusted:}{}
#'     \item{FC_HL:}{Fold change of half-life}
#'     \item{FC_intensity:}{Fold change of intensity}
#'     \item{FC_HL_adapted:}{Fold change of half-life/ fold change of intensity,
#'     position of the half-life fragment is adapted to intensity fragment}
#'     \item{FC_HL_FC_intensity:}{Fold change of half-life/ fold change of
#'     intensity}
#'     \item{event_position:}{}
#'     \item{velocity_ratio:}{}
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{segment_1:}{}
#'     \item{segment_2:}{}
#'     \item{event_duration:}{}
#'     \item{gap_fragments:}{}
#'     \item{features:}{}
#'     }
#'   }
#'   \item{events_HL_int_df:}{all information regarding events related to 
#'   half-life and intensity:
#'   \describe{
#'     \item{event:}{}
#'     \item{p_value:}{}
#'     \item{p_adjusted:}{}
#'     \item{FC_HL:}{}
#'     \item{FC_intensity:}{}
#'     \item{FC_HL_adapted:}{Fold change of half-life/ fold change of intensity,
#'     position of the half-life fragment is adapted to intensity fragment}
#'     \item{FC_HL_FC_intensity:}{Fold change of half-life/ fold change of
#'     intensity}
#'     \item{event_position:}{}
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{segment_1:}{}
#'     \item{segment_2:}{}
#'     \item{event_duration:}{}
#'     \item{gap_fragments:}{}
#'     \item{features:}{}
#'     }
#'   }
#'   \item{events_ps_itss_df:}{all information regarding events related to 
#'   pausing sites and iTSS_I:
#'   \describe{
#'     \item{event:}{}
#'     \item{p_value:}{}
#'     \item{p_adjusted:}{}
#'     \item{event_position:}{}
#'     \item{velocity_ratio:}{}
#'     \item{FC_HL_adapted:}{}
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{segment_1:}{}
#'     \item{segment_2:}{}
#'     \item{event_duration:}{}
#'     \item{gap_fragments:}{}
#'     \item{features:}{}
#'     }
#'   }
#'   \item{events_velocity_df:}{all information regarding events related to 
#'   velocity:
#'   \describe{
#'     \item{event:}{}
#'     \item{p_value:}{}
#'     \item{p_adjusted:}{}
#'     \item{event_position:}{}
#'     \item{velocity_ratio:}{}
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{segment_1:}{}
#'     \item{segment_2:}{}
#'     \item{event_duration:}{}
#'     \item{gap_fragments:}{}
#'     \item{features:}{}
#'     }
#'   }
#'   \item{TI_df:}{all information regarding TI:
#'   \describe{
#'     \item{event:}{}
#'     \item{TI_fragment:}{}
#'     \item{TI_termination_factor:}{}
#'     \item{p_value:}{}
#'     \item{p_adjusted:}{}
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{features:}{}
#'     \item{event_position:}{}
#'     \item{position_1:}{}
#'     \item{position_2:}{}
#'     }
#'   }
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(summary_minimal)
#'
"summary_minimal"

#' The result of rifi_summary for Synechocystis 6803 example data
#' A list containing the output from rifi_summary, including the fragment
#' based data frame, bin based data frame, event data frame and the TI
#' dataframe.
#' @format A list of 4 data frames with 3000 rows and 11 variables, 297
#' rows and 11 variables, 486 rows and 18 variables, and 10 rows and 14
#' variables:
#' \describe{
#'   \item{bin_df:}{all information regarding bins:
#'   \describe{
#'     \item{ID:}{}
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{position:}{}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{segment:}{The segment the bin/probe belongs to}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{delay_fragment:}{The delay fragment the bin/probe belongs to}
#'     \item{delay:}{The delay of the bin/probe}
#'     \item{HL_fragment:}{The half-life fragment the bin/probe belongs to}
#'     \item{half_life:}{The half-life of the bin/probe}
#'     \item{intensity_fragment:}{The intensity fragment the bin/probe belongs
#'      to}
#'     \item{intensity:}{The relative intensity at time point 0}
#'     \item{flag:}{The flag of the bin/probe(TI, PDD)}
#'     \item{TI_termination_factor:}{The TI_termination_factor of the bin/probe 
#'     (in case TI is detected)}
#'     }
#'   }
#'   \item{frag_df:}{all information regarding fragments:
#'   \describe{
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{first_position_frg:}{The first position of the fragment on the
#'      genome}
#'     \item{last_position_frg:}{The last position of the fragment on the
#'      genome}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{segment:}{The segment the fragment belongs to}
#'     \item{delay_fragment:}{The delay fragment of the fragment}
#'     \item{HL_fragment:}{The half-life fragment of the fragment}
#'     \item{half_life:}{The half-life mean of the fragment}
#'     \item{HL_SD:}{The half-life standard deviation of the fragment}
#'     \item{HL_SE:}{The half-life standard error of the fragment}
#'     \item{intensity_fragment:}{The intensity_fragment of the fragment}
#'     \item{intensity:}{The relative intensity at time point 0}
#'     \item{intensity_SD:}{The intensity standard deviation of the fragment}
#'     \item{intensity_SE:}{The intensity standard error of the fragment}
#'     \item{velocity:}{The velocity value of the respective delay fragment}
#'     }
#'   }
#'   \item{event_df:}{all information regarding events:
#'   \describe{
#'     \item{event:}{}
#'     \item{p_value:}{}
#'     \item{p_adjusted:}{}
#'     \item{FC_HL:}{Fold change of half-life}
#'     \item{FC_intensity:}{Fold change of intensity}
#'     \item{FC_HL_adapted:}{Fold change of half-life/ fold change of intensity,
#'     position of the half-life fragment is adapted to intensity fragment}
#'     \item{FC_HL_FC_intensity:}{Fold change of half-life/ fold change of
#'     intensity}
#'     \item{event_position:}{}
#'     \item{velocity_ratio:}{}
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{segment_1:}{}
#'     \item{segment_2:}{}
#'     \item{event_duration:}{}
#'     \item{gap_fragments:}{}
#'     \item{features:}{}
#'     }
#'   }
#'   \item{events_HL_int_df:}{all information regarding events related to 
#'   half-life and intensity:
#'   \describe{
#'     \item{event:}{}
#'     \item{p_value:}{}
#'     \item{p_adjusted:}{}
#'     \item{FC_HL:}{}
#'     \item{FC_intensity:}{}
#'     \item{FC_HL_adapted:}{Fold change of half-life/ fold change of intensity,
#'     position of the half-life fragment is adapted to intensity fragment}
#'     \item{FC_HL_FC_intensity:}{Fold change of half-life/ fold change of
#'     intensity}
#'     \item{event_position:}{}
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{segment_1:}{}
#'     \item{segment_2:}{}
#'     \item{event_duration:}{}
#'     \item{gap_fragments:}{}
#'     \item{features:}{}
#'     }
#'   }
#'   \item{events_ps_itss_df:}{all information regarding events related to 
#'   pausing sites and iTSS_I:
#'   \describe{
#'     \item{event:}{}
#'     \item{p_value:}{}
#'     \item{p_adjusted:}{}
#'     \item{event_position:}{}
#'     \item{velocity_ratio:}{}
#'     \item{FC_HL_adapted:}{}
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{segment_1:}{}
#'     \item{segment_2:}{}
#'     \item{event_duration:}{}
#'     \item{gap_fragments:}{}
#'     \item{features:}{}
#'     }
#'   }
#'   \item{events_velocity_df:}{all information regarding events related to 
#'   velocity:
#'   \describe{
#'     \item{event:}{}
#'     \item{p_value:}{}
#'     \item{p_adjusted:}{}
#'     \item{event_position:}{}
#'     \item{velocity_ratio:}{}
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{segment_1:}{}
#'     \item{segment_2:}{}
#'     \item{event_duration:}{}
#'     \item{gap_fragments:}{}
#'     \item{features:}{}
#'     }
#'   }
#'   \item{TI_df:}{all information regarding TI:
#'   \describe{
#'     \item{event:}{}
#'     \item{TI_fragment:}{}
#'     \item{TI_termination_factor:}{}
#'     \item{p_value:}{}
#'     \item{p_adjusted:}{}
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{features:}{}
#'     \item{event_position:}{}
#'     \item{position_1:}{}
#'     \item{position_2:}{}
#'     }
#'   }
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(summary_synechocystis_6803)
#'
"summary_synechocystis_6803"
#' 
#' 
#' The result of rifi_wrapper for E.coli example data 
#' A list of SummarizedExperiment containing the output of rifi_wrapper. The 
#' list contains 6 elements of SummarizedExperiment output of rifi_preprocess, 
#' rifi_fit, rifi_penalties, rifi_fragmentation, rifi_stats and rifi_summary. 
#' The plot is generated from rifi_visualization. for more detail, please refer 
#' to each function separately. 
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(wrapper_e_coli)
#'
"wrapper_e_coli"

#' The result of rifi_wrapper for E.coli artificial example.
#' A list of SummarizedExperiment containing the output of rifi_wrapper. The 
#' list contains 6 elements of SummarizedExperiment output of rifi_preprocess, 
#' rifi_fit, rifi_penalties, rifi_fragmentation, rifi_stats and rifi_summary. 
#' The plot is generated from rifi_visualization. for more detail, please refer 
#' to each function separately.
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(wrapper_minimal)
#'
"wrapper_minimal"

#' 
#' The result of rifi_wrapper for summary_synechocystis_6803 example data 
#' A list of SummarizedExperiment containing the output of rifi_wrapper. The 
#' list contains 6 elements of SummarizedExperiment output of rifi_preprocess, 
#' rifi_fit, rifi_penalties, rifi_fragmentation, rifi_stats and rifi_summary. 
#' The plot is generated from rifi_visualization. for more detail, please refer 
#' to each function separately. 
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(wrapper_summary_synechocystis_6803)
#'
"wrapper_summary_synechocystis_6803"

#' The result of event_dataframe for E.coli artificial example.
#' A data frame combining the processed genome annotation and a
#' SummarizedExperiment data from rifi_stats. The dataframe is 
#'
#' @format A list with 2 items:
#' \describe{
#'   \item{region:}{the region from the gff file}
#'   \item{gene:}{the annotated gene name}
#'   \item{locus_tag:}{the annotated locus tag}
#'   \item{strand:}{the strand of the annotation}
#'   \item{TU:}{The overarching transcription unit}
#'   \item{position:}{The bin/probe specific position}
#'   \item{FC_fragment_intensity:}{}
#'   \item{FC_intensity:}{}
#'   \item{p_value_intensity:}{}
#'   \item{FC_fragment_HL:}{}
#'   \item{FC_HL:}{}
#'   \item{p_value_HL:}{}
#'   \item{FC_HL_intensity_fragment:}{}
#'   \item{FC_HL_intensity:}{}
#'   \item{FC_HL_adapted:}{}
#'   \item{p_value_Manova:}{}
#'   \item{synthesis_ratio:}{}
#'   \item{synthesis_ratio_event:}{}
#'   \item{pausing_site:}{}
#'   \item{iTSS_I:}{}
#'   \item{event_ps_itss_p_value_Ttest:}{}
#'   \item{ps_ts_fragment:}{}
#'   \item{event_position:}{}
#'   \item{event_duration:}{}
#'   \item{delay_frg_slope:}{}
#'   \item{p_value_slope:}{}
#'   \item{delay:}{}
#'   \item{half_life:}{}
#'   \item{intensity:}{}
#' }
#' 
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(res_minimal)
#'
"res_minimal"

#' A list corresponding to a gff file for E.coli example data
#' A list containing all necessary information from a gff file for rifi_summary 
#' and visualization. The list is stored as metadata in rifi_stats output.
#'
#' @format A list with 2 items:
#' \describe{
#'   \item{data annotation:}{a data frame with 4452 rows and 6 variables
#'     \describe{
#'       \item{region:}{the region from the gff file}
#'       \item{start:}{the start of the annotation}
#'       \item{end:}{the end of the annotation}
#'       \item{strand:}{the strand of the annotation}
#'       \item{gene:}{the annotated gene name}
#'       \item{locus_tag:}{the annotated locus tag}
#'     }
#'   }
#'   \item{genome length:}{a numeric vector containing the length of the genome}
#' }