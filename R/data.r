#' A list corresponding to a gff file for E.coli example data
#' A list containing all necessary information from a gff file for
#' visualization
#' and final output.
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
#'
#' @source \url{https://github.com/CyanolabFreiburg/Transcriptome_wide_decay}
#'
#' @usage data(annot_g_e_coli)
#'
"annot_g_e_coli"

#' A list corresponding to an artificial gff file for the minimal example data
#' A list containing all necessary information from a gff file for visualization
#' and final output.
#' 
#' @format A list with 2 items:
#' \describe{
#'   \item{data annotation:}{a data frame with 4 rows and 6 variables
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
#' @source \url{https://github.com/CyanolabFreiburg/Transcriptome_wide_decay}
#' @usage data(annot_g_minimal)
#'
"annot_g_minimal"

#' A data frame corresponding to a gff file for Synechocystis 6803 example data
#' A list containing all necessary information from a gff file for visualization
#' and final output.
#'
#' @format A list with 2 items:
#' \describe{
#'   \item{data annotation:}{a data frame with 5853 rows and 6 variables:
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
#'
#' @source \url{https://github.com/CyanolabFreiburg/Transcriptome_wide_decay}
#'
#' @usage data(annot_g_synechocystis_6803)
#'
"annot_g_synechocystis_6803"

#' An example input data frame from E. coli.#'
#' A data set from RNAseq data containing information about the intensities at
#' all time points, as well as position, strand and ID information.
#'
#' @format A data frame with 1878 rows and 14 variables:
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
#'   \item{ID:}{unique IDs}
#'   \item{position:}{genome positions}
#'   \item{strand:}{strand information}
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/Transcriptome_wide_decay}
#'
#' @usage data(example_input_e_coli)
#'
"example_input_e_coli"

#' An artificial example input data frame.
#' A data set containing information about the intensities at all time points,
#' as well as position, strand and ID information.
#' @format A data frame with 12 rows and 14 variables:
#' \describe{
#'   \item{0:}{relative intensities at 0 min}
#'   \item{1:}{relative intensities at 1 min}
#'   \item{2:}{relative intensities at 2 min}
#'   \item{3:}{relative intensities at 3 min}
#'   \item{4:}{relative intensities at 4 min}
#'   \item{5:}{relative intensities at 5 min}
#'   \item{6:}{relative intensities at 6 min}
#'   \item{8:}{relative intensities at 8 min}
#'   \item{10:}{relative intensities at 10 min}
#'   \item{15:}{relative intensities at 15 min}
#'   \item{20:}{relative intensities at 20 min}
#'   \item{ID:}{unique IDs}
#'   \item{position:}{genome positions}
#'   \item{strand:}{strand information}
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/Transcriptome_wide_decay}
#'
#' @usage data(example_input_minimal)
#'
"example_input_minimal"

#' An example input data frame from Synechocystis PCC 6803.#'
#' A data set from microarrays data containing information about the intensities
#' at all time points, as well as position, strand and ID information.
#' 
#' @format A data frame with 3000 rows and 10 variables:#' 
#' \describe{
#'   \item{0:}{relative intensities at 0 min}
#'   \item{2:}{relative intensities at 2 min}
#'   \item{4:}{relative intensities at 4 min}
#'   \item{8:}{relative intensities at 8 min}
#'   \item{16:}{relative intensities at 16 min}
#'   \item{32:}{relative intensities at 32 min}
#'   \item{64:}{relative intensities at 64 min}
#'   \item{ID:}{unique IDs}
#'   \item{position:}{genome positions}
#'   \item{strand:}{strand information}
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/Transcriptome_wide_decay}
#'
#' @usage data(example_input_synechocystis_6803)
#'
"example_input_synechocystis_6803"

#' The result of rifi_fit for E.coli example data
#' A list containing the output from rifi_fit, including the probe_df and both
#' fit objects.
#' 
#' @format A list of 3 data frames with 290 rows and 10 variables, 155 rows
#' and 5 variables, and 135 rows and 9 variables:
#' \describe{
#'   \item{probe_df:}{the probe dataframe:
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
#'
#' @source \url{https://github.com/CyanolabFreiburg/Transcriptome_wide_decay}
#'
#' @usage data(fit_e_coli)
#'
"fit_e_coli"

#' The artificial result of rifi_fit for artificial example data
#' A data frame containing the output from rifi_fit.
#' @format A data frame with 24 rows and 10 variables:
#'   \describe{
#'     \item{ID:}{The bin/probe specific ID}
#'     \item{position:}{The bin/probe specific position}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{intensity:}{The relative intensity at time point 0}
#'     \item{probe_TI:}{An internal value to determine which fitting model
#'     is applied}
#'     \item{flag:}{Information on which fitting model is applied}
#'     \item{position_segment:}{The position based segment}
#'     \item{delay:}{The delay value of the bin/probe}
#'     \item{half_life:}{The half-life of the bin/probe}
#'     \item{TI_termination_factor:}{The termination factor of the bin/probe}
#'     }
#' @source \url{https://github.com/CyanolabFreiburg/Transcriptome_wide_decay}
#'
#' @usage data(fit_minimal)
#'
"fit_minimal"

#' The result of rifi_fit for Synechocystis 6803 example data
#' A list containing the output from rifi_fit, including the probe_df and both
#' fit objects.
#' @format A list of 3 data frames with 3000 rows and 10 variables, 2811 rows
#' and 5 variables, and 189 rows and 9 variables:
#' \describe{
#'   \item{probe_df:}{the probe dataframe:
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
#' @source \url{https://github.com/CyanolabFreiburg/Transcriptome_wide_decay}
#'
#' @usage data(fit_synechocystis_6803)
#'
"fit_synechocystis_6803"

#' The result of rifi_fragmentation for E.coli example data
#' A data frame containing the output from rifi_fragmentation, the probe_df
#' @format A data frame with 290 rows and 22 variables:
#' \describe{
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{position:}{The bin/probe specific position}
#'   \item{strand:}{The bin/probe specific strand}
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
#' @source \url{https://github.com/CyanolabFreiburg/Transcriptome_wide_decay}
#'
#' @usage data(fragmentation_e_coli)
#'
"fragmentation_e_coli"

#' The result of rifi_fragmentation for artificial example data
#' A data frame containing the output from rifi_fragmentation, the probe_df
#' @format A data frame with 24 rows and 22 variables:
#' \describe{
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{position:}{The bin/probe specific position}
#'   \item{strand:}{The bin/probe specific strand}
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
#' }
#' 
#' @source \url{https://github.com/CyanolabFreiburg/Transcriptome_wide_decay}
#'
#' @usage data(fragmentation_minimal)
#'
"fragmentation_minimal"

#' The result of rifi_fragmentation for Synechocystis 6803 example data
#' A data frame containing the output from rifi_fragmentation, the probe_df
#' @format A data frame with 3000 rows and 22 variables:
#' \describe{
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{position:}{The bin/probe specific position}
#'   \item{strand:}{The bin/probe specific strand}
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
#' @source \url{https://github.com/CyanolabFreiburg/Transcriptome_wide_decay}
#'
#' @usage data(fragmentation_synechocystis_6803)
#'
"fragmentation_synechocystis_6803"

#' The result of rifi_penalties for E.coli example data
#' A list containing the output from rifi_penalties, including the logbook
#' and the four penalty objects.
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
#' @source \url{https://github.com/CyanolabFreiburg/Transcriptome_wide_decay}
#'
#' @usage data(penalties_e_coli)
#'
"penalties_e_coli"

#' The result of rifi_penalties for artificial example data
#' A vector containing the output from rifi_penalties, the logbook.
#' @format A numeric vector of length 8:
#' \describe{
#'   \item{logbook:}{The logbook vector containing all penalty information}
#' }
#' @source \url{https://github.com/CyanolabFreiburg/Transcriptome_wide_decay}
#'
#' @usage data(penalties_minimal)
#'
"penalties_minimal"

#' The result of rifi_penalties for Synechocystis 6803 example data
#' A list containing the output from rifi_penalties, including the logbook
#' and the four penalty objects.
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
#' @source \url{https://github.com/CyanolabFreiburg/Transcriptome_wide_decay}
#'
#' @usage data(penalties_synechocystis_6803)
#'
"penalties_synechocystis_6803"

#' The result of rifi_preprocess for E.coli example data
#' A list containing the output from rifi_preprocess, including the probe_df
#' and the modified input_df.
#' @format A list of 2 data frames with 290 rows and 7 variables, and 870 rows
#' and 15 variables:
#' \describe{
#'   \item{probe_df:}{the probe dataframe:
#'   \describe{
#'     \item{ID:}{The bin/probe specific ID}
#'     \item{position:}{The bin/probe specific position}
#'     \item{strand:}{The bin/probe specific strand}
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
#'     \item{strand:}{strand information}
#'     \item{filtration:}{indicator wether the replicate is filtered or not}
#'     }
#'   }
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/Transcriptome_wide_decay}
#'
#' @usage data(preprocess_e_coli)
#'
"preprocess_e_coli"

#' The result of rifi_preprocess for artificial example data
#' A list containing the output from rifi_preprocess, including the probe_df
#' and the modified input_df
#' @format A list of 2 data frames with 4 rows and 7 variables, and 12 rows
#' and 15 variables:
#' \describe{
#'   \item{probe_df:}{the probe dataframe:
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
#'     \item{1:}{relative intensities at 1 min}
#'     \item{2:}{relative intensities at 2 min}
#'     \item{3:}{relative intensities at 3 min}
#'     \item{4:}{relative intensities at 4 min}
#'     \item{5:}{relative intensities at 5 min}
#'     \item{6:}{relative intensities at 6 min}
#'     \item{8:}{relative intensities at 8 min}
#'     \item{10:}{relative intensities at 10 min}
#'     \item{15:}{relative intensities at 15 min}
#'     \item{20:}{relative intensities at 20 min}
#'     \item{ID:}{unique IDs}
#'     \item{position:}{genome positions}
#'     \item{strand:}{strand information}
#'     \item{filtration:}{indicator wether the replicate is filtered or not}
#'     }
#'   }
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/Transcriptome_wide_decay}
#'
#' @usage data(preprocess_minimal)
#'
"preprocess_minimal"

#' The result of rifi_preprocess for Synechocystis 6803 example data
#' A list containing the output from rifi_preprocess, including the probe_df
#' and the modified input_df.
#' @format A list of 2 data frames with 3000 rows and 7 variables, and 3000
#' rows and 11 variables:
#' \describe{
#'   \item{probe_df:}{the probe dataframe:
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
#'     \item{strand:}{strand information}
#'     \item{filtration:}{indicator wether the replicate is filtered or not}
#'     }
#'   }
#' }
#' 
#' @source \url{https://github.com/CyanolabFreiburg/Transcriptome_wide_decay}
#'
#' @usage data(preprocess_synechocystis_6803)
#'
"preprocess_synechocystis_6803"

#' The result of rifi_stats for E.coli example data
#' A data frame containing the output from rifi_stats, the probe_df
#' @format A data frame with 290 rows and 45 variables:
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
#' @source \url{https://github.com/CyanolabFreiburg/Transcriptome_wide_decay}
#'
#' @usage data(stats_e_coli)
#'
"stats_e_coli"

#' The result of rifi_stats for artificial example data
#' A data frame containing the output from rifi_stats, the probe_df
#' @format A data frame with 24 rows and 45 variables:
#' \describe{
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{position:}{The bin/probe specific position}
#'   \item{strand:}{The bin/probe specific strand}
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
#' @source \url{https://github.com/CyanolabFreiburg/Transcriptome_wide_decay}
#'
#' @usage data(stats_minimal)
#'
"stats_minimal"

#' The result of rifi_stats for Synechocystis 6803 example data
#' A data frame containing the output from rifi_stats, the probe_df
#' @format A data frame with 3000 rows and 44 variables:
#' \describe{
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{position:}{The bin/probe specific position}
#'   \item{strand:}{The bin/probe specific strand}
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
#' @source \url{https://github.com/CyanolabFreiburg/Transcriptome_wide_decay}
#'
#' @usage data(stats_synechocystis_6803)
#'
"stats_synechocystis_6803"

#' The result of rifi_summary for E.coli example data
#' A list containing the output from rifi_summary, including the fragment
#' based data frame, bin based data frame, event data frame and the TI
#' dataframe
#' @format A list of 4 data frames with 290 rows and 11 variables, 36 rows
#' and 11 variables, 57 rows and 18 variables, and 8 rows and 14 variables:
#' \describe{
#'   \item{bin_df:}{all information regarding bins:
#'   \describe{
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{delay_fragment:}{The delay fragment the bin belongs to}
#'     \item{HL_fragment:}{The half-life fragment the bin belongs to}
#'     \item{half_life:}{The half-life of the bin/probe}
#'     \item{intensity_fragment:}{The intensity fragment the bin belongs to}
#'     \item{intensity:}{The relative intensity at time point 0}
#'     \item{velocity:}{The velocity value of the bin}
#'     }
#'   }
#'   \item{frag_df:}{all information regarding fragments:
#'   \describe{
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{delay_fragment:}{The delay fragment the bin belongs to}
#'     \item{HL_fragment:}{The half-life fragment the bin belongs to}
#'     \item{half_life:}{The half-life of the fragment}
#'     \item{intensity_fragment:}{The intensity fragment the bin belongs to}
#'     \item{intensity:}{The relative intensity at time point 0}
#'     \item{velocity:}{The velocity value of the respective delay fragment}
#'     }
#'   }
#'   \item{event_df:}{all information regarding events:
#'   \describe{
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
#' @source \url{https://github.com/CyanolabFreiburg/Transcriptome_wide_decay}
#'
#' @usage data(summary_e_coli)
#'
"summary_e_coli"
#' The result of rifi_summary for artificial example data
#' A list containing the output from rifi_summary, including the fragment
#' based data frame, bin based data frame, event data frame and the TI
#' dataframe.
#'
#' @format A list of 4 data frames with 290 rows and 11 variables, 36 rows
#' and 11 variables, 57 rows and 18 variables, and 8 rows and 14 variables:
#' \describe{
#'   \item{bin_df:}{all information regarding bins:
#'   \describe{
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{delay_fragment:}{The delay fragment the bin belongs to}
#'     \item{HL_fragment:}{The half-life fragment the bin belongs to}
#'     \item{half_life:}{The half-life of the bin/probe}
#'     \item{intensity_fragment:}{The intensity fragment the bin belongs to}
#'     \item{intensity:}{The relative intensity at time point 0}
#'     \item{velocity:}{The velocity value of the bin}
#'     }
#'   }
#'   \item{frag_df:}{all information regarding fragments:
#'   \describe{
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{delay_fragment:}{The delay fragment the bin belongs to}
#'     \item{HL_fragment:}{The half-life fragment the bin belongs to}
#'     \item{half_life:}{The half-life of the fragment}
#'     \item{intensity_fragment:}{The intensity fragment the bin belongs to}
#'     \item{intensity:}{The relative intensity at time point 0}
#'     \item{velocity:}{The velocity value of the respective delay fragment}
#'     }
#'   }
#'   \item{event_df:}{all information regarding events:
#'   \describe{
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
#' @source \url{https://github.com/CyanolabFreiburg/Transcriptome_wide_decay}
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
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{delay_fragment:}{The delay fragment the bin belongs to}
#'     \item{HL_fragment:}{The half-life fragment the bin belongs to}
#'     \item{half_life:}{The half-life of the bin/probe}
#'     \item{intensity_fragment:}{The intensity fragment the bin belongs to}
#'     \item{intensity:}{The relative intensity at time point 0}
#'     \item{velocity:}{The velocity value of the bin}
#'     }
#'   }
#'   \item{frag_df:}{all information regarding fragments:
#'   \describe{
#'     \item{feature_type:}{}
#'     \item{gene:}{}
#'     \item{locus_tag:}{}
#'     \item{strand:}{The bin/probe specific strand}
#'     \item{TU:}{The overarching transcription unit}
#'     \item{delay_fragment:}{The delay fragment the bin belongs to}
#'     \item{HL_fragment:}{The half-life fragment the bin belongs to}
#'     \item{half_life:}{The half-life of the fragment}
#'     \item{intensity_fragment:}{The intensity fragment the bin belongs to}
#'     \item{intensity:}{The relative intensity at time point 0}
#'     \item{velocity:}{The velocity value of the respective delay fragment}
#'     }
#'   }
#'   \item{event_df:}{all information regarding events:
#'   \describe{
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
#' @source \url{https://github.com/CyanolabFreiburg/Transcriptome_wide_decay}
#'
#' @usage data(summary_synechocystis_6803)
#'
"summary_synechocystis_6803"
