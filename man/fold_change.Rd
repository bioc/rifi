% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fold_change.r
\name{fold_change}
\alias{fold_change}
\title{\subsection{=========================================================================
fold_change}{

fold_change sets a fold-change ratio between the neighboring fragments of
Half-life (HL) and intensity
}}
\usage{
fold_change(inp)
}
\arguments{
\item{inp}{SummarizedExperiment: the input data frame with correct format.}
}
\value{
the SummarizedExperiment with the columns regarding statistics:
\describe{
\item{ID:}{The bin/probe specific ID.}
\item{position:}{The bin/probe specific position.}
\item{strand:}{The bin/probe specific strand.}
\item{intensity:}{The relative intensity at time point 0.}
\item{probe_TI:}{An internal value to determine which fitting model
is applied.}
\item{flag:}{Information on which fitting model is applied.}
\item{position_segment:}{The position based segment.}
\item{delay:}{The delay value of the bin/probe.}
\item{half_life:}{The half-life of the bin/probe.}
\item{TI_termination_factor:}{String, the factor of TI fragment.}
\item{delay_fragment:}{The delay fragment the bin belongs to.}
\item{velocity_fragment:}{The velocity value of the respective delay
fragment.}
\item{intercept:}{The vintercept of fit through the respective delay
fragment.}
\item{slope:}{The slope of the fit through the respective delay fragment.}
\item{HL_fragment:}{The half-life fragment the bin belongs to.}
\item{HL_mean_fragment:}{The mean half-life value of the respective
half-life fragment.}
\item{intensity_fragment:}{The intensity fragment the bin belongs to.}
\item{intensity_mean_fragment:}{The mean intensity value of the respective
intensity fragment.}
\item{TU:}{The overarching transcription unit.}
\item{TI_termination_fragment:}{The TI fragment the bin belongs to.}
\item{TI_mean_termination_factor:}{The mean termination factor of the
respective TI fragment.}
\item{seg_ID:}{The combined ID of the fragment.}
\item{pausing_site:}{presence of pausing site indicated by +/-.}
\item{iTSS_I:}{presence of iTSS_I indicated by +/-.}
\item{ps_ts_fragment:}{The fragments involved in pausing site or iTSS_I.}
\item{event_duration:}{Integer, the duration between two delay fragments.}
\item{event_ps_itss_p_value_Ttest:}{p_value of pausing site or iTSS_I.}
\item{p_value_slope:}{Integer, the p_value added to the inp.}
\item{delay_frg_slope:}{Integer, the slope value of the fit through the
respective delay fragment.}
\item{velocity_ratio:}{Integer, the ratio value of velocity from 2 delay
fragments.}
\item{event_position:}{Integer, position of the event added to the input.}
\item{FC_HL:}{Integer, the fold change value of 2 HL fragments.}
\item{FC_fragment_HL:}{String, the fragments corresponding to HL fold
change.}
\item{p_value_HL:}{Integer, the p_value added to the input of 2 HL
fragments.}
\item{FC_intensity:}{Integer, the fold change value of 2 intensity
fragments.}
\item{FC_fragment_intensity:}{String, the fragments corresponding to
intensity fold change.}
\item{p_value_intensity:}{Integer, the p_value added to the input of 2
intensity fragments.}
\item{synthesis_ratio:}{Integer, the value correspomding to synthesis rate.}
\item{synthesis_ratio_event:}{String, the event assigned by synthesis
rate either Termination or iTSS.}
\item{FC_HL_intensity:}{Integer, the value corresponding to HL and
intensity fold change.}
\item{FC_HL_intensity_fragment:}{String, the fragments corresponding to
intensity and HL fold change.}
\item{FC_HL_adapted:}{Integer, the fold change of half-life/ fold change
of intensity,position of the half-life fragment is adapted to intensity
fragment.}
}
}
\description{
fold_change sets fold change on intensity and fold change HL fragments of
two successive fragments. Two intensity fragments could belong to one HL
fragment.
This function sets first the borders using the position and applies the fold
change ratio between the neighboring fragments of HL and those from intensity
log2(intensity frgA/intensity frgB/half-life frgA/half-life frgB). All
grepped fragments are from the same TU excluding outliers.
}
\details{
The function used is:
synthesis_r_Function: assigns events depending on the ratio between HL and
intensity of two consecutive fragments.
intensity(int) = synthesis rate(k)/decay(deg) (steady state),
int1/int2 = k1/deg1*deg2/k2
int1 * (deg1/int2) * deg2 = k1/k2 => synthesis ratio.
In case of synthesis ratio is:
synthesis ratio > 0 -> New start
synthesis ratio < 0 -> Termination
}
\examples{
data(stats_minimal)
fold_change(inp = stats_minimal)

}
