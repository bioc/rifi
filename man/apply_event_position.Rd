% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apply_event_position.r
\name{apply_event_position}
\alias{apply_event_position}
\title{\subsection{=========================================================================
apply_event_position}{

apply_event_position extracts event time duration for pausing site or
iTSS
}}
\usage{
apply_event_position(inp)
}
\arguments{
\item{inp}{SummarizedExperiment: the input data frame with correct
format.}
}
\value{
The SummarizedExperiment with the columns regarding statistics:
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
}
}
\description{
apply_event_position is a short version of apply_Ttest_delay function to
extract event time duration for pausing site or iTSS. Its adds a new column
with the duration.
}
\examples{
data(stats_minimal)
apply_event_position(inp = stats_minimal)

}
