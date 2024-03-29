% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nls2_fit.r
\name{nls2_fit}
\alias{nls2_fit}
\title{\subsection{=========================================================================
nls2_fit}{

nls2_fit estimates decay for each probe or bin
}}
\usage{
nls2_fit(
  inp,
  cores = 1,
  decay = seq(0.01, 0.11, by = 0.02),
  delay = seq(0, 10, by = 0.1),
  k = seq(0.1, 1, 0.2),
  bg = 0.2
)
}
\arguments{
\item{inp}{SummarizedExperiment: the input with correct format.}

\item{cores}{integer: the number of assigned cores for the task.}

\item{decay}{numeric vector: A sequence of starting values for the decay.
Default is seq(.08, 0.11, by=.02)}

\item{delay}{numeric vector: A sequence of starting values for the delay.
Default is seq(0,10, by=.1)}

\item{k}{numeric vector: A sequence of starting values for the synthesis
rate. Default is seq(0.1,1,0.2)}

\item{bg}{numeric vector: A sequence of starting values. Default is 0.2.}
}
\value{
the SummarizedExperiment object: with delay and decay added to the
rowRanges. The full fit data is saved in the metadata as "fit_STD".
\describe{
\item{delay:}{Integer, the delay value of the bin/probe}
\item{half_life:}{Integer, the half-life of the bin/probe}
}
}
\description{
nls2_fit uses nls2 function to fit a probe or bin using intensities of the
time series data from different time point. nls2 uses different starting
values through expand grid and selects the best fit. Different filters could
be applied prior fitting to the model.
}
\details{
To apply nls2_fit function, prior filtration could applied.
\enumerate{
\item generic_filter_BG: filter probes with intensities below background using
threshold. Those probes are filtered.
\item filtration_below_backg: additional functions exclusive to microarrays
could be applied. Its very strict to the background (not recommended in
usual case).
\item filtration_above_backg: selects probes with a very high intensity and
above the background (recommended for special transcripts). Probes are
flagged with "\emph{ABG}".
Those transcripts are usually related to a specific function in bacteria.
This filter selects all probes with the same ID, the mean is applied,
the last time point is selected and compared to the threshold.
}

The model used estimates the delay, decay, intensity of the first time
point (synthesis rate/decay) and the background.
The coefficients are gathered in vectors with the corresponding IDs.
Absence of the fit or a very bad fit are assigned with NA.
In case of probes with very high intensities and above the background,
the model used makes abstinence of background coefficient.
The output of all coefficients is saved in the metadata.
The fits are plotted using the function_plot_fit.r through rifi_fit.
}
\examples{
data(preprocess_minimal)
nls2_fit(inp = preprocess_minimal, cores = 2)

}
