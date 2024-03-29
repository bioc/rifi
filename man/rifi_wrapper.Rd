% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rifi_wrapper.r
\name{rifi_wrapper}
\alias{rifi_wrapper}
\title{\subsection{=========================================================================
rifi_wrapper}{

rifi_wrapper wraps conveniently all functions included on rifi workflow
}}
\usage{
rifi_wrapper(inp, cores, path, bg, restr)
}
\arguments{
\item{inp}{data frame: the input data frame with correct format.}

\item{cores}{integer: the number of assigned cores for the task.}

\item{path}{path: path to an annotation file in gff format.}

\item{bg}{numeric: threshold over which the last time point has to be to be
fitted with the above background mode.}

\item{restr}{numeric: a parameter that restricts the freedom of the fit to
avoid wrong TI-term_factors, ranges from 0 to 0.2}
}
\value{
All intermediate objects
}
\description{
rifi_wrapper wraps the functions:
\enumerate{
\item rifi_preprocess
\item rifi_fit
\item rifi_penalties
\item rifi_fragmentation
\item rifi_stats
\item rifi_summary
\item rifi_visualization.
}
}
\examples{
data(example_input_minimal)
rifi_wrapper(inp = example_input_minimal, cores = 2, path = 
gzfile(system.file("extdata", "gff_e_coli.gff3.gz", package = "rifi")), 
bg = 0, restr = 0.01)

}
\seealso{
\code{rifi_preprocess}

\code{rifi_fit}

\code{rifi_penalties}

\code{rifi_fragmentation}

\code{rifi_stats}

\code{rifi_summary}

\code{rifi_visualization}
}
