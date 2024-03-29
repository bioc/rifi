% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gff3_preprocessing.r
\name{gff3_preprocess}
\alias{gff3_preprocess}
\title{\subsection{=========================================================================
gff3_preprocess}{

gff3_preprocess process gff3 file from database for multiple usage
}}
\usage{
gff3_preprocess(path)
}
\arguments{
\item{path}{path: path to the directory containing the gff3 file.}
}
\value{
A list with 2 items:
\describe{
\item{data annotation:}{
\describe{
\item{region:}{String, the region from the gff file}
\item{start:}{Integer, the start of the annotation}
\item{end:}{Integer, the end of the annotation}
\item{strand:}{Boolean, the strand of the annotation}
\item{gene:}{String, the annotated gene name}
\item{locus_tag:}{String, the annotated locus tag}
}
}
\item{genome length:}{a numeric vector containing the length of the genome}
}
}
\description{
gff3_preprocess processes the gff3 file extracting gene names and locus_tag
from all coding regions (CDS), UTRs/ncRNA/asRNA are also extracted if
available.
}
\details{
The resulting dataframe contains region, positions, strand, gene and
locus_tag.
}
\examples{
gff3_preprocess(
path = gzfile(system.file("extdata", "gff_e_coli.gff3.gz", package = "rifi"))
)
}
