# =========================================================================
# gff3_preprocess    Process gff3 file from database for multiple usage
# -------------------------------------------------------------------------
#'
#'
#' gff3_preprocess processes the gff3 file extracting gene names and locus_tag
#' from all coding regions (CDS), UTRs/ncRNA/asRNA are also extracted if
#' available.

#' The resulting dataframe contains region, positions, strand, gene and
#' locus_tag.
#'
#' @param path path: path to the directory containing the gff3 file.
#' 
#' @return A list with 2 items:
#' \describe{
#'   \item{data annotation:}{
#'     \describe{
#'       \item{region:}{String, the region from the gff file}
#'       \item{start:}{Integer, the start of the annotation}
#'       \item{end:}{Integer, the end of the annotation}
#'       \item{strand:}{Boolean, the strand of the annotation}
#'       \item{gene:}{String, the annotated gene name}
#'       \item{locus_tag:}{String, the annotated locus tag}
#'     }
#'   }
#'   \item{genome length:}{a numeric vector containing the length of the genome}
#' }
#'
#' @examples
#' gff3_preprocess(
#' path = gzfile(system.file("extdata", "gff_e_coli.gff3.gz", package = "rifi"))
#' )
#' @export

gff3_preprocess <- function(path) {
  #grep the line containing the genome length
  ge_size <- genomeSize(path)
  #select columns region, start, end, strand and annotation
  inp <- readGFF(path)
  tmp <- inp[, c("type", "start", "end", "strand", "gene", "locus_tag")]
  tmp <-
    tmp[grep("^CDS$|UTR|asRNA|antisense_RNA|ncRNA|^tRNA$",
             tmp$type,
             invert = FALSE), ]
  colnames(tmp)[1] <- "region"
   #replace gene with NA with locus_tag in case NA is not recognized as NA
  tmp[which(tmp$gene == "NA"), "gene"] <-
    tmp[which(tmp$gene == "NA"), "locus_tag"]
  #replace gene with NA with locus_tag
  tmp[is.na(tmp$gene), "gene"] <- tmp[is.na(tmp$gene), "locus_tag"]
  if (length(which(tmp$region %in% "antisense_RNA")) != 0) {
    tmp[which(tmp$region == "antisense_RNA"), "region"] <- "asRNA"
  }
  return(list(tmp, ge_size))
}
