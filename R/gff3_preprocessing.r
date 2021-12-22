#' gff3_preprocess: process gff3 file from database for multiple usage.
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
#' @examples
#' gff3_preprocess(
#' path = gzfile(system.file("extdata", "gff_e_coli.gff3.gz", package = "rifi"))
#' )
#' gff3_preprocess(
#' path = gzfile(system.file("extdata", "gff_synechocystis_6803.gff.gz",
#'package = "rifi")))
#'
#' @export


gff3_preprocess <- function(path) {
  inp <- readLines(path)
  #grep the line containing the genome length
  ge_size <- inp[grep("##sequence-region", inp)]
  #replace "\t" in the end of the line
  ge_size <- gsub("\t", "", ge_size)
  #extract the genome size
  ge_size <- as.numeric(last(unlist(strsplit(ge_size, "\\s+"))))
  #select columns region, start, end, strand and annotation
  inp <-
    read.delim2(path,
                header = FALSE,
                sep = "\t",
                comment.char = "#")
  tmp <- inp[, c(3:5, 7, 9)]
  tmp <- na.omit(tmp)
  tmp <-
    tmp[grep("^CDS$|UTR|asRNA|antisense_RNA|ncRNA|^tRNA$",
             tmp$V3,
             invert = FALSE), ]
  colnames(tmp) <-
    c("region", "start", "end", "strand", "annotation")
  #grep lines with gene annotation
  gene <- str_extract(tmp$annotation, "\\gene=\\w+")
  gene <- gsub("gene=", "", gene)
  #grep lines with locus_tag annotation
  locus_tag <- str_extract(tmp$annotation, "\\locus_tag=\\w+")
  locus_tag <- gsub("locus_tag=", "", locus_tag)
  tmp <- cbind.data.frame(tmp[, -5], gene, locus_tag)
  tmp[, -c(2:3)] <- apply(tmp[, -c(2:3)], 2, as.character)
  tmp[, c(2:3)] <- apply(tmp[, c(2:3)], 2, as.numeric)
  #replace gene with NA with locus_tag in case NA is not recognized as NA
  tmp[which(tmp$gene == "NA"), "gene"] <-
    tmp[which(tmp$gene == "NA"), "locus_tag"]
  #replace gene with NA with locus_tag
  tmp[is.na(tmp$gene), "gene"] <- tmp[is.na(tmp$gene), "locus_tag"]
  tmp <- tmp[!duplicated(tmp), ]
  if (length(which(tmp$region %in% "antisense_RNA")) != 0) {
    tmp[which(tmp$region == "antisense_RNA"), "region"] <- "asRNA"
  }
  return(list(tmp, ge_size))
}
