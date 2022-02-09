genomeSize <- function(path){
    inp <- readLines(path)
    #grep the line containing the genome length
    ge_size <- inp[grep("##sequence-region", inp)]
    #replace "\t" in the end of the line
    ge_size <- gsub("\t", "", ge_size)
    #extract the genome size
    ge_size <- as.numeric(last(unlist(strsplit(ge_size, "\\s+"))))
    return(ge_size)
}






