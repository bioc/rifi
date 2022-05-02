annotation_function_df <-
  function(feature, pos, strand, data_annotation) {
    data_annotation[, "region"]  <-
      as.character(data_annotation[, "region"])
    data_annotation <-
      data_annotation[which(data_annotation[, "strand"] %in% strand),]
    positions <- which(as.numeric(pos) >= data_annotation[, "start"] &
                         as.numeric(pos) <= data_annotation[, "end"])
    if (length(positions) == 1) {
      features  <- data_annotation[positions, feature]
    } else{
      features  <-
        paste(data_annotation[positions, feature], collapse = "|")
    }
    return(features)
  }

se <- function(x) {
  sqrt(var(x) / length(x))
}
