annotation_function_event <-
  function(x, input, strand, input_annotation) {
    input_annotation[, "region"]  <-
      as.character(input_annotation[, "region"])
    input <- input[which(input[, "strand"] == strand), ]
    input_annotation <-
      input_annotation[which(input_annotation[, "strand"] == strand), ]
    for (j in seq_len(nrow(input))) {
      pos <-
        which((as.numeric(input$position[j]) >= input_annotation[, "start"]) &
                (as.numeric(input$position[j]) <= input_annotation[, "end"]))
      if (is_empty(pos)) {
        next ()
      }
      if (length(input_annotation[pos, x]) == 1) {
        input[which(as.numeric(input$position[j]) ==
                      as.numeric(input$position)), x] <-
          input_annotation[pos, x]
      } else{
        input[which(as.numeric(input$position[j]) ==
                      as.numeric(input$position)), x] <-
          paste(input_annotation[pos, x], collapse = ";")
      }
    }
    input[, x]
  }

position_function <- function(x, input, argument, column) {
  for (i in seq_len(nrow(input))) {
    if (input[i, x] == "+") {
      spl <- unlist(str_split(input[i, "ps_ts_fragment"], ":"))
      spl.1 <-
        last(input[which(input$delay_fragment == spl[1]), argument])
      spl.2 <-
        input[which(input$delay_fragment == spl[2]), argument][1]
      spl.position <- round(sum(spl.1, spl.2) / 2)
      input[i, column] <- spl.position
    }
  }
  input
}



