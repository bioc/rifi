#' rifi_wrapper: conveniently wraps all functions included on rifi workflow.
#' Wraps the functions: check_input, rifi_preprocess, rifi_fit, rifi_penalties,
#' rifi_fragmentation, rifi_stats, rifi_summary and rifi_visualization.
#'
#' @param inp data frame: the input data frame with correct format.
#' @param cores integer: the number of assigned cores for the task.
#' @param gff path: path to an annotation file in gff format
#' @param bg numeric: threshold over which the last time point has to be to be
#' fitted with the above background mode.
#' @param restr numeric: a parameter that restricts the freedom of the fit to
#' avoid wrong TI-term_factors, ranges from 0 to 0.2
#'
#' @return All intermediate objects
#'
#' @seealso `check_input`
#' @seealso `rifi_preprocess`
#'
#' @export


rifi_wrapper <- function(inp, cores, gff, bg, restr) {
  
  #run preprocess step
  prepro <- rifi_preprocess(
    inp = inp,
    cores = cores,
    bg = bg,
    rm_FLT = T,
    thrsh_check = 10,
    dista = 300,
    run_PDD = T
  )
  
  #fit the data
  probe <- rifi_fit(
    inp = prepro,
    cores = cores,
    viz = T,
    restr = restr
  )
  
  #calculate penalties
  pen <- rifi_penalties(
    inp = probe,
    details = T,
    viz = T,
    top_i = 25,
    cores = cores,
    dpt = 1,
    smpl_min = 10,
    smpl_max = 100
  )
  
  # run fragmentation
  probe_fra <- rifi_fragmentation(
    inp = probe,
    cores = cores
    )
  
  #extract the annotation from gff file
  annot <- gff3_preprocess(gff)
  
  #run statistics
  probe_sta <- rifi_stats(
    inp = probe_fra, 
    dista = 300)
  
  #run summary
  probe_summary <-
    rifi_summary(
      probe = probe_sta,
      data_annotation = annot[[1]]
      )
  
  #run visualization
  rifi_visualization(data = probe_sta,
                     genomeLength = annot[[2]],
                     annot = annot[[1]])
  
  res <-
    list(prepro, probe, pen, probe_fra, probe_sta, probe_summary)
  
  return(res)
}
