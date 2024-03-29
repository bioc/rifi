#' =========================================================================
#' rifi_wrapper
#' -------------------------------------------------------------------------
#' rifi_wrapper wraps conveniently all functions included on rifi workflow
#'
#' rifi_wrapper wraps the functions: 
#' 1. rifi_preprocess
#' 2. rifi_fit
#' 3. rifi_penalties
#' 4. rifi_fragmentation
#' 5. rifi_stats
#' 6. rifi_summary
#' 7. rifi_visualization.
#'
#' @param inp data frame: the input data frame with correct format.
#' @param cores integer: the number of assigned cores for the task.
#' @param path path: path to an annotation file in gff format.
#' @param bg numeric: threshold over which the last time point has to be to be
#' fitted with the above background mode.
#' @param restr numeric: a parameter that restricts the freedom of the fit to
#' avoid wrong TI-term_factors, ranges from 0 to 0.2
#'
#' @return All intermediate objects
#'
#' @seealso `rifi_preprocess`
#' @seealso `rifi_fit`
#' @seealso `rifi_penalties`
#' @seealso `rifi_fragmentation`
#' @seealso `rifi_stats`
#' @seealso `rifi_summary`
#' @seealso `rifi_visualization`
#' 
#' @examples
#' data(example_input_minimal)
#' rifi_wrapper(inp = example_input_minimal, cores = 2, path = 
#' gzfile(system.file("extdata", "gff_e_coli.gff3.gz", package = "rifi")), 
#' bg = 0, restr = 0.01)
#' 
#' @export


rifi_wrapper <- function(inp, cores, path, bg, restr) {
  
  #run preprocess step
  prepro <- rifi_preprocess(
    inp = inp,
    cores = cores,
    bg = bg,
    rm_FLT = TRUE,
    thrsh_check = 10,
    dista = 300,
    run_PDD = TRUE
  )
  
  #fit the data
  probe <- rifi_fit(
    inp = prepro,
    cores = cores,
    viz = TRUE,
    restr = restr
  )
  
  #calculate penalties
  pen <- rifi_penalties(
    inp = probe,
    details = TRUE,
    viz = TRUE,
    top_i = 25,
    cores = cores,
    dpt = 1,
    smpl_min = 10,
    smpl_max = 100
  )
  
  # run fragmentation
  probe_fra <- rifi_fragmentation(
    inp = pen,
    cores = cores
    )

  #run statistics
  probe_sta <- rifi_stats(
    inp = probe_fra, 
    dista = 300,
    path = path
    )
  
  #run summary
  probe_summary <-
    rifi_summary(
      inp = probe_sta,
      data_annotation = metadata(probe_sta)$annot[[1]]
      )
  
  #run visualization
  rifi_visualization(data = probe_sta,
                     genomeLength = metadata(probe_sta)$annot[[2]],
                     annot = metadata(probe_sta)$annot[[1]])
  
  res <-
    list(prepro, probe, pen, probe_fra, probe_sta, probe_summary)
  
  return(res)
}
