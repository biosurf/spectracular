optimal_combi_spectral <- function(k, fluorochromes_spectra, instrument, fridge = NULL, constraints = NULL) {

  # ------------------------------------------------------------
  # Calculate spectral fingerprints
  # ------------------------------------------------------------

  instrument <- validate_instrument(instrument)
  detector_signal <- calculate_detector_response(
    fluorochromes_spectra,
    instrument
  )$signal


  # ------------------------------------------------------------
  # Calculate cosine similarity matrix
  # ------------------------------------------------------------

  cosine <- matrix(
    0,
    nrow = nrow(detector_signal),
    ncol = nrow(detector_signal),
    dimnames = list(rownames(detector_signal), rownames(detector_signal))
  )

  for(i in seq_len(nrow(detector_signal))) {
    for(ii in seq_len(nrow(detector_signal))) {
      x <- detector_signal[i,]
      y <- detector_signal[ii,]
      cosine[i,ii] <- sum(x * y) / sqrt(sum(x^2) * sum(y^2))
    }
  }


  # ------------------------------------------------------------
  # Generate table of unique edges
  # ------------------------------------------------------------

  edge_index <- which(lower.tri(cosine) & !is.na(cosine), arr.ind = TRUE)
  cosine_long <- data.frame(from = rownames(cosine)[edge_index[,1]], to = colnames(cosine)[edge_index[,2]], value = cosine[edge_index])
  cosine_long <- cosine_long[order(cosine_long$value),]


  # ------------------------------------------------------------
  # Add edges until a valid k-clique is found
  # ------------------------------------------------------------

  nodes <- data.frame(node = rownames(cosine), fluorochrome = rownames(cosine))
  panel_graph <- constrain_panel_graph(k, cosine_long, nodes, constraints, fridge)
  edges <- panel_graph$edges
  original <- setNames(panel_graph$nodes$original, panel_graph$nodes$node)
  combi_list <- list()
  thresholds <- if(k == 1) 0 else sort(unique(edges$value))

  for(threshold in thresholds) {
    edgelist <- edges[edges$value <= threshold,1:2]
    cliques <- find_k_cliques(edgelist, panel_graph$nodes$node, k, limit = 2)
    combi_list <- lapply(cliques, function(x) unname(original[x]))

    if(length(combi_list) > 0) break
  }
  if(length(combi_list) == 0) stop("No feasible panel satisfies the requested size and constraints.")


  # ------------------------------------------------------------
  # Resolve ties by total similarity
  # ------------------------------------------------------------

  if(length(combi_list) > 1) {
    best_nodes <- minimum_sum_clique(edges[edges$value <= threshold,,drop = FALSE], panel_graph$nodes$node, k)
    best_combi <- unname(original[best_nodes])
  } else {
    best_combi <- combi_list[[1]]
  }

  return(unname(best_combi))
}
