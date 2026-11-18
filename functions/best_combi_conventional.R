optimal_combi_conventional <- function(k, fluorochromes_spectra, instrument, mer, fridge = NULL, constraints = NULL) {

  # ------------------------------------------------------------
  # Prepare instrument and detector responses
  # ------------------------------------------------------------

  detectors <- validate_instrument(instrument)
  detectors$detector <- instrument_detector_ids(detectors)

  if(k > nrow(detectors)) {
    stop("In conventional mode, k cannot exceed the number of detectors.")
  }

  response <- calculate_detector_response(fluorochromes_spectra, detectors)
  detector_signal <- response$signal

  # A fluorochrome-detector assignment is eligible only if both
  # excitation and emission within the detector window exceed MER.
  eligible <- response$excitation_ratio >= mer &
    response$emission_peak >= mer


  # ------------------------------------------------------------
  # Generate fluorochrome-detector nodes
  # ------------------------------------------------------------

  node_index <- which(eligible, arr.ind = TRUE)

  if(nrow(node_index) < k) stop("Too few eligible detector assignments for the requested panel size.")

  nodes <- data.frame(
    node = paste0("n", seq_len(nrow(node_index))),
    fluorochrome = rownames(eligible)[node_index[,1]],
    detector = colnames(eligible)[node_index[,2]],
    f = node_index[,1],
    d = node_index[,2]
  )


  # ------------------------------------------------------------
  # Generate weighted edges between compatible assignments
  # ------------------------------------------------------------

  spillover <- matrix(
    NA_real_,
    nrow = nrow(nodes),
    ncol = nrow(nodes),
    dimnames = list(nodes$node, nodes$node)
  )

  edge_list <- list()
  n <- 1

  for(a in seq_len(nrow(nodes) - 1)) {
    for(b in (a + 1):nrow(nodes)) {

      # A fluorochrome and detector can each only be selected once.
      if(nodes$fluorochrome[a] == nodes$fluorochrome[b]) next
      if(nodes$detector[a] == nodes$detector[b]) next

      # Spillover in both directions relative to assigned detectors.
      ab <- detector_signal[nodes$f[a], nodes$d[b]] /
        detector_signal[nodes$f[a], nodes$d[a]]

      ba <- detector_signal[nodes$f[b], nodes$d[a]] /
        detector_signal[nodes$f[b], nodes$d[b]]

      value <- max(ab, ba)
      spillover[a,b] <- spillover[b,a] <- value

      edge_list[[n]] <- data.frame(
        from = nodes$node[a],
        to = nodes$node[b],
        value = value
      )

      n <- n + 1
    }
  }

  edges <- if(length(edge_list)) do.call(rbind, edge_list) else data.frame(from = character(), to = character(), value = numeric())
  edges <- edges[order(edges$value),]


  # ------------------------------------------------------------
  # Add edges until a valid k-clique is found
  # ------------------------------------------------------------

  panel_graph <- constrain_panel_graph(k, edges, nodes, constraints, fridge)
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
  # Resolve ties by total spillover
  # ------------------------------------------------------------

  if(length(combi_list) > 1) {
    selected <- minimum_sum_clique(edges[edges$value <= threshold,,drop = FALSE], panel_graph$nodes$node, k)
    best_nodes <- unname(original[selected])
  } else {
    best_nodes <- combi_list[[1]]
  }


  # ------------------------------------------------------------
  # Return fluorochromes named by optimized detector assignment
  # ------------------------------------------------------------

  x <- match(best_nodes, nodes$node)
  x <- x[order(match(nodes$detector[x], detectors$detector))]

  best_combi <- nodes$fluorochrome[x]
  names(best_combi) <- nodes$detector[x]

  return(best_combi)
}
