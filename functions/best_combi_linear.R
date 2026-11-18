optimal_combi_linear <- function(k, fluorochromes_spectra, lasers, fridge = NULL, constraints = NULL) {
  
  # ------------------------------------------------------------
  # Recalculate emission spectra based on laser configuration
  # ------------------------------------------------------------
  
  spectra_emission <- fluorochromes_spectra
  
  for(i in seq_along(spectra_emission)) {
    x <- spectra_emission[[i]]
    x$Emission[is.na(x$Emission) | x$Emission < 0] <- 0
    x$Excitation[is.na(x$Excitation) | x$Excitation < 0] <- 0
    
    # Excitation ratios at instrument lasers
    lasers_er <- approx(x$Wavelength, x$Excitation, xout = lasers, rule = 2)$y
    
    # Scale emission by the most exciting laser
    x$Emission <- x$Emission * max(lasers_er, na.rm = TRUE) / 100
    spectra_emission[[i]] <- x
  }
  
  
  # ------------------------------------------------------------
  # Calculate SSC matrix
  # ------------------------------------------------------------
  
  ssc <- matrix(
    NA_real_,
    nrow = length(spectra_emission),
    ncol = length(spectra_emission),
    dimnames = list(names(spectra_emission), names(spectra_emission))
  )
  
  for(i in seq_along(spectra_emission)) {
    for(ii in seq_along(spectra_emission)) {
      s1 <- spectra_emission[[i]]$Emission
      s2 <- spectra_emission[[ii]]$Emission
      
      if(sum(s1) > 0 && sum(s2) > 0) {
        ssc[i,ii] <- sum(pmin(s1, s2)) / min(sum(s1), sum(s2))
      }
    }
  }
  
  
  # ------------------------------------------------------------
  # Generate table of unique edges
  # ------------------------------------------------------------
  
  edge_index <- which(lower.tri(ssc) & !is.na(ssc), arr.ind = TRUE)
  ssc_long <- data.frame(from = rownames(ssc)[edge_index[,1]], to = colnames(ssc)[edge_index[,2]], value = ssc[edge_index])
  ssc_long <- ssc_long[order(ssc_long$value),]
  
  
  # ------------------------------------------------------------
  # Add edges until a valid k-clique is found
  # ------------------------------------------------------------
  
  nodes <- data.frame(node = rownames(ssc), fluorochrome = rownames(ssc))
  panel_graph <- constrain_panel_graph(k, ssc_long, nodes, constraints, fridge)
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
  # If several solutions have the same maximum SSC,
  # select the one with the smallest total SSC
  # ------------------------------------------------------------
  
  if(length(combi_list) > 1) {
    best_nodes <- minimum_sum_clique(edges[edges$value <= threshold,,drop = FALSE], panel_graph$nodes$node, k)
    best_combi <- unname(original[best_nodes])
  } else {
    best_combi <- combi_list[[1]]
  }
  
  return(unname(best_combi))
}
