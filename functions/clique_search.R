# Enumerate every size-k clique using standard R. Bounds only discard branches
# that cannot reach k; all feasible choices remain for the minimum-sum tie-break.
# A finite limit returns the same traversal's prefix, for feasibility/local choices.
find_k_cliques <- function(edges, nodes, k, limit = Inf) {
  if(k == 1) return(as.list(nodes[seq_len(min(length(nodes), limit))]))
  indices <- matrix(match(as.matrix(edges[,1:2]), nodes), ncol = 2)
  adjacent <- matrix(FALSE, length(nodes), length(nodes))
  adjacent[indices] <- TRUE
  adjacent <- adjacent | t(adjacent)

  # Every vertex in a k-clique has at least k-1 neighbors in the retained graph.
  candidates <- which(rowSums(adjacent) >= k - 1)
  repeat {
    keep <- rowSums(adjacent[candidates,candidates,drop = FALSE]) >= k - 1
    if(all(keep)) break
    candidates <- candidates[keep]
  }
  if(length(candidates) < k) return(list())
  candidates <- candidates[order(rowSums(adjacent[candidates,candidates,drop = FALSE]), decreasing = TRUE)]
  found <- list()

  search <- function(chosen, candidates) {
    needed <- k - length(chosen)
    if(length(candidates) < needed) return(FALSE)
    if(needed == 1) {
      for(v in candidates) {
        found[[length(found) + 1L]] <<- nodes[c(chosen, v)]
        if(length(found) >= limit) return(TRUE)
      }
      return(FALSE)
    }

    # Each color class is independent. A prefix ending in color c can contain
    # at most c clique members, giving a separate bound for each branch.
    remaining <- candidates
    ordered <- colors <- integer(length(candidates))
    color <- position <- 0L
    while(length(remaining)) {
      color <- color + 1L
      available <- remaining
      while(length(available)) {
        v <- available[1]
        position <- position + 1L
        ordered[position] <- v
        colors[position] <- color
        remaining <- remaining[remaining != v]
        available <- available[available != v & !adjacent[v,available]]
      }
    }

    # Extend through earlier neighbors in this order, visiting each clique once.
    for(i in rev(seq_along(ordered))) {
      if(colors[i] < needed) break
      v <- ordered[i]
      rest <- ordered[seq_len(i - 1L)]
      if(search(c(chosen, v), rest[adjacent[v,rest]])) return(TRUE)
    }
    FALSE
  }
  search(integer(), candidates)
  found
}
