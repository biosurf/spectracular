# Run from the repository root: Rscript tests/clique_search.R
source("functions/clique_search.R")

# Exhaust all 1,024 graphs on five labelled vertices and every requested size.
# This covers isolated vertices, no edges, overlapping maximal cliques and ties.
nodes <- c("A", "with space", "3", "s1_n2", "E")
pairs <- t(combn(seq_along(nodes), 2))
for(mask in 0:(2^nrow(pairs) - 1)) {
  selected <- as.logical(intToBits(mask)[seq_len(nrow(pairs))])
  edge_index <- pairs[selected,,drop = FALSE]
  edges <- data.frame(from = nodes[edge_index[,1]], to = nodes[edge_index[,2]])
  adjacent <- matrix(FALSE, length(nodes), length(nodes))
  adjacent[edge_index] <- TRUE
  adjacent <- adjacent | t(adjacent)
  for(k in seq_along(nodes)) {
    possible <- combn(seq_along(nodes), k, simplify = FALSE)
    expected <- Filter(function(x) all(adjacent[x,x,drop = FALSE][lower.tri(adjacent[x,x,drop = FALSE])]), possible)
    actual <- lapply(find_k_cliques(edges, nodes, k), function(x) sort(match(x, nodes)))
    expected <- sort(vapply(expected, paste, character(1), collapse = ","))
    actual <- vapply(actual, paste, character(1), collapse = ",")
    stopifnot(!anyDuplicated(actual), identical(sort(actual), expected))
  }
}
stopifnot(length(find_k_cliques(edges, nodes, 6)) == 0)
cat("All 5,120 graph/size cases matched exhaustive enumeration without duplicates\n")

# Larger random graphs exercise the bounds at sizes beyond the exhaustive set.
set.seed(431)
nodes <- paste0("node", 1:9)
pairs <- t(combn(seq_along(nodes), 2))
for(density in c(0, .2, .5, .8, 1)) for(repetition in 1:20) {
  edge_index <- pairs[runif(nrow(pairs)) < density,,drop = FALSE]
  edges <- data.frame(from = nodes[edge_index[,1]], to = nodes[edge_index[,2]])
  adjacent <- matrix(FALSE, length(nodes), length(nodes))
  adjacent[edge_index] <- TRUE
  adjacent <- adjacent | t(adjacent)
  for(k in c(2, 3, 5, 8, 9)) {
    possible <- combn(seq_along(nodes), k, simplify = FALSE)
    expected <- Filter(function(x) all(adjacent[x,x,drop = FALSE][lower.tri(adjacent[x,x,drop = FALSE])]), possible)
    actual <- lapply(find_k_cliques(edges, nodes, k), function(x) sort(match(x, nodes)))
    expected <- sort(vapply(expected, paste, character(1), collapse = ","))
    actual <- vapply(actual, paste, character(1), collapse = ",")
    stopifnot(!anyDuplicated(actual), identical(sort(actual), expected))
  }
}
cat("All 500 larger random graph/size cases matched exhaustive enumeration\n")
