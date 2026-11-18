# Run from the repository root: Rscript tests/minimum_sum_clique.R
source("functions/clique_search.R")
source("functions/minimum_sum_clique.R")
source("functions/panel_constraints.R")

# Exhaust all five-vertex graphs, comparing every requested size with all subsets.
# Non-integer and negative weights exercise the bound without relying on positivity.
nodes <- c("A", "with space", "3", "s1_n2", "E")
pairs <- t(combn(seq_along(nodes), 2))
for(mask in 0:(2^nrow(pairs) - 1)) {
  selected <- as.logical(intToBits(mask)[seq_len(nrow(pairs))])
  ij <- pairs[selected,,drop = FALSE]
  edges <- data.frame(from = nodes[ij[,1]], to = nodes[ij[,2]], value = (which(selected) %% 7 - 2) / 3)
  adjacent <- matrix(FALSE, length(nodes), length(nodes))
  adjacent[ij] <- TRUE
  adjacent <- adjacent | t(adjacent)
  weights <- matrix(0, length(nodes), length(nodes))
  weights[ij] <- weights[ij[,2:1,drop = FALSE]] <- edges$value
  for(k in seq_along(nodes)) {
    possible <- Filter(function(p) all(adjacent[p,p,drop = FALSE][lower.tri(adjacent[p,p,drop = FALSE])]), combn(seq_along(nodes), k, simplify = FALSE))
    if(!length(possible)) {
      result <- tryCatch(minimum_sum_clique(edges, nodes, k), error = identity)
      stopifnot(inherits(result, "error"), grepl("No feasible size-k clique", conditionMessage(result)))
      next
    }
    scores <- vapply(possible, function(p) {sub <- weights[p,p,drop = FALSE]; sum(sub[lower.tri(sub)])}, numeric(1))
    panel <- match(minimum_sum_clique(edges, nodes, k), nodes)
    sub <- weights[panel,panel,drop = FALSE]
    stopifnot(length(panel) == k, !anyDuplicated(panel), all(adjacent[panel,panel][lower.tri(sub)]), abs(sum(sub[lower.tri(sub)]) - min(scores)) < 1e-12)
    # Bounded enumeration must be the same prefix as exhaustive enumeration.
    all_cliques <- find_k_cliques(edges, nodes, k)
    for(limit in c(1, 2)) stopifnot(identical(find_k_cliques(edges, nodes, k, limit), head(all_cliques, limit)))
  }
}
cat("All 5,120 graph/size cases matched exhaustive minimum sums; bounded enumeration preserved prefixes.\n")

# Larger graphs, including complete graphs whose maximum clique exceeds k.
set.seed(20261112)
for(density in c(.25, .6, 1)) for(repetition in 1:10) {
  nodes <- paste0("n", 1:9)
  ij <- t(combn(seq_along(nodes), 2))
  ij <- ij[runif(nrow(ij)) < density,,drop = FALSE]
  edges <- data.frame(from = nodes[ij[,1]], to = nodes[ij[,2]], value = runif(nrow(ij), 0, 3))
  for(k in c(2, 3, 5, 8)) {
    panels <- find_k_cliques(edges, nodes, k)
    if(!length(panels)) next
    weights <- matrix(0, 9, 9, dimnames = list(nodes, nodes))
    weights[as.matrix(edges[,1:2])] <- edges$value
    weights <- weights + t(weights)
    totals <- vapply(panels, function(p) {sub <- weights[p,p,drop = FALSE]; sum(sub[lower.tri(sub)])}, numeric(1))
    panel <- minimum_sum_clique(edges, nodes, k)
    sub <- weights[panel,panel,drop = FALSE]
    stopifnot(abs(sum(sub[lower.tri(sub)]) - min(totals)) < 1e-12, identical(panel, minimum_sum_clique(edges[nrow(edges):1,,drop = FALSE], nodes, k)))
  }
}
cat("Random and complete graphs matched exhaustive sums with deterministic results.\n")

# Preserve tiny real differences, zero-weight edges and spillovers above one.
nodes <- LETTERS[1:4]
edges <- data.frame(from = c("A", "A", "B", "A", "B"), to = c("B", "C", "C", "D", "D"), value = c(2, 1, 1, 1, 1 - 5e-13))
stopifnot(identical(minimum_sum_clique(edges, nodes, 3), c("A", "B", "D")))
edges$value <- 0
stopifnot(length(minimum_sum_clique(edges, nodes, 3)) == 3)
edges$value <- Inf
stopifnot(length(minimum_sum_clique(edges, nodes, 3)) == 3)

# Mandatory/group slots and detector assignments remain encoded in the graph.
assignments <- expand.grid(fluorochrome = LETTERS[1:4], detector = paste0("D", 1:4), stringsAsFactors = FALSE)
assignments$node <- paste0("original:", seq_len(nrow(assignments)))
ij <- t(combn(seq_len(nrow(assignments)), 2))
ij <- ij[assignments$fluorochrome[ij[,1]] != assignments$fluorochrome[ij[,2]] & assignments$detector[ij[,1]] != assignments$detector[ij[,2]],]
edges <- data.frame(from = assignments$node[ij[,1]], to = assignments$node[ij[,2]], value = (seq_len(nrow(ij)) %% 11) / 4)
constraints <- list(must_include = "A", cannot_combine = list(c("C", "D")), groups = list(
  list(fluorochromes = c("A", "B", "C"), n = 1), list(fluorochromes = c("B", "C", "D"), n = 1)))
graph <- constrain_panel_graph(3, edges, assignments, constraints)
panels <- find_k_cliques(graph$edges, graph$nodes$node, 3)
weights <- matrix(0, nrow(graph$nodes), nrow(graph$nodes), dimnames = list(graph$nodes$node, graph$nodes$node))
weights[as.matrix(graph$edges[,1:2])] <- graph$edges$value
weights <- weights + t(weights)
totals <- vapply(panels, function(p) {sub <- weights[p,p,drop = FALSE]; sum(sub[lower.tri(sub)])}, numeric(1))
best <- minimum_sum_clique(graph$edges, graph$nodes$node, 3)
chosen <- graph$nodes[match(best, graph$nodes$node),]
sub <- weights[best,best,drop = FALSE]
stopifnot(abs(sum(sub[lower.tri(sub)]) - min(totals)) < 1e-12, !anyDuplicated(chosen$fluorochrome), !anyDuplicated(chosen$detector), !anyDuplicated(chosen$slot), "A" %in% chosen$fluorochrome)
stopifnot(!any(c("igraph", "lpSolve", "Rcpp") %in% loadedNamespaces()))
cat("Slot constraints, detector identities, tiny differences and dependency-free execution passed.\n")

# Compare the integrated optimizers with exhaustive sums on their minimax graphs.
spectra <- readRDS("data/fluorochromes_all.rds")[c("FITC", "PE", "APC", "PerCP", "PE-Cy7", "APC-Cy7")]
cases <- list(NULL, list(must_include = "FITC", cannot_combine = list(c("PE", "PE-Cy7"))),
  list(groups = list(list(fluorochromes = c("FITC", "PE", "APC"), n = 2), list(fluorochromes = c("FITC", "APC", "PerCP"), n = 1))))
for(mode in c("linear", "conventional", "spectral")) {
  backend <- new.env()
  for(file in c("instrument_utils", "panel_constraints", "clique_search", "minimum_sum_clique", paste0("best_combi_", mode))) sys.source(paste0("functions/", file, ".R"), backend)
  backend$enumerate <- backend$find_k_cliques
  backend$captured <- NULL
  backend$find_k_cliques <- function(edges, nodes, k, limit = Inf) {
    panels <- enumerate(edges, nodes, k, limit)
    caller <- parent.frame()
    if(length(panels) && !is.null(caller$panel_graph)) {
      stopifnot(limit == 2, length(panels) <= 2)
      captured <<- list(edges = caller$edges[caller$edges$value <= caller$threshold,,drop = FALSE], nodes = nodes,
        metadata = caller$panel_graph$nodes, panels = enumerate(edges, nodes, k), k = k)
    }
    panels
  }
  environment(backend$find_k_cliques) <- backend
  backend$solve_sum <- backend$minimum_sum_clique
  backend$sum_calls <- 0
  backend$minimum_sum_clique <- function(edges, nodes, k) {
    sum_calls <<- sum_calls + 1
    solve_sum(edges, nodes, k)
  }
  environment(backend$minimum_sum_clique) <- backend
  instrument <- read.delim(if(mode == "spectral") "instruments/cytek_aurora_5laser_64channel.tsv" else "instruments/bd_facscantoII_8color.tsv")
  args <- list(k = 3, fluorochromes_spectra = spectra)
  if(mode == "linear") args$lasers <- c(405, 488, 633) else args$instrument <- instrument
  if(mode == "conventional") args$mer <- 10
  for(constraints in cases) {
    backend$sum_calls <- 0
    do.call(backend[[paste0("optimal_combi_", mode)]], c(args, list(constraints = constraints)))
    graph <- backend$captured
    stopifnot(backend$sum_calls == as.integer(length(graph$panels) > 1))
    weights <- matrix(0, length(graph$nodes), length(graph$nodes), dimnames = list(graph$nodes, graph$nodes))
    weights[as.matrix(graph$edges[,1:2])] <- graph$edges$value
    weights <- weights + t(weights)
    totals <- vapply(graph$panels, function(p) {sub <- weights[p,p,drop = FALSE]; sum(sub[lower.tri(sub)])}, numeric(1))
    best <- minimum_sum_clique(graph$edges, graph$nodes, graph$k)
    sub <- weights[best,best,drop = FALSE]
    selected <- graph$metadata[match(best, graph$metadata$node),]
    stopifnot(abs(sum(sub[lower.tri(sub)]) - min(totals)) < 1e-12, !anyDuplicated(selected$fluorochrome))
    if(mode == "conventional") stopifnot(!anyDuplicated(selected$detector))
    if(!is.null(selected$slot)) stopifnot(!anyDuplicated(selected$slot))
  }
  cat(mode, ": real minimax graphs matched the existing minimum-sum tie-break and retained node identities.\n")
}
