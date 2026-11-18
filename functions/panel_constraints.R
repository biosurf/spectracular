# Shared structure checks; known candidates are optional for editors/settings files.
validate_panel_constraints <- function(constraints, known = NULL) {
  if(is.null(constraints)) constraints <- list()
  if(!is.list(constraints)) stop("Constraints must be a list.")
  if(is.null(constraints$must_include)) constraints$must_include <- character()
  if(is.null(constraints$cannot_combine)) constraints$cannot_combine <- list()
  if(is.null(constraints$groups)) constraints$groups <- list()
  if(!is.character(constraints$must_include) || anyNA(constraints$must_include) ||
     anyDuplicated(constraints$must_include) || any(!nzchar(constraints$must_include))) {
    stop("Must include must contain distinct fluorochrome names.")
  }
  if(!is.list(constraints$cannot_combine) || !is.list(constraints$groups)) {
    stop("Prohibited pairs and groups must be lists.")
  }
  for(pair in constraints$cannot_combine) {
    if(!is.character(pair) || length(pair) != 2 || anyNA(pair) || anyDuplicated(pair) || any(!nzchar(pair))) {
      stop("Each prohibited pair needs two distinct fluorochromes.")
    }
  }
  for(group in constraints$groups) {
    if(!is.list(group)) stop("Each group must contain fluorochromes and n.")
    if(!is.null(group$rule)) {
      stop("Legacy group rules ('exactly'/'at_least') count panel membership and cannot be converted automatically. Recreate groups as distinct required choices using fluorochromes and n, without rule; other panel positions may select additional group members.")
    }
    if(!is.character(group$fluorochromes) || length(group$fluorochromes) == 0 ||
       anyNA(group$fluorochromes) || anyDuplicated(group$fluorochromes) || any(!nzchar(group$fluorochromes))) {
      stop("Each group needs distinct fluorochrome names.")
    }
    if(length(group$n) != 1 || !is.numeric(group$n) || !is.finite(group$n) ||
       group$n < 0 || group$n != floor(group$n) || group$n > length(group$fluorochromes)) {
      stop("Required choices must be a whole number from zero to the number of group members.")
    }
  }
  members <- c(constraints$must_include, unlist(constraints$cannot_combine), unlist(lapply(constraints$groups, `[[`, "fluorochromes")))
  if(!is.null(known) && !all(members %in% known)) stop("Constraints contain unknown fluorochromes: ", paste(setdiff(members, known), collapse = ", "))
  constraints
}


# Expand only when choices need allocation. Original nodes/weights remain the
# source of scientific scores and, in conventional mode, detector assignments.
constrain_panel_graph <- function(k, edges, nodes, constraints = NULL, fridge = NULL) {
  if(length(fridge)) {
    stop("Legacy fridge constraints used exact group-membership counts and cannot be converted automatically. Recreate mandatory candidates in constraints$must_include and groups as distinct required choices with fluorochromes and n; choices cannot be reused across groups.")
  }
  constraints <- validate_panel_constraints(constraints)
  if(length(k) != 1 || !is.numeric(k) || !is.finite(k) || k < 1 || k != floor(k)) stop("Panel size k must be a positive whole number.")
  required <- length(constraints$must_include) + sum(vapply(constraints$groups, `[[`, numeric(1), "n"))
  if(required > k) stop("Required choices need ", required, " panel positions, exceeding k = ", k, ". Mandatory fluorochromes each consume their own position.")
  available <- unique(nodes$fluorochrome)
  missing <- setdiff(constraints$must_include, available)
  if(length(missing)) stop("Required fluorochromes are unavailable or have no eligible assignment: ", paste(missing, collapse = ", "))
  if(length(available) < k) stop("Too few available fluorochromes with eligible assignments for the requested panel size.")
  slots <- as.list(constraints$must_include)
  for(i in seq_along(constraints$groups)) {
    group <- constraints$groups[[i]]
    members <- intersect(group$fluorochromes, available)
    if(length(setdiff(members, constraints$must_include)) < group$n) {
      stop("Group ", i, " has too few available distinct choices after reserving mandatory fluorochromes for their own positions.")
    }
    slots <- c(slots, rep(list(members), group$n))
  }

  # Remove prohibited pairs before copying nodes, so every occurrence is excluded.
  names(edges) <- c("from", "to", "value")
  edges$from <- as.character(edges$from)
  edges$to <- as.character(edges$to)
  from <- nodes$fluorochrome[match(edges$from, nodes$node)]
  to <- nodes$fluorochrome[match(edges$to, nodes$node)]
  keep <- rep(TRUE, nrow(edges))
  for(pair in constraints$cannot_combine) {
    keep <- keep & !((from == pair[1] & to == pair[2]) | (from == pair[2] & to == pair[1]))
  }
  edges <- edges[keep,,drop = FALSE]
  nodes$original <- nodes$node
  if(!length(slots)) return(list(edges = edges, nodes = nodes))

  slots <- c(slots, rep(list(available), k - required))
  copies <- lapply(seq_along(slots), function(i) {
    x <- nodes[nodes$fluorochrome %in% slots[[i]],,drop = FALSE]
    x$slot <- i
    x$node <- paste0("s", i, "_n", seq_len(nrow(x)))
    x
  })
  nodes <- do.call(rbind, copies)
  edges <- merge(edges, data.frame(from = nodes$original, from_node = nodes$node, from_slot = nodes$slot), by = "from", sort = FALSE)
  edges <- merge(edges, data.frame(to = nodes$original, to_node = nodes$node, to_slot = nodes$slot), by = "to", sort = FALSE)
  edges <- edges[edges$from_slot != edges$to_slot, c("from_node", "to_node", "value")]
  names(edges) <- c("from", "to", "value")
  # Base edges already exclude repeated fluorochromes and (where relevant) detectors.
  list(edges = edges, nodes = nodes)
}
