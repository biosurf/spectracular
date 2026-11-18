# Minimum-sum size-k clique in an already thresholded graph. Edges have columns
# from, to and value; nodes are the original slot/assignment IDs, never collapsed.
# Source functions/clique_search.R first. Returns node IDs in input order.
minimum_sum_clique <- function(edges, nodes, k) {
  first <- find_k_cliques(edges, nodes, k, limit = 1)
  if(!length(first)) stop("No feasible size-k clique in the supplied graph.")
  best <- match(first[[1]], nodes)
  if(k == 1 || length(unique(edges$value)) == 1) return(nodes[sort(best)])

  indices <- matrix(match(as.matrix(edges[,1:2]), nodes), ncol = 2)
  adjacent <- matrix(FALSE, length(nodes), length(nodes))
  adjacent[indices] <- TRUE
  adjacent <- adjacent | t(adjacent)
  weights <- matrix(0, length(nodes), length(nodes))
  weights[indices] <- weights[indices[,2:1,drop = FALSE]] <- edges$value
  sub <- weights[best,best,drop = FALSE]
  best_sum <- sum(sub[lower.tri(sub)])
  if(best_sum == 0 && all(edges$value >= 0)) return(nodes[sort(best)])
  # The margin protects lower-bound comparisons; objective values are not rounded.
  roundoff <- 16 * .Machine$double.eps * max(1, abs(edges$value)) * k * (k - 1)

  active <- seq_along(nodes)
  repeat {
    keep <- rowSums(adjacent[active,active,drop = FALSE]) >= k - 1
    if(all(keep)) break
    active <- active[keep]
  }

  # Complement components have every cross-component edge. When their maximum
  # clique sizes sum to k, every solution chooses a maximum clique in each.
  groups <- list()
  unseen <- active
  while(length(unseen)) {
    group <- unseen[1]
    unseen <- unseen[-1]
    position <- 1
    while(position <= length(group)) {
      next_nodes <- unseen[!adjacent[group[position],unseen]]
      group <- c(group, next_nodes)
      unseen <- setdiff(unseen, next_nodes)
      position <- position + 1
    }
    groups[[length(groups) + 1L]] <- group
  }
  compact <- length(groups) > 1 && all(is.finite(edges$value))
  if(compact) {
    local_edges <- lapply(groups, function(g) {
      ij <- which(lower.tri(adjacent[g,g,drop = FALSE]) & adjacent[g,g,drop = FALSE], arr.ind = TRUE)
      data.frame(from = g[ij[,1]], to = g[ij[,2]])
    })
    sizes <- integer(length(groups))
    for(i in seq_along(groups)) {
      lo <- 1
      hi <- min(k, length(groups[[i]]))
      while(lo < hi) {
        mid <- (lo + hi + 1) %/% 2
        if(length(find_k_cliques(local_edges[[i]], groups[[i]], mid, limit = 1))) lo <- mid else hi <- mid - 1
      }
      sizes[i] <- lo
    }
    compact <- sum(sizes) == k
  }
  if(compact) {
    # These are representation budgets, not search limits: larger local families
    # use the exact graph search below, without materializing all complete panels.
    states <- list()
    for(i in seq_along(groups)) {
      states[[i]] <- find_k_cliques(local_edges[[i]], groups[[i]], sizes[i], limit = 100001)
      if(length(states[[i]]) > 100000) {compact <- FALSE; break}
    }
    if(compact) {
      counts <- lengths(states)
      pairs <- outer(counts[counts > 1], counts[counts > 1])
      compact <- sum(pairs[upper.tri(pairs)]) <= 1000000
    }
  }

  if(compact) {
    fixed <- unlist(states[counts == 1], use.names = FALSE)
    choices <- states[counts > 1]
    if(!length(choices)) return(nodes[sort(fixed)])
    unary <- lapply(choices, function(options) vapply(options, function(p) {
      sub <- weights[p,p,drop = FALSE]
      sum(sub[lower.tri(sub)]) + sum(weights[p,fixed,drop = FALSE])
    }, numeric(1)))
    membership <- lapply(choices, function(options) t(vapply(options, function(p) as.numeric(seq_along(nodes) %in% p), numeric(length(nodes)))))
    m <- length(choices)
    pair <- lapply(seq_len(m), function(i) vector("list", m))
    if(m > 1) for(i in seq_len(m - 1)) for(j in (i + 1):m) {
      pair[[i]][[j]] <- membership[[i]] %*% weights %*% t(membership[[j]])
      pair[[j]][[i]] <- t(pair[[i]][[j]])
    }
    sub <- weights[fixed,fixed,drop = FALSE]
    constant <- sum(sub[lower.tri(sub)])

    # Visit influential choices first and seed a feasible incumbent by coordinate
    # descent. All cross-component weights remain in both the score and bound.
    impact <- vapply(seq_len(m), function(i) {
      cost <- unary[[i]]
      for(j in setdiff(seq_len(m), i)) cost <- cost + rowMeans(pair[[i]][[j]])
      diff(range(cost))
    }, numeric(1))
    ordering <- order(impact, decreasing = TRUE)
    choices <- choices[ordering]
    unary <- unary[ordering]
    pair <- lapply(ordering, function(i) pair[[i]][ordering])
    selected <- vapply(unary, which.min, integer(1))
    repeat {
      changed <- FALSE
      for(i in seq_len(m)) {
        cost <- unary[[i]]
        for(j in setdiff(seq_len(m), i)) cost <- cost + pair[[i]][[j]][,selected[j]]
        option <- which.min(cost)
        if(cost[option] < cost[selected[i]] - roundoff) {selected[i] <- option; changed <- TRUE}
      }
      if(!changed) break
    }
    candidate <- c(fixed, unlist(lapply(seq_len(m), function(i) choices[[i]][[selected[i]]]), use.names = FALSE))
    sub <- weights[candidate,candidate,drop = FALSE]
    value <- sum(sub[lower.tri(sub)])
    if(value < best_sum) {best <- candidate; best_sum <- value}

    # Half of each remaining pair cost is bounded at each endpoint. Minimizing
    # the partner independently can only underestimate any complete solution.
    potential <- lapply(seq_len(m), function(depth) lapply(seq_len(m), function(i) {
      cost <- numeric(length(choices[[i]]))
      if(i >= depth) for(j in setdiff(depth:m, i)) cost <- cost + .5 * apply(pair[[i]][[j]], 1, min)
      cost
    }))
    search_choices <- function(depth, costs, total, selected) {
      bound <- total + sum(vapply(depth:m, function(i) min(costs[[i]] + potential[[depth]][[i]]), numeric(1)))
      if(bound > best_sum + roundoff) return(invisible(NULL))
      if(depth == m) {
        selected <- c(selected, which.min(costs[[depth]]))
        candidate <- c(fixed, unlist(lapply(seq_len(m), function(i) choices[[i]][[selected[i]]]), use.names = FALSE))
        sub <- weights[candidate,candidate,drop = FALSE]
        value <- sum(sub[lower.tri(sub)])
        if(value < best_sum) {best <<- candidate; best_sum <<- value}
        return(invisible(NULL))
      }
      for(option in order(costs[[depth]] + potential[[depth]][[depth]])) {
        next_costs <- costs
        for(j in (depth + 1):m) next_costs[[j]] <- costs[[j]] + pair[[depth]][[j]][option,]
        search_choices(depth + 1, next_costs, total + costs[[depth]][option], c(selected, option))
      }
      invisible(NULL)
    }
    search_choices(1, unary, constant, integer())
    return(nodes[sort(best)])
  }

  # General weighted-clique search also covers larger-than-k cliques and graphs
  # without a compact component representation. Only the best panel is retained.
  search <- function(chosen, candidates, costs, total) {
    needed <- k - length(chosen)
    repeat {
      if(length(candidates) < needed) return(invisible(NULL))
      keep <- rowSums(adjacent[candidates,candidates,drop = FALSE]) >= needed - 1
      if(all(keep)) break
      candidates <- candidates[keep]
      costs <- costs[keep]
    }
    if(needed == 1) {
      candidate <- c(chosen, candidates[which.min(costs)])
      sub <- weights[candidate,candidate,drop = FALSE]
      value <- sum(sub[lower.tri(sub)])
      if(value < best_sum) {best <<- candidate; best_sum <<- value}
      return(invisible(NULL))
    }
    # Each future vertex needs needed-1 neighbours; count each future edge half
    # at either endpoint, then optimistically choose the cheapest needed vertices.
    lower <- costs + .5 * vapply(candidates, function(v) {
      sum(sort(weights[v,candidates[adjacent[v,candidates]]])[seq_len(needed - 1)])
    }, numeric(1))
    if(total + sum(sort(lower)[seq_len(needed)]) > best_sum + roundoff) return(invisible(NULL))

    remaining <- candidates[order(lower, decreasing = TRUE)]
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
    for(i in rev(seq_along(ordered))) {
      if(colors[i] < needed) break
      v <- ordered[i]
      rest <- ordered[seq_len(i - 1L)]
      rest <- rest[adjacent[v,rest]]
      search(c(chosen, v), rest, costs[match(rest, candidates)] + weights[v,rest], total + costs[match(v, candidates)])
    }
    invisible(NULL)
  }
  search(integer(), active, numeric(length(active)), 0)
  nodes[sort(best)]
}
