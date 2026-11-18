# Run from the repository root: Rscript tests/panel_constraints.R
# Exhaustive enumeration is independent of the graph and production constraint helpers.
source("functions/instrument_utils.R")
source("functions/panel_constraints.R")
source("functions/clique_search.R")
source("functions/minimum_sum_clique.R")
source("functions/excitation_ratio_exclusion.R")
source("functions/wavelength_to_rgb.R")
for(mode in c("linear", "conventional", "spectral")) {
  source(paste0("functions/best_combi_", mode, ".R"))
  source(paste0("functions/plot_", mode, ".R"))
}
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
grDevices::pdf(file = NULL)

enumerate_panels <- function(mode, spectra, lasers, instrument, k, mer) {
  fluoros <- names(spectra)
  if(mode == "linear") {
    emission <- sapply(spectra, function(x) x$Emission * max(approx(x$Wavelength, x$Excitation, xout = lasers, rule = 2)$y) / 100)
    weights <- outer(seq_along(fluoros), seq_along(fluoros), Vectorize(function(i, j) {
      sum(pmin(emission[,i], emission[,j])) / min(sum(emission[,i]), sum(emission[,j]))
    }))
  } else {
    response <- calculate_detector_response(spectra, instrument)
    signal <- response$signal
    if(mode == "spectral") {
      normalized <- signal / sqrt(rowSums(signal^2))
      weights <- tcrossprod(normalized)
    }
  }
  nodes <- data.frame(f = seq_along(fluoros))
  if(mode == "conventional") {
    nodes <- expand.grid(f = seq_along(fluoros), d = seq_len(nrow(instrument)))
    eligible <- response$excitation_ratio >= mer & response$emission_peak >= mer
    nodes <- nodes[eligible[as.matrix(nodes)],]
  }
  panels <- list()
  scores <- list()
  for(index in combn(seq_len(nrow(nodes)), k, simplify = FALSE)) {
    x <- nodes[index,,drop = FALSE]
    if(anyDuplicated(x$f)) next
    panel <- fluoros[x$f]
    if(mode == "conventional") {
      if(anyDuplicated(x$d)) next
      cross_signal <- signal[x$f,x$d]
      directional <- cross_signal / diag(cross_signal)
      pairwise <- pmax(directional, t(directional))
      names(panel) <- colnames(signal)[x$d]
    } else {
      pairwise <- weights[x$f,x$f]
    }
    values <- pairwise[lower.tri(pairwise)]
    panels[[length(panels) + 1L]] <- panel
    scores[[length(scores) + 1L]] <- c(max(values), sum(values))
  }
  list(panels = panels, scores = do.call(rbind, scores))
}

# Allocate requirements to different selected candidates without using the graph.
# This recursive matching oracle also catches shortages across overlapping groups.
can_allocate <- function(choices, remaining) {
  if(!length(choices)) return(TRUE)
  for(candidate in intersect(choices[[1]], remaining)) {
    if(can_allocate(choices[-1], setdiff(remaining, candidate))) return(TRUE)
  }
  FALSE
}
allowed_panel <- function(panel, constraints) {
  for(pair in constraints$cannot_combine) if(all(pair %in% panel)) return(FALSE)
  choices <- as.list(constraints$must_include)
  for(group in constraints$groups) choices <- c(choices, rep(list(group$fluorochromes), group$n))
  can_allocate(choices, unname(panel))
}

pool <- c("FITC", "PE", "APC", "PerCP", "PE-Cy7", "APC-Cy7")
spectra_all <- readRDS("data/fluorochromes_all.rds")[pool]
k <- 3
mer <- 10
cases <- list(
  none = NULL,
  empty = list(),
  mandatory = list(must_include = "FITC"),
  backbone = list(must_include = c("FITC", "APC")),
  all_mandatory = list(must_include = pool[1:3]),
  prohibited = list(cannot_combine = list(c("FITC", "PE"), c("APC", "APC-Cy7"))),
  reversed_pair = list(cannot_combine = list(c("PE", "FITC"), c("APC-Cy7", "APC"))),
  one_choice = list(groups = list(list(fluorochromes = pool[1:3], n = 1))),
  repeated_choices = list(groups = list(list(fluorochromes = pool[1:4], n = 2))),
  all_positions = list(groups = list(list(fluorochromes = pool[1:4], n = 3))),
  extra_members_allowed = list(groups = list(list(fluorochromes = pool, n = 1))),
  overlap = list(groups = list(list(fluorochromes = pool[1:3], n = 2), list(fluorochromes = pool[c(1, 3)], n = 1))),
  mandatory_overlap = list(must_include = "FITC", groups = list(
    list(fluorochromes = c("FITC", "PE"), n = 1), list(fluorochromes = c("FITC", "APC"), n = 1))),
  mandatory_not_reused = list(must_include = "FITC", groups = list(list(fluorochromes = "FITC", n = 1))),
  overlap_shortage = list(groups = list(list(fluorochromes = pool[1:2], n = 2), list(fluorochromes = pool[1:2], n = 1))),
  singleton_overlap = list(groups = list(list(fluorochromes = "FITC", n = 1), list(fluorochromes = "FITC", n = 1))),
  unavailable_alternative = list(groups = list(list(fluorochromes = c("FITC", "absent"), n = 1))),
  zero_choices = list(groups = list(list(fluorochromes = pool, n = 0))),
  mixed = list(must_include = "FITC", cannot_combine = list(c("PE", "APC")), groups = list(
    list(fluorochromes = pool[1:3], n = 1), list(fluorochromes = pool[3:6], n = 1))),
  conflict = list(must_include = c("FITC", "PE"), cannot_combine = list(c("FITC", "PE"))),
  no_edges = list(cannot_combine = combn(pool, 2, simplify = FALSE)),
  missing_required = list(must_include = "absent"),
  too_many_required = list(must_include = pool),
  too_many_choices = list(groups = list(list(fluorochromes = pool, n = 4))),
  combined_overflow = list(must_include = "FITC", groups = list(list(fluorochromes = pool, n = 3))),
  unavailable_group = list(groups = list(list(fluorochromes = "absent", n = 1)))
)

for(mode in c("linear", "conventional", "spectral")) {
  instrument <- read.delim(if(mode == "spectral") "instruments/cytek_aurora_5laser_64channel.tsv" else "instruments/bd_facscantoII_8color.tsv")
  lasers <- unique(instrument$laser)
  if(mode == "spectral") stopifnot(nrow(instrument) == 64)
  spectra <- excitation_ratio_exclusion(lasers, spectra_all, mer)
  args <- list(k = k, fluorochromes_spectra = spectra)
  if(mode == "linear") args$lasers <- lasers else args$instrument <- instrument
  if(mode == "conventional") args$mer <- mer
  optimizer <- get(paste0("optimal_combi_", mode))
  exhaustive <- enumerate_panels(mode, spectra, lasers, instrument, k, mer)

  for(case in names(cases)) {
    constraints <- cases[[case]]
    allowed <- vapply(exhaustive$panels, allowed_panel, logical(1), constraints = constraints)
    result <- tryCatch(do.call(optimizer, c(args, list(constraints = constraints))), error = identity)
    if(!any(allowed)) {
      if(!inherits(result, "error")) stop(mode, "/", case, ": infeasible request returned a panel")
      next
    }
    if(inherits(result, "error")) stop(mode, "/", case, ": ", conditionMessage(result))
    stopifnot(length(result) == k, !anyDuplicated(result))
    if(mode == "conventional") stopifnot(!anyDuplicated(names(result))) else stopifnot(is.null(names(result)))
    matches <- vapply(exhaustive$panels, function(panel) {
      if(mode == "conventional") {
        setequal(paste(names(panel), panel), paste(names(result), result))
      } else setequal(panel, result)
    }, logical(1))
    stopifnot(sum(matches) == 1, allowed[matches])
    scores <- exhaustive$scores[allowed,,drop = FALSE]
    best <- scores[order(scores[,1], scores[,2])[1],]
    stopifnot(isTRUE(all.equal(unname(exhaustive$scores[matches,]), unname(best), tolerance = 1e-12)))
  }

  # Exclusion is upstream subsetting, including excluded alternatives in a group.
  subset_args <- args
  subset_args$fluorochromes_spectra <- spectra[setdiff(names(spectra), "PE")]
  excluded <- do.call(optimizer, c(subset_args, list(constraints = list(groups = list(list(fluorochromes = c("PE", "APC"), n = 1))))))
  stopifnot(!"PE" %in% excluded, "APC" %in% excluded)
  shortage_args <- args
  shortage_args$fluorochromes_spectra <- spectra[1:2]
  stopifnot(inherits(tryCatch(do.call(optimizer, shortage_args), error = identity), "error"))

  # Legacy membership semantics must be rejected, including mixed old/new input.
  for(legacy_constraints in list(NULL, list(must_include = "PerCP"))) {
    legacy <- tryCatch(do.call(optimizer, c(args, list(fridge = list("FITC", c("PE", "APC")), constraints = legacy_constraints))), error = identity)
    stopifnot(inherits(legacy, "error"), grepl("Legacy fridge", conditionMessage(legacy)))
  }
  for(rule in c("exactly", "at_least")) {
    legacy <- tryCatch(do.call(optimizer, c(args, list(constraints = list(groups = list(list(fluorochromes = pool, n = 1, rule = rule)))))), error = identity)
    stopifnot(inherits(legacy, "error"), grepl("Legacy group rules", conditionMessage(legacy)))
  }
  for(bad in list(
    list(must_include = c("FITC", "FITC")),
    list(cannot_combine = list("FITC")),
    list(cannot_combine = list(c("FITC", "FITC"))),
    list(groups = list(list(fluorochromes = c("FITC", "FITC"), n = 1))),
    list(groups = list(list(fluorochromes = pool, n = 1.5))),
    list(groups = list(list(fluorochromes = pool, n = -1))),
    list(groups = list(list(fluorochromes = pool, n = 7)))
  )) stopifnot(inherits(tryCatch(do.call(optimizer, c(args, list(constraints = bad))), error = identity), "error"))
  cat(mode, ":", length(cases), "exhaustive cases, pool exclusion, compatibility and invalid-input checks passed\n")

  # Exercise the reference script itself with a small pool and all four requirements.
  expressions <- parse("spectracular_v1.5.R")
  env <- new.env()
  for(expression in expressions) {
    if(is.call(expression) && identical(expression[[1]], as.name("<-"))) {
      variable <- as.character(expression[[2]])
      if(variable == "mode") expression[[3]] <- mode
      if(variable == "k") expression[[3]] <- k
      if(variable == "mer") expression[[3]] <- mer
      if(variable == "fluorochromes_subset") expression[[3]] <- pool
      if(variable == "constraints") expression[[3]] <- cases$mixed
    }
    eval(expression, env)
  }
  stopifnot(inherits(env$p_spectra, "patchwork"), allowed_panel(env$best_combi, cases$mixed))
  grob <- patchwork::patchworkGrob(env$p_spectra)
  stopifnot(inherits(grob, "gtable"))
  if(mode == "conventional") {
    plots <- attr(env$p_spectra, "components")$detectors
    for(i in seq_len(nrow(env$instrument))) {
      bins <- plots[[i]]$layers[[3]]$data
      assigned <- bins[bins$Assigned,]
      stopifnot(all(names(env$best_combi)[match(assigned$Sample, env$best_combi)] == assigned$Detector))
    }
  }
  cat(mode, ": reference script optimization and plot build passed\n")
}

# Check every possible edge in an expanded conventional graph. All and only
# compatible cross-slot copies must retain the exact original assignment weights.
nodes <- expand.grid(fluorochrome = LETTERS[1:4], detector = c("d1", "d2", "d3"), stringsAsFactors = FALSE)
nodes$node <- paste0("original:", seq_len(nrow(nodes)))
indices <- t(combn(seq_len(nrow(nodes)), 2))
indices <- indices[nodes$fluorochrome[indices[,1]] != nodes$fluorochrome[indices[,2]] &
  nodes$detector[indices[,1]] != nodes$detector[indices[,2]],]
edges <- data.frame(from = nodes$node[indices[,1]], to = nodes$node[indices[,2]], value = seq_len(nrow(indices)) / 7)
constraints <- list(must_include = "A", cannot_combine = list(c("B", "C")), groups = list(
  list(fluorochromes = c("A", "B", "C"), n = 1), list(fluorochromes = c("B", "C", "D"), n = 1)))
expanded <- constrain_panel_graph(4, edges, nodes, constraints)
stopifnot(setequal(expanded$nodes$slot, 1:4))
base_weights <- matrix(0, nrow(nodes), nrow(nodes), dimnames = list(nodes$node, nodes$node))
base_weights[as.matrix(edges[,1:2])] <- edges$value
base_weights <- base_weights + t(base_weights)
actual_weights <- matrix(0, nrow(expanded$nodes), nrow(expanded$nodes), dimnames = list(expanded$nodes$node, expanded$nodes$node))
actual_weights[as.matrix(expanded$edges[,1:2])] <- expanded$edges$value
actual_weights <- actual_weights + t(actual_weights)
for(pair in combn(seq_len(nrow(expanded$nodes)), 2, simplify = FALSE)) {
  a <- expanded$nodes[pair[1],]
  b <- expanded$nodes[pair[2],]
  compatible <- a$slot != b$slot && a$fluorochrome != b$fluorochrome && a$detector != b$detector &&
    !setequal(c(a$fluorochrome, b$fluorochrome), c("B", "C"))
  expected <- if(compatible) base_weights[a$original, b$original] else 0
  stopifnot(identical(actual_weights[a$node, b$node], expected))
}
for(constraints in list(NULL, list(cannot_combine = list(c("B", "C"))), list(groups = list(list(fluorochromes = "A", n = 0))))) {
  simple <- constrain_panel_graph(3, edges, nodes, constraints)
  stopifnot(identical(simple$nodes$node, nodes$node), is.null(simple$nodes$slot))
}
cat("Expanded graph edge invariants, unchanged weights and unexpanded fast path passed\n")

# Construct identical known edge weights in each mode. ABC and ABD first become
# feasible at 0.25, but ABC has the smaller sum (0.625 versus 0.75).
target <- matrix(c(0, .25, .25, .25, .25, 0, .125, .25,
                   .25, .125, 0, .375, .25, .25, .375, 0), 4, 4)
dimnames(target) <- list(LETTERS[1:4], LETTERS[1:4])
members <- list()
for(i in 1:3) for(j in (i + 1):4) {
  for(repeat_index in seq_len(8 * target[i,j])) members[[length(members) + 1L]] <- c(i, j)
}
for(i in 1:4) {
  for(repeat_index in seq_len(8 - 8 * sum(target[i,]))) members[[length(members) + 1L]] <- i
}
for(mode in c("linear", "conventional", "spectral")) {
  instrument <- data.frame(laser = 405, channel = paste0("D", seq_along(members)),
    lower = 450 + 3 * seq_along(members), upper = 451 + 3 * seq_along(members))
  if(mode == "conventional") {
    instrument <- data.frame(laser = 405, channel = paste0("D", 1:4), lower = 450 + 20 * (1:4), upper = 460 + 20 * (1:4))
  }
  synthetic <- setNames(vector("list", 4), LETTERS[1:4])
  for(i in 1:4) {
    x <- data.frame(Wavelength = 300:900, Excitation = 100, Emission = 0)
    for(d in seq_len(nrow(instrument))) {
      value <- if(mode == "conventional") if(i == d) 100 else 100 * target[i,d] else 100 * (i %in% members[[d]])
      x$Emission[x$Wavelength >= instrument$lower[d] & x$Wavelength <= instrument$upper[d]] <- value
    }
    synthetic[[i]] <- x
  }
  args <- list(k = 3, fluorochromes_spectra = synthetic)
  if(mode == "linear") args$lasers <- 405 else args$instrument <- instrument
  if(mode == "conventional") args$mer <- 50
  optimizer <- get(paste0("optimal_combi_", mode))
  oracle <- enumerate_panels(mode, synthetic, 405, instrument, 3, 50)
  expected <- c(.625, .75, .875, .75)
  stopifnot(isTRUE(all.equal(oracle$scores[,2], expected, tolerance = 1e-12)))
  stopifnot(setequal(do.call(optimizer, args), c("A", "B", "C")))
  stopifnot(setequal(do.call(optimizer, c(args, list(constraints = list(must_include = "D")))), c("A", "B", "D")))
  stopifnot(setequal(do.call(optimizer, c(args, list(constraints = list(cannot_combine = list(c("B", "C")))))), c("A", "B", "D")))
  # Repeated and overlapping group choices preserve the sum tie-break too.
  constrained <- list(groups = list(list(fluorochromes = c("A", "B", "D"), n = 2), list(fluorochromes = c("A", "C"), n = 1)))
  stopifnot(setequal(do.call(optimizer, c(args, list(constraints = constrained))), c("A", "B", "C")))
  # A one-choice group does not prohibit selecting its other members elsewhere.
  stopifnot(setequal(do.call(optimizer, c(args, list(constraints = list(groups = list(list(fluorochromes = c("A", "B", "C"), n = 1)))))), c("A", "B", "C")))
  single_args <- args
  single_args$k <- 1
  single <- do.call(optimizer, c(single_args, list(constraints = list(must_include = "D"))))
  stopifnot(identical(unname(single), "D"))
  cat(mode, ": tied thresholds, minimum-sum tie-break, group allocation and singleton panel passed\n")
}
# Detector eligibility can make distinct fluorochromes compete for one detector.
# A/B can use only D1, C only D2; D has no eligible assignment at this MER.
instrument <- data.frame(laser = 405, channel = c("D1", "D2"), lower = c(480, 580), upper = c(490, 590))
synthetic <- setNames(lapply(1:4, function(i) {
  x <- data.frame(Wavelength = 300:900, Excitation = 100, Emission = 0)
  if(i < 4) {
    d <- if(i < 3) 1 else 2
    x$Emission[x$Wavelength >= instrument$lower[d] & x$Wavelength <= instrument$upper[d]] <- 100
  }
  x
}), LETTERS[1:4])
result <- optimal_combi_conventional(2, synthetic, instrument, 50, constraints = list(must_include = "A"))
stopifnot(setequal(result, c("A", "C")), !anyDuplicated(names(result)))
for(constraints in list(list(must_include = c("A", "B")), list(groups = list(list(fluorochromes = c("A", "B"), n = 2))))) {
  failure <- tryCatch(optimal_combi_conventional(2, synthetic, instrument, 50, constraints = constraints), error = identity)
  stopifnot(inherits(failure, "error"), grepl("No feasible panel", conditionMessage(failure)))
}
failure <- tryCatch(optimal_combi_conventional(2, synthetic, instrument, 50, constraints = list(must_include = "D")), error = identity)
stopifnot(inherits(failure, "error"), grepl("no eligible assignment: D", conditionMessage(failure)))
cat("Conventional detector conflicts and unavailable required assignments passed\n")
invisible(grDevices::dev.off())
