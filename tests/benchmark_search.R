# Run from the repository root: Rscript tests/benchmark_search.R
# Spectral mode, no constraints. Edit these cases and settings before running.
cases <- data.frame(n = c(120, 120, 120), k = c(20, 40, 50))
repetitions <- 1
timeout_seconds <- 600 # Whole child process; 0 disables the limit.
seed <- 1             # Repetitions use the same pool, including its ordering.
mer <- 0
instrument_file <- "instruments/cytek_aurora_5laser_64channel.tsv"

args <- commandArgs(trailingOnly = TRUE)
if(length(args) && args[1] == "--worker") {
  n <- as.integer(args[2]); k <- as.integer(args[3]); method <- args[4]
  for(f in c("instrument_utils", "panel_constraints", "clique_search", "minimum_sum_clique", "excitation_ratio_exclusion", "best_combi_spectral")) source(paste0("functions/", f, ".R"))

  # Swap only clique enumeration; weights, thresholds and tie-breaking stay shared.
  if(method == "igraph") {
    invisible(loadNamespace("igraph"))
    find_k_cliques <- function(edges, nodes, k, limit = Inf) {
      graph <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
      # This igraph API enumerates before truncating; standard R can stop early.
      cliques <- lapply(igraph::cliques(graph, min = k, max = k), function(x) nodes[match(names(x), as.character(nodes))])
      cliques[seq_len(min(length(cliques), limit))]
    }
  }
  instrument <- read.delim(instrument_file)
  spectra <- excitation_ratio_exclusion(unique(instrument$laser), readRDS("data/fluorochromes_all.rds"), mer)
  stopifnot(n >= k, k >= 2, n <= length(spectra))
  set.seed(seed)
  if(n < length(spectra)) spectra <- spectra[sample(names(spectra), n)]

  # Time the full optimizer, excluding setup and the independent score check.
  elapsed <- system.time(panel <- optimal_combi_spectral(k, spectra, instrument))["elapsed"]
  stopifnot(length(panel) == k, !anyDuplicated(panel))
  signal <- calculate_detector_response(spectra[panel], instrument)$signal
  cosine <- tcrossprod(signal / sqrt(rowSums(signal^2)))
  weights <- cosine[lower.tri(cosine)]
  saveRDS(data.frame(seconds = unname(elapsed), maximum = max(weights), total = sum(weights)), args[5])
} else {
  if(!requireNamespace("igraph", quietly = TRUE)) stop("Install the benchmark dependency first: install.packages('igraph')")
  cat(R.version.string, "| igraph", as.character(packageVersion("igraph")), "\n")
  results <- data.frame()
  for(i in seq_len(nrow(cases))) for(r in seq_len(repetitions)) {
    methods <- if(r %% 2) c("igraph", "exact_R") else c("exact_R", "igraph")
    for(method in methods) {
      cat(sprintf("n=%d, k=%d, repetition=%d, %s ... ", cases$n[i], cases$k[i], r, method)); flush.console()
      output <- tempfile(); log <- tempfile()
      command <- c("tests/benchmark_search.R", "--worker", cases$n[i], cases$k[i], method, shQuote(output))
      status <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), command, stdout = log, stderr = log, timeout = timeout_seconds))
      row <- data.frame(n = cases$n[i], k = cases$k[i], repetition = r, method = method, seconds = NA_real_, maximum = NA_real_, total = NA_real_, status = "ok", objectives_agree = NA)
      if(status == 0) {
        row[c("seconds", "maximum", "total")] <- readRDS(output)
        cat(sprintf("%.3f seconds\n", row$seconds))
      } else {
        row$status <- if(status == 124) paste0("timeout (", timeout_seconds, "s)") else "error"
        cat(row$status, "\n")
        if(status != 124) cat(readLines(log, warn = FALSE), sep = "\n")
      }
      results <- rbind(results, row)
      unlink(c(output, log))
    }
    pair <- (nrow(results) - 1):nrow(results)
    if(all(results$status[pair] == "ok")) {
      agree <- all(abs(c(diff(results$maximum[pair]), diff(results$total[pair]))) <= 1e-12)
      results$objectives_agree[pair] <- agree
      if(!agree) warning("Objective mismatch for n=", cases$n[i], ", k=", cases$k[i])
    }
  }
  # Console output only; 'results' remains available when run with source().
  print(results, row.names = FALSE, digits = 6)
}
