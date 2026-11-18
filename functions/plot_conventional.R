plot_conventional <- function(best_combi, spectra, instrument) {

  # ------------------------------------------------------------
  # Prepare instrument, assignments, spectra and colours
  # ------------------------------------------------------------

  detectors <- validate_instrument(instrument)
  detectors$detector <- instrument_detector_ids(detectors)
  lasers <- unique(detectors$laser)

  # best_combi is named by the optimized detector assignment.
  assigned_detector <- setNames(names(best_combi), unname(best_combi))
  best_combi <- unname(best_combi)

  wl_min <- min(sapply(spectra, function(x) min(x$Wavelength)))
  wl_max <- max(sapply(spectra, function(x) max(x$Wavelength)))

  detectors$plot_lower <- pmax(detectors$lower, wl_min)
  detectors$plot_upper <- pmin(detectors$upper, wl_max)
  detectors$midpoint <- (detectors$plot_lower + detectors$plot_upper) / 2

  fluoro_colors <- sapply(best_combi, function(sample) {
    x <- spectra[[sample]]
    wavelength_to_rgb(x$Wavelength[which.max(x$Emission)])
  })
  names(fluoro_colors) <- best_combi

  spectra_reference <- spectra[best_combi]

  for(i in seq_along(spectra_reference)) {
    x <- spectra_reference[[i]]
    x$Emission[is.na(x$Emission) | x$Emission < 0] <- 0
    x$Excitation[is.na(x$Excitation) | x$Excitation < 0] <- 0
    spectra_reference[[i]] <- x
  }


  # ------------------------------------------------------------
  # A. Excitation and emission spectra
  # ------------------------------------------------------------

  linear_data <- do.call(rbind, lapply(best_combi, function(sample) {
    x <- spectra_reference[[sample]]

    rbind(
      data.frame(
        Wavelength = x$Wavelength,
        Intensity = x$Excitation,
        Type = "Excitation",
        Sample = sample
      ),
      data.frame(
        Wavelength = x$Wavelength,
        Intensity = x$Emission,
        Type = "Emission",
        Sample = sample
      )
    )
  }))

  label_data <- do.call(rbind, lapply(best_combi, function(sample) {
    x <- spectra_reference[[sample]]
    data.frame(
      Wavelength = x$Wavelength[which.max(x$Emission)],
      Sample = sample
    )
  }))

  A <- ggplot(
    linear_data,
    aes(
      x = Wavelength,
      y = Intensity,
      group = interaction(Sample, Type),
      color = Sample,
      fill = Sample
    )
  ) +
    geom_rect(
      data = detectors,
      aes(xmin = plot_lower, xmax = plot_upper, ymin = -Inf, ymax = Inf),
      fill = "grey50",
      colour = "black",
      alpha = 0.15,
      inherit.aes = FALSE
    ) +
    geom_line(
      data = subset(linear_data, Type == "Excitation"),
      linetype = "dotted",
      linewidth = 0.3
    ) +
    geom_ribbon(
      data = subset(linear_data, Type == "Emission"),
      aes(ymin = 0, ymax = Intensity),
      alpha = 0.5,
      color = NA
    ) +
    geom_vline(
      xintercept = lasers,
      linetype = "dashed",
      color = "darkgrey",
      linewidth = 0.5
    ) +
    geom_text(
      data = label_data,
      aes(x = Wavelength, y = 50, label = Sample),
      angle = 90,
      size = 3,
      color = "black",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = data.frame(Wavelength = lasers),
      aes(x = Wavelength, y = 90, label = paste(Wavelength, "nm")),
      angle = 45,
      size = 3,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = detectors,
      aes(x = midpoint, y = 90, label = channel),
      size = 3,
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = fluoro_colors) +
    scale_fill_manual(values = fluoro_colors) +
    scale_x_continuous(limits = c(300, 900), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    labs(
      title = "Excitation and emission spectra",
      x = "Wavelength (nm)",
      y = "Relative intensity (%)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(face = "bold")
    )


  # ------------------------------------------------------------
  # B. Detector-specific signals
  # ------------------------------------------------------------

  response <- calculate_detector_response(spectra_reference, detectors)
  detector_signal <- response$signal

  detector_spectra <- list()
  detector_bins <- list()
  n <- 1

  for(d in seq_len(nrow(detectors))) {
    detector <- detectors$detector[d]

    for(sample in best_combi) {
      x <- spectra_reference[[sample]]
      er <- response$excitation_ratio[sample, detector]
      emission_excited <- x$Emission * er / 100

      lower <- detectors$plot_lower[d]
      upper <- detectors$plot_upper[d]
      signal <- detector_signal[sample, detector]

      detector_spectra[[n]] <- data.frame(
        Wavelength = x$Wavelength,
        Intensity = emission_excited,
        Sample = sample,
        Detector = detector
      )

      detector_bins[[n]] <- data.frame(
        lower = lower,
        upper = upper,
        midpoint = (lower + upper) / 2,
        Signal = signal / (upper - lower),
        Sample = sample,
        Detector = detector,
        Assigned = assigned_detector[sample] == detector
      )

      n <- n + 1
    }
  }

  detector_spectra <- do.call(rbind, detector_spectra)
  detector_bins <- do.call(rbind, detector_bins)

  detector_plots <- lapply(seq_len(nrow(detectors)), function(d) {
    detector <- detectors$detector[d]
    spec <- detector_spectra[detector_spectra$Detector == detector,]
    bins <- detector_bins[detector_bins$Detector == detector,]

    ggplot(spec, aes(x = Wavelength, y = Intensity, color = Sample, fill = Sample)) +
      geom_ribbon(aes(ymin = 0, ymax = Intensity), alpha = 0.35, color = NA) +
      geom_line(alpha = 0.25, linewidth = 0.3) +
      geom_rect(
        data = bins,
        aes(xmin = lower, xmax = upper, ymin = 0, ymax = Signal, fill = Sample),
        alpha = 0.35,
        color = NA,
        inherit.aes = FALSE
      ) +
      geom_segment(
        data = bins,
        aes(
          x = lower,
          xend = upper,
          y = Signal,
          yend = Signal,
          color = Sample,
          linewidth = Assigned
        ),
        inherit.aes = FALSE
      ) +
      scale_linewidth_manual(
        values = c(`FALSE` = 0.8, `TRUE` = 1.5),
        guide = "none"
      ) +
      scale_color_manual(values = fluoro_colors) +
      scale_fill_manual(values = fluoro_colors) +
      scale_x_continuous(
        limits = c(detectors$plot_lower[d], detectors$plot_upper[d]),
        expand = c(0, 0)
      ) +
      scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
      labs(
        title = paste0(detectors$channel[d], " (", detectors$laser[d], " nm)"),
        x = "Wavelength (nm)",
        y = "Relative intensity (%)"
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.title = element_text(face = "bold")
      )
  })

  detector_grid <- patchwork::wrap_plots(detector_plots, ncol = 4)


  # ------------------------------------------------------------
  # C. Predicted spillover for optimized assignments
  # ------------------------------------------------------------

  assigned_index <- match(
    assigned_detector[best_combi],
    colnames(detector_signal)
  )

  assigned_signal <- detector_signal[
    cbind(seq_along(best_combi), assigned_index)
  ]

  directional_spillover <- sapply(assigned_index, function(d) {
    detector_signal[,d] / assigned_signal
  })

  rownames(directional_spillover) <- best_combi
  colnames(directional_spillover) <- best_combi

  spillover <- pmax(directional_spillover, t(directional_spillover))

  pairwise_spillover <- spillover[lower.tri(spillover)]
  max_spillover <- max(pairwise_spillover)
  mean_spillover <- mean(pairwise_spillover)

  df <- as.data.frame(as.table(spillover))
  colnames(df) <- c("y", "x", "spillover")
  df$x <- factor(df$x, levels = best_combi)
  df$y <- factor(df$y, levels = rev(best_combi))

  C <- ggplot(df, aes(x = x, y = y, fill = spillover)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(
      aes(label = sprintf("%.2f", spillover), color = spillover > 0.55 * max(spillover)),
      size = 3.5
    ) +
    scale_fill_gradient(
      low = "white",
      high = "darkred",
      limits = c(0, max(spillover)),
      name = "Spillover"
    ) +
    scale_color_manual(values = c("black", "white"), guide = "none") +
    coord_fixed() +
    labs(
      title = "Pairwise maximum predicted spillover",
      subtitle = sprintf(
        "Maximum spillover = %.2f    Mean spillover = %.2f",
        max_spillover,
        mean_spillover
      ),
      x = NULL,
      y = NULL
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(hjust = 1),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      plot.title = element_text(face = "bold")
    )


  # ------------------------------------------------------------
  # Combine plots
  # ------------------------------------------------------------

  n_detector_rows <- ceiling(nrow(detectors) / 4)

  figure <- A / detector_grid / C +
    patchwork::plot_layout(heights = c(1, n_detector_rows, 2))

  # Named components let clients display panels separately; the composite remains the default.
  attr(figure, "components") <- list(spectra = A, detectors = detector_grid, matrix = C)

  return(figure)
}
