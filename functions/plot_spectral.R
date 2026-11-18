plot_spectral <- function(best_combi, spectra, instrument) {

  # ------------------------------------------------------------
  # Prepare instrument, spectra and colours
  # ------------------------------------------------------------

  instrument <- validate_instrument(instrument)
  instrument$detector <- instrument_detector_ids(instrument)
  lasers <- unique(instrument$laser)

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
    scale_color_manual(values = fluoro_colors) +
    scale_fill_manual(values = fluoro_colors) +
    scale_x_continuous(limits = c(250, 900), expand = c(0, 0)) +
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
  # B. Predicted spectral signatures
  # ------------------------------------------------------------

  detector_signal <- calculate_detector_response(
    spectra_reference,
    instrument
  )$signal

  detector_signal_norm <- detector_signal / apply(detector_signal, 1, max)

  signature_df <- as.data.frame(as.table(detector_signal_norm))
  colnames(signature_df) <- c("Sample", "detector", "signal")

  channel_meta <- instrument[,c("detector", "channel", "laser")]
  channel_meta$laser <- factor(channel_meta$laser, levels = lasers)

  signature_df <- merge(signature_df, channel_meta, by = "detector")
  signature_df$Sample <- factor(signature_df$Sample, levels = best_combi)
  signature_df$channel <- factor(
    signature_df$channel,
    levels = unique(instrument$channel)
  )

  B <- ggplot(
    signature_df,
    aes(
      x = channel,
      y = signal,
      color = Sample,
      group = Sample
    )
  ) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.5) +
    facet_wrap(~laser, scales = "free_x") +
    scale_color_manual(values = fluoro_colors) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(
      title = "Predicted spectral signatures",
      x = "Detector channel",
      y = "Normalized detector signal"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title = element_text(face = "bold")
    )


  # ------------------------------------------------------------
  # C. Cosine similarity matrix
  # ------------------------------------------------------------

  cosine <- matrix(
    NA_real_,
    nrow = nrow(detector_signal),
    ncol = nrow(detector_signal),
    dimnames = list(rownames(detector_signal), rownames(detector_signal))
  )

  for(i in seq_len(nrow(detector_signal))) {
    for(ii in seq_len(nrow(detector_signal))) {
      x <- detector_signal[i,]
      y <- detector_signal[ii,]
      cosine[i,ii] <- sum(x * y) / sqrt(sum(x^2) * sum(y^2))
    }
  }

  pairwise_cosine <- cosine[lower.tri(cosine)]
  max_cosine <- max(pairwise_cosine, na.rm = TRUE)
  mean_cosine <- mean(pairwise_cosine, na.rm = TRUE)

  df <- as.data.frame(as.table(cosine))
  colnames(df) <- c("y", "x", "similarity")
  df$x <- factor(df$x, levels = best_combi)
  df$y <- factor(df$y, levels = rev(best_combi))

  C <- ggplot(df, aes(x = x, y = y, fill = similarity)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(
      aes(label = sprintf("%.2f", similarity), color = similarity > 0.55),
      size = 3.5
    ) +
    scale_fill_gradient(
      low = "white",
      high = "darkred",
      limits = c(0, 1),
      name = "Similarity"
    ) +
    scale_color_manual(values = c("black", "white"), guide = "none") +
    coord_fixed() +
    labs(
      title = "Cosine similarity",
      subtitle = sprintf(
        "Maximum similarity = %.2f    Mean similarity = %.2f",
        max_cosine,
        mean_cosine
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

  figure <- A / B / C +
    patchwork::plot_layout(heights = c(1, 1, 2))

  # Named components let clients display panels separately; the composite remains the default.
  attr(figure, "components") <- list(spectra = A, signatures = B, matrix = C)

  return(figure)
}
