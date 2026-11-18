plot_linear <- function(best_combi, spectra, lasers) {
  
  # ------------------------------------------------------------
  # Recalculate emission spectra based on laser configuration
  # ------------------------------------------------------------
  
  spectra_emission <- spectra
  
  for(i in seq_along(spectra_emission)) {
    x <- spectra_emission[[i]]
    x$Emission[is.na(x$Emission) | x$Emission < 0] <- 0
    x$Excitation[is.na(x$Excitation) | x$Excitation < 0] <- 0
    
    lasers_er <- approx(x$Wavelength, x$Excitation, xout = lasers, rule = 2)$y
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
  # Prepare spectra for plotting
  # ------------------------------------------------------------
  
  plot_data <- do.call(rbind, lapply(best_combi, function(sample) {
    x <- spectra_emission[[sample]]
    color <- wavelength_to_rgb(x$Wavelength[which.max(x$Emission)])
    
    rbind(
      data.frame(
        Wavelength = x$Wavelength,
        Intensity = x$Excitation,
        Type = "Excitation",
        Sample = sample,
        Color = color
      ),
      data.frame(
        Wavelength = x$Wavelength,
        Intensity = x$Emission,
        Type = "Emission",
        Sample = sample,
        Color = color
      )
    )
  }))
  
  # Label fluorochromes at their emission maxima
  label_data <- do.call(rbind, lapply(best_combi, function(sample) {
    x <- spectra_emission[[sample]]
    data.frame(
      Wavelength = x$Wavelength[which.max(x$Emission)],
      Sample = sample
    )
  }))
  
  
  # ------------------------------------------------------------
  # Plot spectra
  # ------------------------------------------------------------
  
  A <- ggplot(
    plot_data,
    aes(
      x = Wavelength,
      y = Intensity,
      group = interaction(Sample, Type),
      color = Color,
      fill = Color
    )
  ) +
    geom_line(
      data = subset(plot_data, Type == "Excitation"),
      linetype = "dotted",
      linewidth = 0.3
    ) +
    geom_ribbon(
      data = subset(plot_data, Type == "Emission"),
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
      aes(x = Wavelength, y = 20, label = Sample),
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
    scale_color_identity() +
    scale_fill_identity() +
    scale_x_continuous(limits = c(250, 900), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    labs(
      title = "Linear spectra",
      x = "Wavelength (nm)",
      y = "Relative intensity (%)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.margin = margin(1, 1, 1, 1, "cm"),
      plot.title = element_text(face = "bold")
    )
  
  
  # ------------------------------------------------------------
  # Plot SSC matrix
  # ------------------------------------------------------------
  
  combi_ssc <- ssc[best_combi,best_combi,drop = FALSE]
  pairwise_ssc <- combi_ssc[lower.tri(combi_ssc)]
  
  max_ssc <- max(pairwise_ssc, na.rm = TRUE)
  mean_ssc <- mean(pairwise_ssc, na.rm = TRUE)
  
  df <- as.data.frame(as.table(combi_ssc))
  colnames(df) <- c("y", "x", "overlap")
  df$x <- factor(df$x, levels = best_combi)
  df$y <- factor(df$y, levels = rev(best_combi))
  
  B <- ggplot(df, aes(x = x, y = y, fill = overlap)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(
      aes(label = sprintf("%.2f", overlap), color = overlap > 0.55),
      size = 3.5
    ) +
    scale_fill_gradient(
      low = "white",
      high = "darkred",
      limits = c(0, 1),
      name = "SSC"
    ) +
    scale_color_manual(values = c("black", "white"), guide = "none") +
    coord_fixed() +
    labs(
      title = "Spectral overlap (SSC)",
      subtitle = sprintf(
        "Maximum SSC = %.2f    Mean SSC = %.2f",
        max_ssc,
        mean_ssc
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
  
  figure <- A / B + patchwork::plot_layout(heights = c(1, 2))
  
  # Named components let clients display panels separately; the composite remains the default.
  attr(figure, "components") <- list(spectra = A, matrix = B)

  return(figure)
}