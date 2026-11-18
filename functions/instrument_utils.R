# ------------------------------------------------------------
# Instrument configuration helpers
# ------------------------------------------------------------

validate_instrument <- function(instrument) {
  required <- c("laser", "channel", "lower", "upper")
  missing <- setdiff(required, names(instrument))

  if(nrow(instrument) == 0) {
    stop("Instrument table contains no detector rows.")
  }

  if(length(missing) > 0) {
    stop("Instrument table is missing required columns: ", paste(missing, collapse = ", "))
  }

  instrument$laser <- as.numeric(instrument$laser)
  instrument$channel <- as.character(instrument$channel)
  instrument$lower <- as.numeric(instrument$lower)
  instrument$upper <- as.numeric(instrument$upper)

  if(anyNA(instrument[,required])) {
    stop("Instrument table contains missing or non-numeric laser/filter values.")
  }

  if(any(instrument$channel == "")) {
    stop("Instrument channel names cannot be empty.")
  }

  if(any(instrument$lower >= instrument$upper)) {
    stop("Each detector must have lower < upper.")
  }

  detector_id <- paste0(instrument$laser, "_", instrument$channel)
  if(anyDuplicated(detector_id)) {
    stop("Each laser/channel combination must be unique.")
  }

  instrument
}


instrument_detector_ids <- function(instrument) {
  paste0(instrument$laser, "_", instrument$channel)
}


integrate_window <- function(wavelength, intensity, lower, upper) {
  lower <- max(lower, min(wavelength))
  upper <- min(upper, max(wavelength))

  if(lower >= upper) return(0)

  wl <- sort(unique(c(
    lower,
    wavelength[wavelength > lower & wavelength < upper],
    upper
  )))

  intensity_interp <- approx(wavelength, intensity, xout = wl)$y
  sum(diff(wl) * (head(intensity_interp, -1) + tail(intensity_interp, -1)) / 2)
}


calculate_detector_response <- function(fluorochromes_spectra, instrument) {
  instrument <- validate_instrument(instrument)
  detector_id <- instrument_detector_ids(instrument)

  detector_signal <- matrix(
    0,
    nrow = length(fluorochromes_spectra),
    ncol = nrow(instrument),
    dimnames = list(names(fluorochromes_spectra), detector_id)
  )

  excitation_ratio <- detector_signal
  emission_peak <- detector_signal

  for(i in seq_along(fluorochromes_spectra)) {
    x <- fluorochromes_spectra[[i]]
    x$Emission[is.na(x$Emission) | x$Emission < 0] <- 0
    x$Excitation[is.na(x$Excitation) | x$Excitation < 0] <- 0

    for(d in seq_len(nrow(instrument))) {
      er <- approx(
        x$Wavelength,
        x$Excitation,
        xout = instrument$laser[d],
        rule = 2
      )$y

      emission_excited <- x$Emission * er / 100

      detector_signal[i,d] <- integrate_window(
        x$Wavelength,
        emission_excited,
        instrument$lower[d],
        instrument$upper[d]
      )

      in_window <- x$Wavelength >= instrument$lower[d] &
        x$Wavelength <= instrument$upper[d]

      excitation_ratio[i,d] <- er
      emission_peak[i,d] <- if(any(in_window)) max(x$Emission[in_window], na.rm = TRUE) else 0
    }
  }

  list(
    signal = detector_signal,
    excitation_ratio = excitation_ratio,
    emission_peak = emission_peak
  )
}
