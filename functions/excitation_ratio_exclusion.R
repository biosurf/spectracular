excitation_ratio_exclusion <- function(lasers, fluorochromes_spectra, mer) {
  
  fluorochromes_spectra_mer_included <- list()
  
  for (i in seq_along(fluorochromes_spectra)) {
    
    x <- fluorochromes_spectra[[i]]
    
    # Excitation ratio at each instrument laser
    lasers_er <- approx(
      x = x$Wavelength,
      y = x$Excitation,
      xout = lasers,
      rule = 2
    )$y
    
    # Include if sufficiently excited by at least one laser
    if (max(lasers_er, na.rm = TRUE) >= mer) {
      
      fluorochromes_spectra_mer_included[
        names(fluorochromes_spectra)[i]
      ] <- list(x)
      
    } else {
      
      print(
        paste0(
          names(fluorochromes_spectra)[i],
          " omitted. Maximum excitation ratio is below ",
          mer,
          "%"
        )
      )
    }
  }
  
  return(fluorochromes_spectra_mer_included)
}