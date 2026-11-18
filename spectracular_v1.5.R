#################################################
#                                               #
#               Spectracular v1.5               #
#                   11/11-2026                  #
#                 Developed by:                 #
#                Lars Ronn Olsen                #
#            Mike Bogetofte Barnkob             #
#                                               #
#         compensation is for suckers           #
#                                               #
#################################################

## Turn status comments on/off
verbose <- TRUE

## Run from the repository root; all paths below are relative to that directory.


# ------------------------------------------------------------
# Input parameters
# ------------------------------------------------------------

# Mode: linear, conventional, or spectral
mode <- "conventional"

# Desired number of fluorochromes in panel
k <- 8

# Minimal excitation ratio
mer <- 25

# Instrument configuration
if(mode == "linear") {
  lasers <- c(355,405,488,532,628)
}

if(mode == "conventional") {
  instrument <- read.delim('./instruments/bd_facscantoII_8color.tsv')
  lasers <- unique(instrument$laser)
}

if(mode == "spectral") {
  instrument <- read.delim('./instruments/cytek_aurora_5laser_64channel.tsv')
  lasers <- unique(instrument$laser)
}

# Optional panel constraints. Exclude unavailable fluorochromes via fluorochromes_subset below.
# Each mandatory fluorochrome and each requested group choice occupies its own position.
# Choices must be distinct, including across overlapping groups and mandatory candidates.
# length(must_include) + sum(group n) must not exceed k; remaining positions are unrestricted.
# Additional group members may appear in other positions. Legacy rule/fridge inputs are rejected.
constraints <- NULL
# constraints <- list(
#   must_include = c("FITC", "APC"),
#   cannot_combine = list(c("PE", "PE-Cy7"), c("APC", "APC-Cy7")),
#   groups = list(
#     list(fluorochromes = c("PE", "PE-CF594", "PE-Cy5"), n = 1),
#     list(fluorochromes = c("BUV395", "BUV496", "BUV563"), n = 2)
#   )
# )


# ------------------------------------------------------------
# Load fluorochrome data
# ------------------------------------------------------------

fluorochromes_all <- readRDS('./data/fluorochromes_all.rds')

fluorochromes_subset <- c("Alexa Fluor 350", "Alexa Fluor 405", "Alexa Fluor 488", "Alexa Fluor 514", "Alexa Fluor 610", "Alexa Fluor 633", "Alexa Fluor 700", "Alexa Fluor 790", "APC-H7", "APC-R700", "APC", "Brilliant Blue 515", "Brilliant Blue 700", "Brilliant Violet 421", "Brilliant Violet 480", "Brilliant Violet 510", "Brilliant Violet 605", "Brilliant Violet 650", "Brilliant Violet 711", "Brilliant Violet 750", "Brilliant Violet 786", "BUV395", "BUV496", "BUV563", "BUV661", "BUV737", "BUV805", "FITC", "Horizon V450", "Horizon V500", "PE-CF594", "PE-Cy5.5", "PE-Cy5", "PE-Cy7", "PE", "PerCP-Cy5.5", "PerCP-eFluor 710", "PerCP", "Alexa Fluor 680", "APC-Cy7", "APC-Fire 810", "Brilliant Violet 570", "BUV615", "PE-Vio 615", "PE-Vio770", "PerCP-Vio 700", "Super Bright 436", "Super Bright 600", "Super Bright 645", "Super Bright 702", "Super Bright 780", "VioBlue", "VioBright 515", "BYG584")

fluorochromes_spectra <- fluorochromes_all[
  names(fluorochromes_all) %in% fluorochromes_subset
]


# ------------------------------------------------------------
# Load packages and functions
# ------------------------------------------------------------

suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))

source('./functions/instrument_utils.R')
source('./functions/panel_constraints.R')
source('./functions/clique_search.R')
source('./functions/minimum_sum_clique.R')
source('./functions/excitation_ratio_exclusion.R')
source('./functions/wavelength_to_rgb.R')
source('./functions/best_combi_linear.R')
source('./functions/best_combi_conventional.R')
source('./functions/best_combi_spectral.R')
source('./functions/plot_linear.R')
source('./functions/plot_conventional.R')
source('./functions/plot_spectral.R')


# ------------------------------------------------------------
# Filter fluorochromes by excitation
# ------------------------------------------------------------

fluorochromes_spectra <- excitation_ratio_exclusion(
  lasers,
  fluorochromes_spectra,
  mer
)


# ------------------------------------------------------------
# Calculate optimal combination and plot results
# ------------------------------------------------------------

if(mode == "linear") {
  best_combi <- optimal_combi_linear(
    k,
    fluorochromes_spectra,
    lasers,
    constraints = constraints
  )

  p_spectra <- plot_linear(
    best_combi,
    fluorochromes_spectra,
    lasers
  )
}

if(mode == "conventional") {
  best_combi <- optimal_combi_conventional(
    k,
    fluorochromes_spectra,
    instrument,
    mer,
    constraints = constraints
  )

  p_spectra <- plot_conventional(
    best_combi,
    fluorochromes_spectra,
    instrument
  )
}

if(mode == "spectral") {
  best_combi <- optimal_combi_spectral(
    k,
    fluorochromes_spectra,
    instrument,
    constraints = constraints
  )

  p_spectra <- plot_spectral(
    best_combi,
    fluorochromes_spectra,
    instrument
  )
}
