# Spectracular

**[Launch the Spectracular web app](https://lronn-spectracular.share.connect.posit.cloud/)** — hosted on Posit Connect Cloud; no local installation needed.

Spectracular optimizes fluorochrome panels for conventional and spectral flow cytometry, with an additional full-spectra analysis mode. This repository contains the v1.5 standalone R workflow and shared optimization functions. Use the hosted Shiny app above for interactive analyses; its interface and deployment files are maintained separately.

## Analysis modes

| Mode | Instrument input | Pairwise quantity minimized |
| --- | --- | --- |
| **Conventional** | Laser and detector configuration | Maximum predicted spillover in either direction, with fluorochrome–detector assignments optimized jointly |
| **Spectral** | Laser and detector configuration | Cosine similarity between full detector-response vectors |
| **Full spectra** (`linear` in R) | Laser wavelengths | Emission overlap (Szymkiewicz–Simpson coefficient, SSC) after excitation scaling |

In every mode, the optimization algorithm minimizes the largest pairwise value in the selected panel, then minimizes the sum of pairwise values to break ties. Conventional results retain the optimized detector assignment for each fluorochrome.

## Using the web app

1. Choose an analysis mode.
2. Select an instrument preset or upload an instrument TSV. For Full spectra, enter the laser wavelengths directly.
3. Choose the available fluorochromes, panel size and minimum excitation ratio. Add any panel constraints.
4. Run the analysis and inspect the spectra, selected panel, pairwise matrix and excitation table.
5. Use **Save configuration** to download your settings. **Load configuration** restores them for a later analysis.

Panel constraints support mandatory fluorochromes, prohibited pairs and a specified number of choices from each group. Each mandatory fluorochrome and each required group choice occupies a separate position. Overlapping groups require distinct choices, and a fluorochrome can appear only once. Other positions may select additional group members. Exclude candidates through **Available fluorochromes**, or by subsetting the spectra in R.

See the [user guide](help.md) for detailed instructions, configuration compatibility and interpretation of the results.

## Running locally

Download or clone this repository, keeping its directory structure intact. Open R with the repository root as the working directory. Install the standalone optimizer and plotting dependencies:

```r
install.packages(c("ggplot2", "patchwork"))
```

The shared optimizers use standard R for clique feasibility checks, exact minimum-sum branch-and-bound and matrix-to-edge conversion; no external graph package or `reshape2` is required. Plotting uses `ggplot2` and `patchwork`. Scripts that source the backend directly should also source `functions/clique_search.R` and `functions/minimum_sum_clique.R`, alongside the existing instrument and constraint helpers. The minimax objective and minimum-sum tie-break are unchanged.

Open [spectracular_v1.5.R](spectracular_v1.5.R) and set the analysis mode, instrument, panel size, excitation threshold, candidate pool and optional constraints. All file paths are relative to the repository root. Run the script in R:

```r
source("spectracular_v1.5.R")
```

The selected panel is returned in `best_combi`, with detector assignments in its names for Conventional mode. Display the resulting figure with `print(p_spectra)`. The standalone workflow was checked with R 4.5.1.

## Checks

Run the backend checks from the repository root:

```sh
Rscript tests/clique_search.R
Rscript tests/minimum_sum_clique.R
Rscript tests/panel_constraints.R
```

The optional `tests/benchmark_search.R` compares the standard-R search with igraph and requires installing igraph separately. Edit its candidate counts, panel sizes and time limit before running it; the supplied large cases can take several minutes. This comparison dependency is not required by the optimizer or the checks above.

## Instrument configurations

Conventional and Spectral modes use TSV files with one row per detector and these required columns:

```text
laser    channel    lower    upper
```

Laser wavelengths and detector bounds are in nanometres. The [instrument library](instruments/README.md) documents the bundled configurations, sources and modelling assumptions. Presets are starting points: check them against the lasers and filters installed on your instrument before use.

Detector responses are estimates based on reference spectra and detector windows. The optimized pairwise quantities do not capture all biological requirements, measured spreading effects or whole-panel unmixing quality. See the [scientific context](docs/SCIENTIFIC_CONTEXT.md) for the model and its assumptions.

## Documentation

- [User guide](help.md)
- [Instrument library and sources](instruments/README.md)
- [Scientific context](docs/SCIENTIFIC_CONTEXT.md)
- [v1.5 changes](docs/V1.5_CHANGES.md)
- [Development changelog](docs/CHANGELOG.md)

## Authors and attribution

Developed by Lars Ronn Olsen and Mike Bogetofte Barnkob. When reporting analyses, reference the [Spectracular repository](https://github.com/biosurf/spectracular) and record the version or commit and instrument configuration used.

## License and commercial licensing

Original Spectracular code is provided under the **[PolyForm Noncommercial License 1.0.0](LICENSE)**. It permits use, modification and redistribution for the purposes allowed by that license. Its permissions expressly include use by educational institutions and public research organizations regardless of funding source; the full license text governs.

**Commercial use outside those permissions requires a separate written license from the relevant copyright holders.** The copyright holders retain the ability to offer separate commercial licenses for their original Spectracular code. For commercial licensing enquiries, contact the maintainers through the [GitHub repository](https://github.com/biosurf/spectracular).

Third-party packages and data retain their own licensing terms. Neither this license nor a commercial agreement for Spectracular replaces those terms or grants rights the Spectracular licensors do not hold.
