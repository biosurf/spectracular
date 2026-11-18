# AGENTS.md

## Purpose

This repository contains the public standalone release of Spectracular, an R implementation for fluorochrome panel optimization across linear, conventional, and spectral flow-cytometry settings. The hosted Shiny interface is maintained separately.

Before changing algorithmic behavior, read `docs/SCIENTIFIC_CONTEXT.md`. Do not alter scientific assumptions, optimization criteria, or terminology unless the task explicitly asks for it.

## Repository architecture

Use the root-level directories as the single source of truth:

```text
Spectracular/
├── AGENTS.md
├── spectracular_v1.5.R
├── functions/
├── data/
├── instruments/
├── help.md
├── tests/
└── docs/
    ├── SCIENTIFIC_CONTEXT.md
    ├── CHANGELOG.md
    └── V1.5_CHANGES.md
```

Key rules:

- `functions/` is the canonical implementation of all optimization, instrument, and plotting functions.
- `spectracular_v1.5.R` is the canonical standalone/reference workflow and must source functions from `functions/`.
- Run the standalone workflow and tests from the repository root, using explicit relative paths without project-root discovery or machine-specific working directories.
- `data/` contains shared data such as `fluorochromes_all.rds`.
- `instruments/` contains reusable instrument configuration tables and example presets.
- `docs/` contains scientific context, implementation plans, and development history intended to guide Codex and human contributors.
- Manuscript-specific figure and analysis scripts live outside this repository and should source the canonical functions from this project.
- Do not create an `analyses/` directory unless a future need clearly justifies it.
- Do not convert the project into an R package unless explicitly requested.

## Canonical interfaces

Keep function inputs faithful to the scientific model rather than forcing artificial uniformity.

### Linear

Linear mode needs only laser wavelengths.

```r
optimal_combi_linear(k, fluorochromes_spectra, lasers, fridge = NULL, constraints = NULL)
plot_linear(best_combi, spectra, lasers)
```

### Conventional

Conventional mode uses a detector-level instrument table.

```r
optimal_combi_conventional(k, fluorochromes_spectra, instrument, mer, fridge = NULL, constraints = NULL)
plot_conventional(best_combi, spectra, instrument)
```

The optimizer returns fluorochrome names as values and optimized detector assignments as names. Preserve this mapping.

### Spectral

Spectral mode uses the same detector-level instrument table format.

```r
optimal_combi_spectral(k, fluorochromes_spectra, instrument, fridge = NULL, constraints = NULL)
plot_spectral(best_combi, spectra, instrument)
```

`constraints` optionally contains `must_include`, `cannot_combine`, and `groups` with `fluorochromes` and a non-negative integer `n`. Each mandatory fluorochrome and each requested group choice occupies a separate selection slot; overlapping groups require distinct fluorochromes. Remaining positions are unrestricted and may select additional group members. The standalone workflow and Shiny source `functions/panel_constraints.R`. Legacy `fridge` and `group$rule` inputs are rejected with compatibility guidance rather than reinterpreted. Availability/exclusion is handled by subsetting the input spectra.

## Instrument configuration

Conventional and spectral instrument tables use at least these columns:

```text
laser    channel    lower    upper
```

- `laser`: excitation wavelength in nm.
- `channel`: detector/channel name.
- `lower`, `upper`: detector wavelength bounds in nm.
- Additional columns may be preserved, but code should not require them unless explicitly needed.

Shared detector-response logic belongs in `functions/instrument_utils.R` rather than being reimplemented inside optimizers or plotting functions.

Expected shared helpers include:

```r
validate_instrument()
instrument_detector_ids()
integrate_window()
calculate_detector_response()
```

Linear mode may derive a laser vector from an instrument table at the caller/UI level using `unique(instrument$laser)`, but its core functions should continue to accept `lasers` directly.

## Algorithmic invariants

Do not change these without explicit instruction:

- Optimization is based on technology-specific pairwise edge weights.
- Edges are added by **unique weight threshold**, so all tied edges at a threshold enter simultaneously.
- The first threshold containing a valid `k`-clique gives the minimum possible maximum pairwise weight.
- If several solutions share that maximum, break the tie by minimizing the sum of pairwise weights within the solution.
- Selection requirements allocate distinct choices to `k` slots. No edges connect nodes in the same slot, copies of the same fluorochrome, or prohibited pairs. Conventional graphs also exclude repeated detectors. No completed-panel membership-count rules or cross-group reuse are supported.
- Keep the `fridge` argument only to give a clear compatibility error for nonempty legacy inputs; never silently reinterpret old group rules or saved settings.
- Conventional mode cannot select more fluorochromes than there are detectors.
- Conventional detector assignments are part of the optimized solution. Plotting must use those assignments rather than replacing them with the strongest detector after optimization.
- All modes use the standard-R search in `functions/clique_search.R`, retaining at most two cliques per threshold to distinguish a unique solution from a tie. Ties use exact minimum-sum branch-and-bound in `functions/minimum_sum_clique.R`. Prune only branches that cannot contain a size-k clique or whose sum lower bound cannot improve the incumbent. Include size-k subsets of larger cliques; do not approximate the sum or add lexicographic refinement.
- Do not cap conventional spillover values at 1.
- Do not change the current conventional detector-eligibility/MER rule unless explicitly requested.

## Public release scope

The hosted Shiny app is a client of the shared backend, maintained in the private development project. This public release contains the standalone workflow, backend, data, instruments, backend tests and documentation.

- Keep Shiny source, deployment files and Shiny-dependent tests out of this release.
- Preserve the hosted-app link and its user guide in `help.md`.
- Prefer reusable TSV files for instrument configurations and fluorochrome inventories.

## Coding style

Keep R code lean, readable, and well commented.

- Prefer short functions that reflect the actual algorithmic complexity.
- Comment scientific or non-obvious logic; do not narrate every trivial line.
- Avoid unnecessary wrappers, abstractions, type coercions, and defensive checks when inputs are already controlled by the caller.
- Reuse shared helpers instead of duplicating calculations.
- Preserve existing naming unless a rename materially improves consistency.
- Avoid splitting simple `[[...]]` expressions across multiple lines.
- Do not introduce parallelization unless explicitly requested.
- Do not silently change formulas, thresholds, tie-breaking, or returned object structure while refactoring.

## Change discipline

For substantial tasks:

1. Read `AGENTS.md` and the relevant documents in `docs/` first.
2. Inspect the current implementation before editing.
3. Make the smallest coherent change that satisfies the task.
4. Keep `spectracular_v1.5.R` on the shared backend in `functions/`.
5. If the current code conflicts with these instructions or the scientific context, report the mismatch rather than silently choosing one interpretation.
6. At the end, report which files changed, checks performed, and any unresolved issues.
7. If the task changed the repository, append a concise entry to `docs/CHANGELOG.md`.

## Change log

`docs/CHANGELOG.md` is the persistent development record for completed implementation cycles. It complements Git history by recording intent and context rather than line-level diffs.

After every completed implementation cycle that changes the repository, append a new entry containing:

- Date.
- Short description of the requested change.
- Files changed.
- Important behavioral or architectural changes.
- Tests or checks performed.
- Unresolved issues or follow-up decisions.

A cycle is one user-requested implementation task or one coherent set of closely related changes.

Keep entries concise. Do not use the changelog as a detailed diff; Git remains the source of truth for exact code changes. Do not rewrite previous entries unless correcting a factual error.

## Validation expectations

After backend changes, check the three modes independently through `spectracular_v1.5.R` where practical:

- Linear: optimization and plotting with a laser vector.
- Conventional: optimization and plotting with a conventional instrument preset; confirm returned detector assignments are used in the plot.
- Spectral: optimization and plotting with a spectral instrument preset; confirm full detector-response vectors are used for cosine similarity.

After path or release-layout changes, run the standalone workflow and backend tests from a separate copy of the repository to check that no files outside the release are required.
