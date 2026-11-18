# Spectracular scientific context

## Purpose and scope

Spectracular is a fluorescence-cytometry panel-design implementation built around a common optimization idea: represent candidate choices as a weighted graph and identify a set of `k` fluorochromes whose worst pairwise interference or similarity is as small as possible.

The important scientific point is that the **optimization framework is shared, but the pairwise quantity being optimized is technology-specific**. Linear spectra, conventional flow cytometry, and spectral flow cytometry therefore use different representations and edge weights.

Prefer to describe Spectracular as the implementation/software. The underlying graph procedure can be referred to more generically as the optimization algorithm.

## Core optimization framework

For each mode, Spectracular constructs a weighted graph in which lower edge weights are better.

The optimizer adds edges in increasing order of **unique edge-weight thresholds**. All edges tied at the same threshold are added together. The first threshold at which a valid `k`-clique exists is therefore the minimum achievable maximum pairwise edge weight.

If several valid `k`-cliques appear at the same minimum-maximum threshold, Spectracular selects the one with the smallest sum of pairwise weights within the clique.

This is a min-max optimization with a secondary sum-based tie-break.

The shared `functions/clique_search.R` implementation searches for exact size-`k` cliques in standard R. It repeatedly removes vertices with fewer than `k - 1` remaining neighbors, then orders each remaining candidate set by a greedy proper coloring. The color number bounds the clique size within each prefix of that ordering. Partial choices extend through earlier common neighbors, and branches are discarded only when they cannot reach `k`. The optimizers retain at most two cliques at each threshold: one alone is a unique solution, while two establish a tie. At the first feasible threshold, ties are resolved by exact minimum-sum branch-and-bound over the complete thresholded graph. Size-`k` subsets of larger cliques remain eligible. Scientific weights, the sequential unique-threshold search and both optimization criteria are unchanged. When both objectives tie, the chosen equally optimal panel can differ between search implementations.

### Callable minimum-sum tie-break

`functions/minimum_sum_clique.R` supplies `minimum_sum_clique(edges, nodes, k)` for exact secondary optimization. Pass the weighted graph at the already established optimal minimax threshold: `edges` has `from`, `to` and `value` columns, and `nodes` contains the graph's node IDs. The result is one minimum-sum size-`k` clique, returned as node IDs in their input order. An infeasible graph raises an error.

```r
source("functions/clique_search.R")
source("functions/minimum_sum_clique.R")
optimal_edges <- panel_graph$edges[panel_graph$edges$value <= threshold, ]
best_nodes <- minimum_sum_clique(optimal_edges, panel_graph$nodes$node, k)
```

Keep slot and fluorochrome-detector node identities when calling the helper; map the returned IDs through the existing node metadata afterward. All graph restrictions and cross-component weights remain part of the optimization.

The helper uses compact choices from complement components when their maximum clique sizes sum to `k`. A cost lower bound prunes combinations that cannot improve the incumbent. Other graphs use a direct weighted-clique branch-and-bound search, including size-`k` subsets of larger cliques. Local-state and pair-table budgets select between these two exact search routes; they do not cap the search or approximate the objective. Floating-point margins make pruning conservative without rounding the weights or treating distinct sums as ties. No additional package is required.

All three optimizers call this helper when a tie is detected; the reference workflow and Shiny source the same implementation. `find_k_cliques()` accepts an optional `limit`, returning a deterministic prefix of its traversal; its default remains exhaustive enumeration for callers that need all cliques. The optimizers use `limit = 2`, so they do not enumerate every tie before secondary optimization. No lexicographic refinement or approximate secondary objective is introduced.

### Selection constraints

The standalone workflow and Shiny accept a `constraints` list:

- `must_include`: distinct required fluorochrome names, each represented by a singleton selection slot.
- `cannot_combine`: a list of prohibited fluorochrome pairs; all edges between their occurrences are removed.
- `groups`: entries containing `fluorochromes` and a non-negative integer `n`; each contributes `n` selection slots populated by its eligible members.

Candidate exclusion is performed by subsetting the spectra before optimization. Construct `k` slots: mandatory singletons, grouped choices, then unrestricted slots for any remaining positions. A group with `n = 0` contributes no slots and does not exclude its members. Reject requests needing more than `k` slots, unavailable mandatory candidates, insufficient group choices, or no admissible solution.

For example, two choices from `{A, B, C}` and one from `{A, C}` produce three slots: `{A, B, C}`, `{A, B, C}`, `{A, C}`. Each selected fluorochrome can fill only one slot. Mandatory singletons consume their own positions: if `A` is mandatory and a group `{A, B}` requires one choice, that group must supply `B`.

There are no edges within a slot, between copies of the same fluorochrome, or between mutually excluded fluorochromes. All other eligible edges retain their original scientific weights. Thus a `k`-clique fills all `k` slots with distinct fluorochromes, including across overlapping groups. Conventional vertices additionally retain the fluorochrome–detector assignment, and edges between assignments sharing a detector are absent. Optimized detector assignments remain in the returned names.

These are requirements for distinct choices allocated to positions, **not exact counts of membership in the completed panel**. Other positions may select additional members of any group. The rules concern fluorochromes, not marker–antibody assignments. When there are no mandatory or positive group requirements, retain the original graph representation, with prohibited-pair edges removed as needed.

The shared implementation is `functions/panel_constraints.R`. The three optimizers keep their technology-specific calculations, unique weight thresholds (including tied edges together), minimax objective and secondary sum-of-weights tie-break. Constraints are enforced by graph topology rather than membership-count filtering of cliques.

### Compatibility with earlier constraints

Earlier `fridge` inputs required exactly one member per group. Earlier `groups` carried `rule = "exactly"` or `"at_least"` and counted completed-panel membership, allowing one fluorochrome to satisfy several groups. These meanings differ from distinct choice allocation. Nonempty `fridge` and any legacy `group$rule` are rejected with instructions to review and recreate requirements. The optimizer signatures retain the `fridge` argument for this explicit compatibility error; `NULL` and empty lists impose no requirement.

New Shiny configurations store `constraint_semantics = "distinct_choices_v1"`. Files with legacy fridge entries or unversioned groups are rejected before session state changes. Files without groups remain loadable because mandatory-candidate and prohibited-pair meanings are unchanged. No automatic migration claims equivalence with the old group rules.

## Linear mode

### Scientific representation

Linear mode treats each fluorochrome primarily through its excitation and emission spectrum rather than through a detector layout.

The production optimizer uses the instrument laser wavelengths to scale each fluorochrome's emission by its strongest available excitation:

```r
lasers_er <- approx(x$Wavelength, x$Excitation, xout = lasers, rule = 2)$y
x$Emission <- x$Emission * max(lasers_er, na.rm = TRUE) / 100
```

The pairwise edge weight is the spectral similarity coefficient (SSC):

```text
SSC(i,j) = sum(min(E_i, E_j)) / min(sum(E_i), sum(E_j))
```

For non-negative spectra, this is equivalent to the previous containment-style formulation.

Lower SSC means less pairwise spectral overlap.

### Figure 1 distinction

The conceptual Figure 1 may show raw/unadjusted spectra to explain the generic graph formulation without instrument-specific effects. That pedagogical choice should not be confused with the production linear optimizer, which can incorporate the available laser wavelengths as above.

### Terminology

Use terms such as:

- spectral overlap
- SSC matrix
- pairwise overlap

Do not call the linear SSC matrix a spillover matrix.

## Conventional flow-cytometry mode

### Scientific representation

Conventional flow cytometry requires explicit fluorochrome-to-detector assignment. A candidate node is therefore not merely a fluorochrome; it is an eligible **fluorochrome-detector assignment**.

The instrument is represented by a detector table containing at least:

```text
laser    channel    lower    upper
```

For fluorochrome `i` and detector `d`, the predicted detector signal is the excitation-scaled emission integrated across the detector wavelength window:

```text
D(i,d) = integral_window E_i(lambda) * X_i(L_d)/100 dlambda
```

The same detector-response calculation should be shared between conventional and spectral modes.

### Current detector eligibility rule

An assignment is currently eligible when both conditions hold:

```r
er >= mer && emission_peak >= mer
```

where:

- `er` is excitation at the detector's laser.
- `emission_peak` is the maximum intrinsic emission value within the detector window.

This rule is intentionally retained for now.

### Known caveat in detector eligibility

The current peak-within-window criterion can make a detector eligible even when it captures little integrated signal, for example when a narrow detector window only clips the fluorochrome spectrum. This can lead to assignments with very weak primary signal and consequently very large spillover ratios.

A likely future alternative is to require predicted detector response relative to that fluorochrome's strongest detector response, for example:

```r
response_ratio <- detector_signal[i,d] / max(detector_signal[i,]) * 100
eligible <- er >= mer && response_ratio >= mer
```

This is an **open scientific decision**. Do not implement it unless explicitly requested.

### Pairwise spillover edge weight

For fluorochrome `i` assigned to detector `d_i` and fluorochrome `j` assigned to detector `d_j`:

```text
S(i -> j) = D(i,d_j) / D(i,d_i)
S(j -> i) = D(j,d_i) / D(j,d_j)
```

Spectracular uses the symmetric pairwise edge weight:

```text
w(i,j) = max(S(i -> j), S(j -> i))
```

The directional ratio is a conventional spillover concept; taking the maximum of the two directions to create a symmetric graph edge is Spectracular-specific.

### Assignment constraints

- A fluorochrome may only be selected once.
- For the supported formulation, a detector may only be assigned once.
- Therefore `k` cannot exceed the number of detectors.

Allowing shared detector/channel assignments would require biological knowledge that Spectracular does not currently model.

### Returned solution

The conventional optimizer returns fluorochromes as values and optimized detector assignments as names, e.g.:

```r
488_FL1   488_FL2   633_FL4
"FITC"    "PE"      "APC"
```

This mapping is part of the optimized result. Plotting and downstream code must preserve it. Do not replace it by recalculating each fluorochrome's strongest detector after optimization.

### Spillover values above 1

Spillover is a ratio, not a percentage bounded by 100%.

- `1` = 100%
- `2.5` = 250%
- `250` = 25,000%

Values above 1 are possible, especially for poor/subdominant detector assignments. Do not cap the underlying values at 1. Log transformation is acceptable for visualization only.

### Interpretation limits

The current model predicts detector responses from published/reference excitation and emission spectra plus instrument wavelength windows. It does not automatically capture all real-world properties of measured single-stain controls, including detector efficiencies, gain settings, spectral response, tandem behavior, fluorochrome brightness in a given conjugate, or biological marker-expression constraints.

Accordingly, predicted spillover should be described as a model-based approximation, not as a complete prediction of compensation/spreading performance.

## Spectral flow-cytometry mode

### Scientific representation

Spectral mode represents each fluorochrome by its predicted response across **all detector channels**.

For each channel, signal is calculated from:

- excitation at that channel's laser,
- fluorochrome emission spectrum,
- detector wavelength range.

For a five-laser Cytek Aurora configuration, the detector responses are concatenated into one full vector spanning all 64 channels. The channels from different lasers are not averaged into separate scores; they form one combined detector-response vector.

### Pairwise edge weight

Pairwise similarity is cosine similarity of the complete detector-response vectors:

```text
cos(i,j) = (x . y) / (||x|| ||y||)
```

Lower cosine similarity means the predicted signatures are more distinct.

The optimization uses the unnormalized detector-response vectors. Normalizing each fluorochrome by a scalar for display does not change cosine similarity, but display normalization should remain conceptually separate from the optimization calculation.

### No per-channel MER

Spectral mode does not apply a per-channel detector-eligibility threshold analogous to conventional mode. The entire multi-laser signature contributes to the cosine calculation.

A global excitation-based exclusion may still be applied upstream when appropriate.

### Interpretation limits

Cosine similarity is a pairwise signature-distinguishability measure. It should not be overinterpreted as a full prediction of spectral unmixing quality or spreading error for an entire panel. Whole-panel collinearity and experimentally measured single-stain reference controls contain information not represented by a simple pairwise cosine score.

## Instrument model

### Shared table format

Conventional and spectral modes should use the same simple detector-level instrument table:

```text
laser    channel    lower    upper
```

This format is intentionally transparent and reusable. It should be easy to create, inspect, edit, share, and version outside Spectracular.

Supplied instrument configurations should be examples/presets, not hard-coded assumptions.

#### Panel constraints

All modes use the shared selection-slot formulation described above. Instrument tables and detector eligibility remain independent of the constraint interface.

## Linear mode

Linear mode only needs a vector of laser wavelengths. At the UI/caller level it is fine to derive that vector from an instrument table:

```r
lasers <- unique(instrument$laser)
```

The linear core functions should nevertheless accept `lasers`, because that accurately reflects their scientific input.

## Software/data architecture decisions

The public release contains the standalone workflow and shared backend. The hosted Shiny interface and its deployment files are maintained separately.

- Root-level `functions/` is the single source of truth for algorithms and plotting.
- Standalone analyses and Shiny must call the same functions.
- Shared detector-response calculations belong in `instrument_utils.R`.
- Instrument configurations should be reusable tabular files, preferably TSV.
- A fluorochrome availability/inventory list may also be TSV.
- Full Spectracular configuration/session state should remain RDS.
- Shiny should not maintain separate algorithm implementations or duplicate function files.

## Terminology and manuscript framing

Prefer technology-specific names for the pairwise matrices:

- Linear: `SSC matrix`, `spectral-overlap matrix`.
- Conventional: `pairwise maximum predicted spillover`, `spillover-derived edge-weight matrix`.
- Spectral: `cosine-similarity matrix`, `predicted spectral-signature similarity`.

Avoid using `spillover matrix` as a generic name for every mode.

Avoid implying that Spectracular captures every consideration used by expert panel designers. Expert choices may incorporate antigen density, fluorochrome brightness, available conjugates, tandem stability, biological co-expression, expected spreading error, and other information outside the current optimization model.

## Current unresolved scientific decisions

The following are intentionally unresolved and should not be changed incidentally during software refactoring:

1. **Conventional detector eligibility.** The current peak-within-window MER criterion may eventually be replaced or supplemented by a relative predicted detector-response criterion.
2. **Measured reference controls.** Future versions could allow measured single-stain detector responses/signatures instead of or alongside predictions from reference spectra.
3. **Whole-panel metrics.** Current optimization is based on pairwise edge weights. Higher-order collinearity/unmixing effects are not explicitly optimized.
4. **Biological panel constraints.** Marker expression, antigen density, fluorochrome brightness, tandem behavior, and shared-channel strategies are not currently modeled.
5. **Custom spectra/signatures.** User-uploaded fluorochrome spectra or measured signatures are plausible future features but are outside the current core update.

When changing the software, preserve these as explicit limitations rather than silently introducing assumptions intended to solve them.
