# Using Spectracular

This guide describes the [hosted Spectracular web app](https://lronn-spectracular.share.connect.posit.cloud/). For the standalone R workflow included in this repository, see the [README](README.md#running-locally).

Choose an analysis mode, configure your instrument and panel settings, then click **Run analysis**.

## Analysis modes

| Mode | Instrument input | Pairwise quantity minimized |
| --- | --- | --- |
| **Full spectra** | Editable laser wavelengths | Spectral overlap (SSC) after excitation scaling |
| **Conventional** | A detector table | Maximum predicted spillover in either direction, using optimized detector assignments |
| **Spectral** | A detector table | Cosine similarity of the full response vector across all detector channels |

Lower pairwise values are better. Optimization minimizes the worst pairwise value, then the sum of pairwise values to break ties.

## Instrument

For **Full spectra**, add or remove laser wavelengths in nm. No detector table is required.

For **Conventional**, the initial preset is BD FACSCanto II 8-color. **Spectral** starts with Cytek Aurora 5-laser / 64-channel. Choose another configuration from the preset menu or upload an instrument TSV. The selector and summary identify the active instrument; expand **Preview instrument table** to inspect it. Choose a bundled preset from the selector to replace an uploaded or restored table.

Presets are starting configurations: installed lasers and filters can differ, even on the same model. Check the exact configuration against your instrument. The FACSymphony 30-color entry is a published BD example for a configurable system. To customize a preset, download its TSV from the repository's [instrument folder](https://github.com/biosurf/spectracular/tree/master/instruments), edit the wavelength windows or detector rows, and upload it here. The [instrument guide](https://github.com/biosurf/spectracular/blob/master/instruments/README.md) records sources, configuration variants, and modelling assumptions, including the 900 nm upper limit used for otherwise unbounded long-pass channels.

### Creating a custom instrument TSV

Create a plain-text, tab-separated file with a header and one row per detector. A spreadsheet editor can export this as a `.tsv` file. Use these required columns:

```text
laser	channel	lower	upper
488	FITC-A	515	545
488	PE-A	564	606
635	APC-A	653	669
```

- **laser:** excitation wavelength in nm.
- **channel:** detector name, unique within each laser.
- **lower**, **upper:** detector wavelength bounds in nm; lower must be smaller than upper.

Use numeric wavelengths without units in the cells. Conventional and Spectral accept the same format. Additional columns are retained but are not required. Select the desired mode, upload the file under **Instrument**, and check the summary and preview before running.

## Panel settings

**Panel size** is the number of fluorochromes to select, with a minimum of two. The maximum follows **Available fluorochromes** and, for Conventional, the detector count. Select the fluorochromes available to your experiment; changing modes does not replace that selection.

**Minimum excitation ratio (MER, %)** is set in whole percentages and excludes fluorochromes whose excitation at every active laser is below MER. The live count shows how many pass this excitation filter. In Conventional, a fluorochrome-detector assignment must also meet MER for excitation at that detector's laser and peak intrinsic emission within its window. Spectral does not apply a per-channel MER cutoff.

### Panel constraints

Click **Edit constraints** to configure the following rules:

- **Must include:** select every fluorochrome that must appear in the panel, for example to extend an existing panel. Each uses its own position and cannot also supply a required group choice.
- **Cannot combine:** select two fluorochromes and click **Add pair**. The pair cannot appear together. Select a table row and click **Remove selected pair** to remove a restriction.
- **Groups:** select the members and the number of **Required choices**, then click **Add group**. Each requested choice must be a different fluorochrome. Other panel positions may include additional members of the same group. Zero choices imposes no requirement; it does not exclude the members.

For example, two choices from {FITC, PE, APC} and one choice from {FITC, APC} require three distinct fluorochromes, so a three-color panel must contain FITC, PE and APC. If FITC is also mandatory, it uses a separate position: the first group would need PE and APC, leaving no unused fluorochrome for the second group, so that request is infeasible. A mandatory FITC plus one choice from {FITC, PE} requires both FITC and PE. To edit a group, select its table row, change its settings, and click **Update selected group**. Use **Remove selected group** to delete it.

Click **Apply constraints** to use the edited rules. **Cancel** discards unapplied changes. **Clear all** clears the editor; apply to remove all constraints. Applying constraints clears previous results, so run the analysis again.

All requirements apply together. Overlapping groups require distinct choices, and each fluorochrome appears only once in the final panel. The number of mandatory fluorochromes plus the sum of group choices must not exceed the panel size. Remaining positions are unrestricted, so group requirements do not set exact membership counts in the completed panel. These requirements concern fluorochromes, not assignments to antibody markers.

Exclude unavailable fluorochromes through **Available fluorochromes**. Missing or excitation-filtered members are not silently substituted; feedback identifies required fluorochromes and groups with too few remaining candidates. Joint constraints and conventional detector eligibility can still make a panel impossible. Adjust the constraints, pool, panel size, MER or instrument when needed.

## Results

- **Results plot:** excitation/emission spectra, plus Conventional detector windows/signals or Spectral response signatures. Taller plots scroll; the heatmap is separate.
- **Selected panel:** selected fluorochromes; Conventional also shows their actual optimized detector assignments.
- **Pairwise matrix:** SSC for Full spectra, pairwise maximum predicted spillover for Conventional, or cosine similarity for Spectral. Conventional ratios can exceed 1: a value of 2.5 means 250% predicted spillover.
- **Excitation table:** excitation ratios at the active lasers for every selected fluorochrome. All rows are available without pagination; scroll as needed.

Detector responses are model-based estimates from reference spectra and detector windows. These pairwise metrics do not capture all biological constraints, measured compensation/spreading effects or whole-panel unmixing quality.

## Saving and restoring settings

**Save configuration** downloads an RDS file containing the active mode, available fluorochromes, panel size, MER, all applied constraints, and either the Full spectra laser vector or the complete detector table. Instrument names are labels; the saved table takes precedence over current preset definitions.

**Load configuration** restores those settings and shows restored detector tables as a saved entry in the selector. Previous results are cleared; click **Run analysis** to regenerate them. New files record the distinct-choice semantics. Older files without groups remain supported. Files with legacy required groups, exactly/at-least rules or unversioned groups are rejected without changing the current session: review and recreate those requirements using **Edit constraints**, then save a new configuration. They cannot be converted automatically because overlapping groups and mandatory members previously shared choices, and exact counts could exclude other members. Older files containing only the retired filter-builder settings require you to select a preset or upload an instrument TSV again.
