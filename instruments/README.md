# Instrument configurations

These TSVs are editable starting configurations for Spectracular, not a guarantee of the optics installed on every instrument bearing the same model name. Check the laser wavelengths, enabled channels, filters and dichroic layout against the actual instrument before use. Defaults and published example layouts are both useful starting points; the labels distinguish them.

Each row is one fluorescence measurement channel, with wavelengths in nm:

```text
laser	channel	lower	upper
488	FITC	515	545
```

Download a TSV, edit the rows or wavelength bounds, and select **Upload instrument settings** in Shiny. In the standalone script, read it with `read.delim("instruments/<filename>.tsv")`. The shared backend uses `laser`, `channel`, `lower` and `upper`; extra columns in the original presets are retained. Channels must be unique within each laser. Labels such as V1/B1 or B527 identify the measurement window; they need not match the names used in an exported FCS file.

## Wavelength conventions

- A nominal bandpass `center/width` becomes `center - width/2` to `center + width/2`. Half-nanometre boundaries are retained.
- For the new BD LSRFortessa, FACSymphony 30-color and A5 SE tables, bandpasses are intersected with the published dichroic path: the channel's long-pass cutoff and the next higher-wavelength branch's cutoff. For example, Fortessa R_B has a 730/45 filter but a 710–750 nm path, so its stored window is 710–750 nm. A5 SE B750 has a 750/60 filter and a 724–770 nm path, giving 724–770 nm.
- The Attune and NovoCyte layouts likewise include the published routing limits where these restrict a bandpass. Attune VL1 is 417–465 nm after its 417 LP blocker. Where the source supplies an explicit detector range (Cytek and several MACSQuant channels), that range is used directly.
- Otherwise, the listed nominal filter bounds are used. These rectangular windows approximate optics; they do not encode filter transmission curves, laser power, detector sensitivity or gains.
- **Unbounded long-pass channels use 900 nm as an explicit modelling endpoint**, following the existing FACSCalibur example. This is not a manufacturer-specified detector cutoff. It affects Accuri FL3, FACSLyric V755LP, MACSQuant Analyzer 10 B4/R2, MACSQuant VYB Y4, and Navios/EX FL5/FL8. Adjust it if the instrument's effective response is known. The shared integration routine also restricts integration to the available spectrum.
- FSC, SSC and imaging channels are excluded. Fluorescence channel counts can therefore differ from an instrument's advertised total parameter count.
- The three original TSVs are unchanged. New Cytek subsets retain the original Aurora's explicitly documented start/end boundaries; they are not reconstructed from rounded centre wavelengths.

## Conventional configurations

Sources checked 2026-11-09. Each filename links to the actual TSV. B = blue, V = violet, Y/YG = yellow-green, R = red.

| TSV | Lasers (nm) | Fluorescence channels | Basis and source |
| --- | --- | ---: | --- |
| [BD FACSCanto II](bd_facscantoII_8color.tsv) | 405, 488, 633 | 8 | Existing 4–2–2 preset; [BD filter guide][canto]. |
| [BD FACSCalibur](bd_facscalibur_4color.tsv) | 488, 635 | 4 | Existing preset; [BD brochure][calibur]. FL3 uses the 900 nm modelling endpoint. |
| [BD Accuri C6 Plus](bd_accuri_c6plus_4color.tsv) | 488, 640 | 4 | Standard installed filter set, not optional replacements; [technical specifications][accuri]. |
| [BD FACSLyric](bd_facslyric_12color.tsv) | 405, 488, 640 | 12 | Full 12-color column in the [BD filter guide][lyric]. The 4/6/8/10-color versions are not represented by this file. |
| [BD LSRFortessa UV/V/B/R](bd_lsrfortessa_4laser_UV_18color.tsv) | 355, 405, 488, 640 | 18 | Manual base configuration: UV2/V6/B6/R4; [user guide][fortessa], printed pp. 136, 141, 144, 147–148. This is the UV configuration, not a 561 nm-equipped Fortessa. |
| [Beckman CytoFLEX B4-R3-V4](beckman_cytoflex_B4_R3_V4_11color.tsv) | 405, 488, 638 | 11 | B4/R3/V4 configuration and filters in the [manufacturer configurator][cytoflex]. |
| [Beckman CytoFLEX S V4-B2-Y4-R3](beckman_cytoflex_s_V4_B2_Y4_R3_13color.tsv) | 405, 488, 561, 638 | 13 | V-B-Y-R series in [Beckman Figure 1][cytoflex_srt]. This version has 660/10 filters on V and R; other installed variants differ. |
| [Beckman CytoFLEX SRT V5-B2-Y5-R3](beckman_cytoflex_srt_V5_B2_Y5_R3_15color.tsv) | 405, 488, 561, 638 | 15 | Full SRT layout in [Beckman Figure 1][cytoflex_srt]; includes the additional violet and yellow-green channels. |
| [Thermo Attune NxT V4-B3-Y4-R3](thermo_attune_nxt_V4_B3_Y4_R3_14color.tsv) | 405, 488, 561, 637 | 14 | [Rockefeller's dated 2020 default optical layout][attune], both pages; filter specifications cross-checked against the [manufacturer filter catalogue][attune_filters]. Not the newer V6 configuration. |
| [Miltenyi MACSQuant Analyzer 10](miltenyi_macsquant_analyzer10_8color.tsv) | 405, 488, 640 | 8 | [Miltenyi manual][macsquant], Table 9.5. Ten total parameters include two scatter channels. |
| [Miltenyi MACSQuant Analyzer 16](miltenyi_macsquant_analyzer16_14color.tsv) | 405, 488, 640 | 14 | [Miltenyi manual][macsquant], Table 9.5; V5/B6/R3 fluorescence channels. |
| [Miltenyi MACSQuant VYB](miltenyi_macsquant_vyb_8color.tsv) | 405, 488, 561 | 8 | [Miltenyi manual][macsquant], Table 9.5; V2/B2/Y4. |
| [Beckman Navios](beckman_navios_10color.tsv) | 405, 488, 638 | 10 | [FDA K183592][navios], optical-filter comparison, with bandwidths from [Beckman's filter list][navios_filters]. FL3 is 620/30; FL4 uses 695/30 rather than the alternative 675/20. |
| [Beckman Navios EX](beckman_navios_ex_10color.tsv) | 405, 488, 638 | 10 | Same documented layout except FL3 is 614/20, as specified in [FDA K183592][navios]. |
| [Agilent NovoCyte 3000](agilent_novocyte_3000_13color.tsv) | 405, 488, 640 | 13 | Default V/B/R configuration in [Agilent guide][novocyte], printed p. 3. Six physical PMTs provide 13 laser-resolved fluorescence channels. Row labels retain laser plus detector number; e.g. V5, B5 and R5 are distinct measurements. Not a Quanteon configuration. |
| [BD FACSymphony 30-color example](bd_facsymphony_5laser_30color_example.tsv) | 355, 405, 488, 561, 637 | 30 | Published five-laser example in [BD configuration guide][symphony], p. 2: UV8/V8/B6/YG5/R3, excluding blue SSC. A starting template for configurable FACSymphony systems, **not a universal A3/A5/S6 factory default**. |

## Spectral configurations

| TSV | Lasers (nm) | Fluorescence channels | Basis and source |
| --- | --- | ---: | --- |
| [Cytek Aurora 5-laser](cytek_aurora_5laser_64channel.tsv) | 355, 405, 488, 561, 640 | 64 | Existing preset; [Aurora guide][aurora], printed pp. 75–76. |
| [Cytek Aurora V/B/R](cytek_aurora_3laser_VBR_38channel.tsv) | 405, 488, 640 | 38 | Corresponding V16/B14/R8 arrays from the [Aurora detector table][aurora]. |
| [Cytek Aurora V/B/YG/R](cytek_aurora_4laser_VBYR_48channel.tsv) | 405, 488, 561, 640 | 48 | Corresponding V16/B14/YG10/R8 arrays from the [Aurora detector table][aurora]; the fourth laser is 561 nm, not UV. |
| [Cytek Northern Lights](cytek_northern_lights_3laser_38channel.tsv) | 405, 488, 640 | 38 | V16/B14/R8 arrays in the [Northern Lights guide][northern], optical configuration and filter-bandwidth tables. |
| [Cytek Aurora CS V/B/R](cytek_aurora_cs_3laser_VBR_38channel.tsv) | 405, 488, 640 | 38 | V16/B14/R8 arrays in [Aurora CS guide][aurora_cs], printed pp. 80–82. |
| [Cytek Aurora CS V/B/YG/R](cytek_aurora_cs_4laser_VBYR_48channel.tsv) | 405, 488, 561, 640 | 48 | V16/B14/YG10/R8 arrays in [Aurora CS guide][aurora_cs]; fourth laser is 561 nm. |
| [Cytek Aurora CS 5-laser](cytek_aurora_cs_5laser_64channel.tsv) | 355, 405, 488, 561, 640 | 64 | Full detector table in [Aurora CS guide][aurora_cs], printed pp. 80–82. |
| [BD FACSymphony A5 SE, UV 355 nm](bd_facsymphony_a5se_355nm_48channel.tsv) | 355, 405, 488, 561, 637 | 48 | [A5 SE guide][a5se], printed pp. 50–53, UV10/V14/B9/YG9/R6; parameter names preserved. Dichroic limits included. |
| [BD FACSymphony A5 SE, UV 349 nm](bd_facsymphony_a5se_349nm_48channel.tsv) | 349, 405, 488, 561, 637 | 48 | Alternative UV wavelength explicitly documented in the [A5 SE guide][a5se], printed p. 51. Detector windows are unchanged; excitation wavelength differs. |

## Not included in this batch

Sony ID7000/SA3800, BD FACSDiscover S8/A8 and Beckman CytoFLEX mosaic still need complete numerical per-channel boundaries confirmed before creating a TSV. An overall wavelength range and a detector count are not enough to infer evenly spaced channels. The Agilent NovoCyte Opteon remains a candidate pending a complete check of its detailed core-facility table against the manufacturer's configuration; no missing bands have been filled by inference.

[canto]: https://www.bdbiosciences.com/content/dam/bdb/marketing-documents/BD_FACSCanto_II_FilterGuide.pdf
[calibur]: https://www.bdbiosciences.com/content/dam/bdb/marketing-documents/eu/23-3090-04_facscalibur_brochure.pdf
[accuri]: https://www.bdbiosciences.com/content/dam/bdb/marketing-documents/BD-Accuri-C6-Plus-Technical-Specs.pdf
[lyric]: https://www.bdbiosciences.com/content/dam/bdb/marketing-documents/products-pdf-folder/instruments/clinical-cell-analyzers/BD-Clinical-FACSLyric-Filter-QR-Guide-%28IVD%29.pdf
[fortessa]: https://www.bdbiosciences.com/content/dam/bdb/marketing-documents/BD_LSRFortessa_cell_analyzer_user_guide.pdf
[cytoflex]: https://www.beckman.com/assets/applications/flowcytometry/cytoflexquoter/cytoflexq.html
[cytoflex_srt]: https://media.beckman.com/-/media/pdf-assets/application-notes/flow-cytometry-application-note-migration-of-panels-designed-on-the-cytoflex-s-to-cytoflex-srt.pdf?rev=a9d94b6909764b4c9548fbfda1005f35
[attune]: https://www.rockefeller.edu/fcrc/uploads/www.rockefeller.edu/sites/147/2021/01/Default-Filters-for-Configuration-ThermoFisher_Attune-NxT_14-COLOR_20201216.pdf
[attune_filters]: https://www.thermofisher.com/order/catalog/product/100022763S
[macsquant]: https://static.miltenyibiotec.com/asset/150655405641/document_nvailh5sa963jcjmuhfelstk41?content-disposition=inline
[navios]: https://www.accessdata.fda.gov/cdrh_docs/reviews/K183592.pdf
[navios_filters]: https://www.beckman.com/support/faq/product/navios-optical-filters
[novocyte]: https://www.agilent.com/cs/library/quickreference/public/qr-guide-flow-cytometer-configuration-guide-novocyte-5994-4185en-agilent.pdf
[symphony]: https://www.bdbiosciences.com/content/dam/bdb/marketing-documents/BD_FACSymphony_Configuration_Guide.pdf
[aurora]: https://welcome.cytekbio.com/hubfs/N9-20006%20Rev.%20E_Cytek%20Aurora%20User%20Guide.pdf
[northern]: https://welcome.cytekbio.com/hubfs/Website%20Downloadable%20Content/Guides%20and%20Manuals/N9-20039US%20Northern%20Lights%20UG.pdf
[aurora_cs]: https://welcome.cytekbio.com/hubfs/Website%20Downloadable%20Content/Guides%20and%20Manuals/N9-20077%20Aurora%20CS%20UG.pdf
[a5se]: https://www.bdbiosciences.com/content/dam/bdb/marketing-documents/products-pdf-folder/instruments/research-cell-analyzers/bd-facsymphony-a5-se-user-guide-ruo-global.pdf
