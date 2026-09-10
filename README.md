# Fungal AST — NanoString Transcriptional Antifungal Susceptibility Pipeline

A complete, runnable record of how we go from raw NanoString nCounter output to the published heatmap + SPR figure, demonstrated end to end on a *Candida albicans* / fluconazole dataset.

The premise: a susceptible fungal isolate mounts a large, stereotyped transcriptional response when it meets an antifungal. A resistant isolate barely reacts. Measuring that response with a small targeted probe panel gives a susceptibility call in hours rather than the days a growth-based MIC needs.

Everything in this repository operates on the data in `data/raw/`, so you can clone it and reproduce the figures without supplying anything of your own. The worked demo is *C. albicans* / fluconazole; three further species/drug pairs are included and run through the same scripts by changing a single line.

---

## Contents

- [Where this sits in the wider method](#where-this-sits-in-the-wider-method)
- [Quick start](#quick-start)
- [Species and drug pairs in this repository](#species-and-drug-pairs-in-this-repository)
- [Repository layout](#repository-layout)
- [How the scripts interact](#how-the-scripts-interact)
- [Script reference](#script-reference)
- [Probe naming convention](#probe-naming-convention)
- [The normalization method in detail](#the-normalization-method-in-detail)
- [Adapting to a different species or drug](#adapting-to-a-different-species-or-drug)
- [Requirements](#requirements)

---

## Where this sits in the wider method

The full experimental method has five stages. **This repository covers stages 3–5.** Stages 1–2 happen once per species/drug pair and produce the physical probe panel that stages 3–5 then read out.

| # | Stage | Where |
|---|-------|-------|
| 1 | **RNAseq on a small set of reference strains.** A handful of well-characterized susceptible and resistant isolates are grown ± drug and sequenced. | Upstream, not in this repo |
| 2 | **Gene selection.** The `geneSelect/` scripts rank genes by how strongly they separate drug-treated from untreated samples and susceptible from resistant strains. This step is run once during panel design. | `geneSelect/` |
| 3 | **Probe panel synthesis and hybridization.** A custom NanoString CodeSet is ordered against the selected genes, plus housekeeping normalizers and the standard ERCC spike-ins. Isolates are grown, RNA is hybridized to the CodeSet, and the nCounter system reads the counts. | Upstream, not in this repo |
| 4 | **Normalization and fold change.** Raw counts are corrected, normalized, and converted to per-strain log2 fold changes. | **This repo — `01`, `02`** |
| 5 | **Figure and susceptibility scoring.** The heatmap + SPR panel. | **This repo — `03`** |

> **Note on stage 2.** `GeneSelect` is deliberately not the main focus of this repo. It is used during panel design, consumes RNAseq counts rather than NanoString data, and is kept separate from the downstream normalization/plotting workflow.

---

## Quick start

```bash
git clone <this-repo>
cd fungal-ast-pipeline
```

Then, from **R with the working directory set to the repository root**:

```r
setwd("~/path/to/fungal-ast-pipeline")   # all paths in the scripts are relative to here

source("R/01_normalize_nanostring.R")    # data/raw/<pair>/  -> data/normalized/<pair>/
source("R/02_compute_logfc.R")           # compiled counts   -> data/logfc/
source("R/03_heatmap_spr.R")             # data/logfc/       -> figures/
```

Each script prints what it did and names the next one. Run them in order — each consumes the previous one's output.

Result: `figures/albicansFluc_heatmap_SPR.svg` (plus `.pdf` and a CSV of the SPR values for supplementary tables). A reference copy of both is committed, so you have something to diff against if your run differs.

**To run a different species/drug pair**, change one line at the top of each of the three scripts:

```r
SPECIES_DRUG <- "aurisVori"      # instead of "albicansFluc"
```

Everything else — input folder, strain sheet, lane suffixes, output names, figure title — follows from that one name via `R/00_species_config.R`.

---

## Species and drug pairs in this repository

*C. albicans* / fluconazole is the worked demo: it is what the numbers in this README and the committed reference figure refer to. The other pairs run through exactly the same three scripts.

| `SPECIES_DRUG` | Pair | Runs | Strains | Lanes | Status |
|---|---|---|---|---|---|
| `albicansFluc` | *C. albicans* / fluconazole | 3 | 18 / 18 | `X4` → `F4` | **Demo.** Reference figure committed. |
| `parapsilosisFluc` | *C. parapsilosis* / fluconazole | 3 | 12 / 12 | `X4` → `F4` | Runs end to end. |
| `aurisVori` | *C. auris* / voriconazole | 5 | 21 / 24 | `X4` → `V4` | **Blocked** on duplicate Identifiers — see data notes. |
| `glabrataMica` | *C. glabrata* / micafungin | 4 | 24 / 26 | `X1` → `M1` | Runs; 1 h timepoint. |
| `glabrataFluc` | *C. glabrata* / fluconazole | — | 0 / 26 | `X4` → `F4`? | Strain sheet only, awaiting exports. |
| `aurisMica` | *C. auris* / micafungin | — | 0 / 26 | `X4` → `M4`? | Strain sheet only, awaiting exports. |

Lane suffixes are not guessable — *C. glabrata* micafungin used a **1 hour** timepoint (`X1`/`M1`) where everything else used 4 hours. The two pairs awaiting data carry unverified suffixes, marked `?`; check them against the first export you drop in.

`data/raw/RUNS.csv` records, for every run, its original filename, the NanoString panel (RLF) it was read on, the date, and its lane labels.

### Data notes worth knowing before using these figures

- **Two *C. auris* strains share Identifiers with two others** in `aurisVori_IDs.csv` and `aurisMica_IDs.csv`: `X` and `AA` are both `1105`, `Y` and `P` are both `390`. Since the Identifier becomes the figure's column label, one of each pair would otherwise silently overwrite the other — which of the two survives depends on lane ordering, and the SPR values move with it. **`02` stops with an error until the sheets are fixed**, so there is currently no *C. auris* figure. This is deliberate: a quietly dropped strain is worse than a failed run.
- **Two dead lanes** in *C. auris* voriconazole, listed in that pair's `exclude_lanes`: `AAVQ4` (an unused condition) and `BV4` (strain B's treated lane, so strain B drops out; there is no rerun of it).
- **Thin normalization** on `parapsilosisFluc_run02` (3 surviving housekeeping probes) and `aurisVori_run02`/`run03` (4–5). Low FOV counts on those cartridges. Their fold changes are usable but rest on a narrow base.
- ***C. auris* spans two panel versions**, `CAUR_VORI_C10397` (runs 01–03) and `CAUR_VORI_2_C4530` (runs 04–05). Probe names match across both.
- **Lanes re-run on a later cartridge** are reported by `01`, which keeps the earliest and writes the rest to `data/duplicate_lanes_<pair>.csv`. Four lanes in *C. parapsilosis* and four in *C. auris* are affected.
- **Non-standard lanes are ignored, not merged**: the 2 h and `VQ` conditions in `aurisVori_run01`, and eight `Fx066 24 / 44`-style lanes in `parapsilosisFluc_run03` whose meaning is not recorded anywhere.

---

## Repository layout

```
fungal-ast-pipeline/
├── README.md
├── R/
│   ├── 00_species_config.R         Registry of species/drug pairs
│   ├── nanostring_helpers.R        Shared normalization functions
│   ├── 01_normalize_nanostring.R   Raw counts  -> normalized counts
│   ├── 02_compute_logfc.R          Normalized  -> log2 fold change
│   ├── 03_heatmap_spr.R            log2FC      -> publication figure
│   └── 04_gene_pathway_overlap.R   Genes recurring across pairs
├── geneSelect/                     Stage-2 panel design; not called by R/
├── metadata/
│   └── <pair>_IDs.csv              Strain MICs, susceptibility calls, RNAseq flags
├── data/
│   ├── raw/
│   │   ├── RUNS.csv                Provenance of every raw file (tracked)
│   │   └── <pair>/<pair>_runNN.csv nSolver exports (tracked)
│   ├── normalized/<pair>/          Generated by 01 (gitignored)
│   └── logfc/                      Generated by 02 (gitignored)
└── figures/                        Generated by 03 (reference copy tracked)
```

Raw exports are renamed to `<pair>_runNN.csv` on the way in, numbered chronologically; `RUNS.csv` maps each back to the filename it arrived with.

The raw data, strain sheets and one reference copy of the demo figure are tracked in git; everything else the scripts generate is ignored, so a fresh clone always reproduces outputs from source.

---

## GeneSelect note

The `geneSelect/` directory contains the gene-selection scripts used during the earlier panel-design step. 

---

## How the scripts interact

```
  data/raw/<pair>/*.csv               metadata/<pair>_IDs.csv
  (nSolver exports, 1 per cartridge)  (MIC, S/R call, RNAseq flag, Identifier)
         │                                          │
         │                                          │
         ▼                                          │
  ┌──────────────────────────┐                      │
  │ 01_normalize_nanostring  │◄── R/nanostring_helpers.R
  │                          │                      │
  │  neg ctrl → 6xSD floor   │                      │
  │  → pos ctrl → probe QC   │                      │
  │  → housekeeping optimize │                      │
  └──────────────────────────┘                      │
         │                                          │
         ▼                                          │
  data/normalized/<pair>/*_normalized.csv           │
  data/compiled_normdata_<pair>.csv                 │
         │                                          │
         ▼                                          │
  ┌──────────────────────────┐                      │
  │ 02_compute_logfc         │◄─────────────────────┤  (SampleID → Identifier)
  │                          │                      │
  │  log2(treated/untreated) │                      │
  └──────────────────────────┘                      │
         │                                          │
         ▼                                          │
  data/logfc/logfc_<pair>_compiled.csv              │
         │                                          │
         ▼                                          │
  ┌──────────────────────────┐                      │
  │ 03_heatmap_spr           │◄─────────────────────┘  (MIC, S/R, derivation strains)
  │                          │
  │  strip probe prefixes    │
  │  SPR projection          │
  │  order by MIC            │
  └──────────────────────────┘
         │
         ▼
  figures/<pair>_heatmap_SPR.svg / .pdf
  figures/<pair>_heatmap_SPR_table.csv
```

The metadata sheet is read by **both** `02` and `03`, for different reasons: `02` uses it only to rename columns from internal sample codes to publication identifiers; `03` uses it for MICs, susceptibility calls, and to work out which strains were the RNAseq derivation strains.

---

## Script reference

Every script has a clearly marked **`SECTION 1: USER CONFIG`** block at the top. That is the only part you should normally need to edit.

### `R/nanostring_helpers.R`

Sourced by `01`. Not run directly. Contains the normalization primitives:

| Function | Purpose |
|---|---|
| `colGeoMeans()` | Geometric mean per lane. Used everywhere counts are scaled, because nCounter variation is multiplicative. |
| `posCtrlCorr()` | Divides each lane by its ERCC positive-control scale factor. |
| `negCtrlCorr()` | Subtracts each lane's negative-control background. |
| `removeFailedProbes()` | Flags (and optionally drops) probes reading below a threshold in any sample. |
| `calcCtrlCoVs()` | Coefficient of variation of each housekeeping probe, relative to the lane geometric mean. |
| `removeCtrlMaxCoV()` | One step of greedy worst-first housekeeping probe pruning. |
| `minValue6xStDev()` | Per-lane detection floor at 6 × SD of the negative controls. |
| `negCtrl_wReplacement()`, `subLowValues()` | Alternative methods, unused by the current pipeline. Retained so older analyses stay reproducible. |

### `R/01_normalize_nanostring.R`

Reads every CSV in `data/raw/<pair>/` and writes one normalized CSV per run, plus compiled tables across all runs and a control-probe QC record.

Key config:

| Setting | Default | Meaning |
|---|---|---|
| `SPECIES_DRUG` | `"albicansFluc"` | Which pair to run. Everything below follows from it. |
| `DIR`, `SUBDIR`, `NEW_SUBDIR` | `data/`, `raw/<pair>/`, `normalized/<pair>/` | Input and output locations, from the config. Replace with absolute paths to run against data elsewhere. |
| `INPUT_FORMAT` | `"nsolver"` | `"nsolver"` for standard exports; `"matrix"` for a pre-assembled probes × samples table. |
| `minCtrl` | `cfg$min_ctrl` (10) | A housekeeping probe reading below this in **any** sample is dropped from **all** samples. Per-pair: control-probe levels differ by more than an order of magnitude between panels. |
| `minResp` | `0` | Same rule for response probes. `0` disables it. |
| `untrOnly` | `FALSE` | Apply the response threshold to the untreated lane only. |
| `limitCoV` | `0.25` | Housekeeping probes are pruned worst-first until every survivor is within its CoV. |

Outputs, all suffixed with the pair name: `data/normalized/<pair>/<run>_normalized.csv`, `data/compiled_rawdata_<pair>.csv`, `data/compiled_normdata_<pair>.csv`, `data/control_probe_QC_<pair>.csv`, and `data/duplicate_lanes_<pair>.csv` when a lane appears in more than one run. The suffix matters: without it, running a second pair would overwrite the first one's compiled tables.

> `control_probe_QC_<pair>.csv` records, per run, which housekeeping probes **passed**, were **removed** by CoV optimization, or **failed** the `minCtrl` threshold. Check it — if a run has very few surviving normalizers, its fold changes are built on a thin foundation. `01` warns below three and stops at zero.

`01` also warns about two things that are otherwise easy to miss: **dead lanes** (every endogenous probe at background — one such lane fails every housekeeping probe in its run, since a normalizer must be reliable in all of them) and **lanes appearing in more than one run**, where it keeps the earliest and records the rest.

### `R/02_compute_logfc.R`

Pairs each strain's treated lane with its untreated lane and takes `log2(treated / untreated)`.

| Setting | Default | Meaning |
|---|---|---|
| `SPECIES_DRUG` | `"albicansFluc"` | Must match the value used in `01`. |
| `UNTREATED_SUFFIX` / `TREATED_SUFFIX` | `"X4"` / `"F4"` | Lane naming convention, from the config. `X` = untreated, `F` = fluconazole, `4` = 4 hours. `V` = voriconazole, `M` = micafungin; the trailing digit is the timepoint. |
| `DROP_CONTROL_PROBES` | `TRUE` | Drop `_C_` housekeeping probes. Their fold change is ~0 by construction. |
| `SAMPLE_ALIASES` | `cfg$sample_aliases` | Explicit fixes where a lane label doesn't match the metadata `SampleID`. |
| `LANE_KEY` / `LABEL_KEY` | `"SampleID"` / `"Identifier"` | Which strain-sheet column the lane prefixes match, and which one names the output columns. |
| `NA_REPLACEMENT` | `0.1` | Floor for missing normalized values, so the log2 is defined. |

Infinite fold changes (zero denominator) are set to `0` — treated as "no measurable change" rather than dropping the strain.

> **On `SAMPLE_ALIASES`:** in the demo data, run 01 labels one lane `087X4`/`087F4` while the strain sheet calls that strain `Fx087`. Without the alias, that strain is silently dropped. The script warns loudly about any unpaired or unmapped sample, so watch the console.

`02` also stops if two strains share one `Identifier`. Output columns are named by it, so a shared value means the second assignment overwrites the first and a strain disappears from the figure with nothing in the console to show for it.

### `R/03_heatmap_spr.R`

Builds the figure. This is the single-run version: one measurement per strain, plain dots, no error bars.

| Setting | Default | Meaning |
|---|---|---|
| `SPECIES_DRUG` | `"albicansFluc"` | Must match the value used in `01` and `02`. Sets the input, output basename and figure title. |
| `STRIP_PROBE_PREFIX` | `TRUE` | Turn `CaFluc3_R_ERG11` into `ERG11` on the row labels. Matched lazily, so two-token panel names (`Caur_Vori_R_...`) strip correctly. |
| `COL_ORDER` | `"MIC_asc"` | Most susceptible on the left. Ties broken by descending SPR. |
| `SUS_COLORS` | black / grey | Susceptibility strip colours. Keys must match your metadata values. |
| `FIG_WIDTH` / `FIG_HEIGHT` | `10` / `8` in | Widen for many strains. |

---

## Probe naming convention

Everything downstream classifies probes by pattern-matching the probe name, so **the naming convention is load-bearing**. Get it wrong when designing a panel and probes will be silently misclassified.

Format: `<panel>_<class>_<gene>`

| Pattern | Class | Role |
|---|---|---|
| `POS_A` … `POS_F` | ERCC positive spike-in | Lane-to-lane scaling. Known concentrations. |
| `NEG_A` … `NEG_F` | ERCC negative spike-in | Background estimate and detection floor. |
| `<panel>_C_<gene>` | Housekeeping / normalizer | e.g. `CaFluc3_C_ACT1`. Used to build the normalization factor; **excluded** from fold-change output. |
| `<panel>_R_<gene>` | Response | e.g. `CaFluc3_R_ERG11`. The genes selected in stage 2. |
| `<panel>_Rup_<gene>` | Response, expected up | e.g. `CaFluc1_Rup_CDR1`. Treated identically to `_R_` throughout. |
| `<panel>_Rdn_`, `_Rb_`, `_Ri_` | Other response subclasses | Also treated as response probes. |
| `<panel>_B_<gene>` | Baseline expression | Constitutive rather than drug-responsive. Grouped with response probes by default; see Appendix C of `01` for the baseline-only variant. |

`NEG_G` / `NEG_H` are unused wells on some cartridge layouts and are dropped.

`03` strips `<panel>_<class>_` for display, so `CaFluc3_R_ERG11` is labelled `ERG11`. Genes without a common name keep their systematic ID — `CaFluc1_R_CAALFM_C301540WA` becomes `CAALFM_C301540WA`.

---

## The normalization method in detail

`01` applies seven steps **in this order**. The order matters: background is subtracted before lane scaling, because the background is additive and the lane effect is multiplicative.

1. **Response probe floor.** Response probes reading below `minResp` are flagged `NA` but kept (`rm = FALSE`). Default `minResp = 0` makes this a no-op; step 3 does the real protection.

2. **Negative control correction.** The mean of the six ERCC negative probes in a lane estimates its non-specific binding. That value is subtracted from every non-control probe in that lane.

3. **Detection floor at 6 × SD of the negative controls.** Anything below this is indistinguishable from noise. Values are clamped **up** to the floor rather than zeroed or dropped, which keeps the probe usable while stopping a near-zero denominator from producing an absurd fold change later.

4. **Positive control correction.** The geometric mean of the six ERCC positive spike-ins in a lane is divided by the average across lanes to give a per-lane scale factor. Every probe in the lane is divided by it, correcting for hybridization efficiency, binding density, and scan quality.

5. **Housekeeping probe QC.** A `_C_` probe reading below `minCtrl` in *any* sample is dropped from *all* samples — a normalizer has to be reliable everywhere or it is not a normalizer.

6. **Housekeeping probe optimization.** Each surviving `_C_` probe is expressed as a ratio to the geometric mean of all `_C_` probes in its lane, which removes the lane effect; the remaining spread is that probe's own instability. The worst probe is dropped, the CoVs recomputed, and the process repeats until every survivor is within `limitCoV` (0.25). Removals are logged to `control_probe_QC.csv`.

7. **Normalize.** Every probe is divided by the geometric mean of the optimized housekeeping set for its lane. Housekeeping probes are normalized the same way, so a well-behaved one sits near 1 and can be sanity-checked.



## Adapting to a different species or drug

Nothing in `01`–`03` is specific to a species. Adding a pair means adding an entry to `R/00_species_config.R` and putting two things where the config expects them.

1. **Drop the nSolver exports** into `data/raw/<pair>/`, named `<pair>_runNN.csv` in run order, and add a row per file to `data/raw/RUNS.csv`.
2. **Write a strain sheet** at `metadata/<pair>_IDs.csv` with the five required columns — `SampleID`, `MIC`, `Susceptibility`, `RNAseq`, `Identifier` — marking derivation strains in `RNAseq`. `SampleID` must match the lane label prefix in the data; `Identifier` is the publication name and must be unique. Extra columns are ignored, so a sheet can carry its own annotations.
3. **Add the entry** to `SPECIES_CONFIGS`:

   ```r
   aurisMica = list(
     label            = "C. auris / micafungin",
     untreated_suffix = "X4",
     treated_suffix   = "M4",
     sample_aliases   = c(),
     fig_title        = "C. auris Micafungin log2 Fold Change - Heatmap + SPR"
   )
   ```

   Every path is derived from the entry's name, so `aurisMica` automatically reads `data/raw/aurisMica/` and `metadata/aurisMica_IDs.csv` and writes `figures/aurisMica_heatmap_SPR.svg`.
4. **Set `SPECIES_DRUG <- "aurisMica"`** at the top of `01`, `02` and `03`, and run them in order.

Four things are worth checking on a first run, because each fails quietly:

- **Lane suffixes.** Read them off the `Sample ID` row of an export rather than assuming. The micafungin panel here used `X1`/`M1`, not `X4`/`M4`.
- **`min_ctrl`.** Control-probe levels differ by more than an order of magnitude between panels, so the default of 10 is not universal. `01` stops with an explanation if nothing survives it.
- **Probe naming.** `<panel>_<class>_<gene>`, per the convention above — the most common source of silent misclassification.
- **Unique `Identifier`s.** Two strains sharing one will stop `02`.

Optional per-pair settings: `sample_aliases` for lane labels that don't match the sheet, `exclude_lanes` for dead lanes, `exclude_files` to leave a whole run out, and `min_ctrl`.

> **If a sheet is arranged the other way round** — lane codes under `Identifier`, strain names under `SampleID` — don't rewrite it. Set `lane_key` and `label_key` in that pair's config entry and the pipeline reads it as it stands, which survives a fresh export from wherever the sheet is maintained. `glabrataMica` is set up this way.

---

## Requirements

R ≥ 4.0. The scripts install anything missing on first run.

| Package | Source | Used by |
|---|---|---|
| `tidyverse` | CRAN | `01` |
| `dplyr` | CRAN | `02`, `03` |
| `ComplexHeatmap` | Bioconductor | `03` |
| `circlize` | Bioconductor | `03` |
| `svglite` | CRAN | `03` |

---

## Demo dataset

`data/raw/albicansFluc/` contains three *C. albicans* fluconazole runs on the `CAFLUC3_2_C4732` panel — 18 strains, each with a paired untreated (`X4`) and fluconazole-treated (`F4`) lane at 4 hours, spanning MICs from 0.0625 to 256 µg/mL. Six strains went through RNAseq and define the signature; the remaining twelve are validation strains.

The raw data for the other pairs sits alongside it under `data/raw/`, one folder per pair, with provenance in `data/raw/RUNS.csv`.
