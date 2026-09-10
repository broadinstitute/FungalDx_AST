# =============================================================================
# 00_species_config.R
#
# The species/drug pairs this pipeline knows how to run.
#
# C. albicans / fluconazole is the worked demo in the README. Everything else
# here runs through exactly the same three scripts - the only thing that
# changes is one line at the top of 01, 02 and 03:
#
#     SPECIES_DRUG <- "aurisVori"
#
# Each entry below holds only what genuinely differs between pairs: the lane
# suffixes, any lane-label fixes, and the figure title. All paths are derived
# from the entry's name by get_species_config(), so a pair called "aurisVori"
# automatically reads:
#
#     data/raw/aurisVori/            raw nSolver exports
#     metadata/aurisVori_IDs.csv     strain sheet
#
# and writes data/normalized/aurisVori/, data/logfc/logfc_aurisVori_compiled.csv
# and figures/aurisVori_heatmap_SPR.{svg,pdf}.
#
# TO ADD A PAIR: drop the exports in data/raw/<name>/, write
# metadata/<name>_IDs.csv, add an entry below. Nothing else changes.
# =============================================================================

SPECIES_CONFIGS <- list(

  # --- Demo pair. 3 runs, 18 strains. See README. ---------------------------
  albicansFluc = list(
    label            = "C. albicans / fluconazole",
    untreated_suffix = "X4",
    treated_suffix   = "F4",
    # Run 01 labels one lane "087X4"/"087F4"; the sheet calls it "Fx087".
    sample_aliases   = c("087" = "Fx087"),
    fig_title        = "C. albicans Fluconazole log2 Fold Change - Heatmap + SPR"
  ),

  # --- 3 runs, 12 strains, all paired. --------------------------------------
  # Run 03 re-runs four strains and also carries eight lanes named
  # Fx06624 / Fx06644 etc. Those do not end in X4 or F4, so they are ignored
  # with a warning rather than silently mixed in. If they turn out to be a
  # condition you want, give them a suffix and add a pair for it.
  parapsilosisFluc = list(
    label            = "C. parapsilosis / fluconazole",
    untreated_suffix = "X4",
    treated_suffix   = "F4",
    sample_aliases   = c(),
    fig_title        = "C. parapsilosis Fluconazole log2 Fold Change - Heatmap + SPR"
  ),

  # --- 5 runs, 24 strains, two panel versions (C10397 and C4530). -----------
  # Runs 04/05 are a rerun pair: 04 is all untreated, 05 all treated, so they
  # only pair up after 01 compiles them together. Run 01 also carries 2 h and
  # "VQ" lanes, which are ignored (they do not end in X4/V4).
  aurisVori = list(
    label            = "C. auris / voriconazole",
    untreated_suffix = "X4",
    treated_suffix   = "V4",
    sample_aliases   = c(),
    # Two lanes are dead - every endogenous probe sits at background (max
    # count 1 and 5 respectively). Left in, they take every housekeeping
    # probe in their run down with them, because a normalizer failing in ANY
    # lane is dropped from ALL lanes.
    #   AAVQ4 (run 01) - an unused "VQ" condition, no downstream cost.
    #   BV4   (run 03) - strain B's TREATED lane. Excluding it leaves BX4
    #                    unpaired, so strain B drops out of the figure with a
    #                    warning from 02. There is no rerun of B in runs
    #                    04/05. Delete this entry if you would rather see the
    #                    strain and judge it yourself.
    exclude_lanes    = c("AAVQ4", "BV4"),
    fig_title        = "C. auris Voriconazole log2 Fold Change - Heatmap + SPR"
  ),

  # --- 4 runs, 24 of 26 strains (no raw data for U / 080). ------------------
  # Note the 1 h timepoint: lanes are X1 / M1, not X4 / M4.
  #
  # This sheet is maintained the other way round from the rest: the lane
  # codes (A, B, ... 076) sit in its Identifier column and the AR03xx strain
  # names in SampleID. Rather than rewrite the sheet - it is exported from
  # elsewhere and a fresh export would silently undo the edit - the two keys
  # below tell the pipeline which column is which.
  glabrataMica = list(
    label            = "C. glabrata / micafungin",
    untreated_suffix = "X1",
    treated_suffix   = "M1",
    sample_aliases   = c(),
    lane_key         = "Identifier",   # column the lane prefixes match
    label_key        = "SampleID",     # column that labels the figure
    fig_title        = "C. glabrata Micafungin log2 Fold Change - Heatmap + SPR"
  ),

  # --- Strain sheet only so far; drop exports in data/raw/glabrataFluc/. ----
  # Suffixes below are the lab default and unverified for this pair.
  glabrataFluc = list(
    label            = "C. glabrata / fluconazole",
    untreated_suffix = "X4",
    treated_suffix   = "F4",
    sample_aliases   = c(),
    fig_title        = "C. glabrata Fluconazole log2 Fold Change - Heatmap + SPR"
  ),

  # --- Strain sheet only so far; drop exports in data/raw/aurisMica/. -------
  # Suffixes unverified. C. glabrata micafungin used a 1 h timepoint (X1/M1),
  # so check the lane labels in your first export before trusting X4/M4.
  aurisMica = list(
    label            = "C. auris / micafungin",
    untreated_suffix = "X4",
    treated_suffix   = "M4",
    sample_aliases   = c(),
    fig_title        = "C. auris Micafungin log2 Fold Change - Heatmap + SPR"
  )
)


#' Fetch a species/drug config and fill in its derived paths.
#'
#' @param species_drug Name of an entry in SPECIES_CONFIGS.
#' @return The entry, plus every path the scripts need.
get_species_config <- function(species_drug) {

  if (!species_drug %in% names(SPECIES_CONFIGS)) {
    stop("Unknown SPECIES_DRUG '", species_drug, "'. Registered pairs: ",
         paste(names(SPECIES_CONFIGS), collapse = ", "),
         ".\n  Add an entry to R/00_species_config.R to register a new one.")
  }

  cfg <- SPECIES_CONFIGS[[species_drug]]
  cfg$key <- species_drug

  # ---- Derived paths -------------------------------------------------------
  cfg$raw_subdir             <- file.path("raw",        species_drug)
  cfg$normalized_subdir      <- file.path("normalized", species_drug)
  cfg$metadata_file          <- file.path("metadata", paste0(species_drug, "_IDs.csv"))
  cfg$compiled_rawdata_file  <- file.path("data", paste0("compiled_rawdata_",  species_drug, ".csv"))
  cfg$compiled_normdata_file <- file.path("data", paste0("compiled_normdata_", species_drug, ".csv"))
  cfg$control_qc_file        <- file.path("data", paste0("control_probe_QC_",  species_drug, ".csv"))
  cfg$duplicate_lanes_file   <- file.path("data", paste0("duplicate_lanes_",   species_drug, ".csv"))
  cfg$logfc_file             <- file.path("data/logfc", paste0("logfc_", species_drug, "_compiled.csv"))
  cfg$out_basename           <- paste0(species_drug, "_heatmap_SPR")

  # ---- Defaults for optional fields ---------------------------------------
  if (is.null(cfg$sample_aliases)) cfg$sample_aliases <- character(0)

  # Raw files to skip for this pair, e.g. c("parapsilosisFluc_run03.csv") to
  # leave a rerun out of the compiled table entirely.
  if (is.null(cfg$exclude_files))  cfg$exclude_files  <- character(0)

  # Individual lanes to drop before normalization, by lane label (e.g "BV4").
  # Use for dead or contaminated lanes. This matters more than it looks: a
  # housekeeping probe that fails in one lane is dropped from every lane in
  # that run, so a single dead lane can wipe out a whole run's normalizers.
  if (is.null(cfg$exclude_lanes))  cfg$exclude_lanes  <- character(0)

  # Minimum count for a housekeeping probe, applied after correction. A probe
  # below this in ANY lane is dropped from ALL lanes. Panels differ by more
  # than an order of magnitude in control-probe level, so this is per-pair.
  if (is.null(cfg$min_ctrl))       cfg$min_ctrl       <- 10

  # Which strain-sheet column the raw lane prefixes match, and which column
  # labels the figure. The usual arrangement is SampleID = internal lane code,
  # Identifier = publication name, and these defaults cover it. Override both
  # when a sheet is kept the other way round (see glabrataMica).
  if (is.null(cfg$lane_key))       cfg$lane_key       <- "SampleID"
  if (is.null(cfg$label_key))      cfg$label_key      <- "Identifier"

  cfg
}
