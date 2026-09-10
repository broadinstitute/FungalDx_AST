# =============================================================================
# 01_normalize_nanostring.R
#
# STEP 1 of the Fungal AST pipeline.
#
# Takes raw NanoString nCounter exports (nSolver "RCC collector" CSV format)
# and produces normalized count tables.
#
# NORMALIZATION ORDER (each step depends on the previous one):
#   1. Response-probe floor      - flag response probes reading below minResp
#   2. Negative control correct  - subtract per-lane background
#   3. Detection floor           - clamp values up to 6 x SD of negative controls
#   4. Positive control correct  - divide each lane by its ERCC scale factor
#   5. Control probe QC          - drop housekeeping probes reading below minCtrl
#   6. Control probe optimization- greedily drop the least stable housekeeping
#                                  probes until all CoVs are within limitCoV
#   7. Normalize                 - divide every probe by the geometric mean of
#                                  the surviving housekeeping probes
#
# INPUT   data/raw/*.csv         raw nSolver exports
# OUTPUT  data/normalized/<run>_normalized.csv   one per input file
#         data/compiled_rawdata.csv              all runs, raw counts
#         data/compiled_normdata.csv             all runs, normalized
#
# NEXT STEP  R/02_compute_logfc.R
#
# =============================================================================


# =============================================================================
# SECTION 0: LIBRARIES AND HELPERS
# =============================================================================

if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
library(tidyverse)

# Path to this repository's R/ directory.
# Set your working directory to the repo root before running, or edit this.
HELPERS_FILE <- "R/nanostring_helpers.R"

if (!file.exists(HELPERS_FILE)) {
  stop("Cannot find ", HELPERS_FILE, ".\n",
       "  Set your working directory to the repository root, e.g.\n",
       '  setwd("~/path/to/fungal-ast-pipeline")')
}
source(HELPERS_FILE)
source("R/00_species_config.R")


# =============================================================================
# SECTION 1: USER CONFIG - edit these for each dataset
# =============================================================================

# ---- Which species/drug pair -----------------------------------------------
# This is the only line that should normally change. See
# R/00_species_config.R for the full list of registered pairs and what to add
# when registering a new one.
SPECIES_DRUG <- "albicansFluc"
cfg <- get_species_config(SPECIES_DRUG)

# ---- Paths --------------------------------------------------------------
# Directory containing the analysis. Paths below are relative to the repo
# root; replace with absolute paths to run against data held elsewhere
# (e.g. a shared Google Drive folder).
DIR         <- "data/"                  # analysis root
SUBDIR      <- paste0(cfg$raw_subdir, "/")         # input:  raw nSolver CSVs live in DIR/SUBDIR
NEW_SUBDIR  <- paste0(cfg$normalized_subdir, "/")  # output: per-run normalized CSVs

# ---- Input format --------------------------------------------------------
# "nsolver" : standard nCounter export. 14 metadata rows, then probe rows;
#             column 2 holds the probe name, data starts at column 4.
# "matrix"  : an already-assembled probes x samples table, probe names in
#             column 1. Use this for hand-combined runs.
INPUT_FORMAT <- "nsolver"

# ---- Probe count thresholds ---------------------------------------------
# Applied to RAW counts, before any correction.

# Minimum for housekeeping / control probes ("_C_"). A control probe reading
# below this in ANY sample is dropped from ALL samples. Should sit comfortably
# above the negative control probes. Per-pair, because control-probe levels
# differ by more than an order of magnitude between panels - set it in
# R/00_species_config.R (`min_ctrl`), not here.
minCtrl <- cfg$min_ctrl

# Minimum for response probes. Set to 0 to disable (the default): response
# probes are informative even when low, and the 6x-SD floor already protects
# against dividing by noise.
minResp <- 0

# If TRUE, apply the response-probe threshold to the FIRST COLUMN only
# (the untreated sample), rather than requiring every sample to pass.
untrOnly <- FALSE

# ---- Control probe optimization -----------------------------------------
# Maximum tolerated coefficient of variation for a housekeeping probe.
# Probes are dropped worst-first until every survivor is within this limit.
limitCoV <- 0.25

# ---- Outlier handling (disabled by default) ------------------------------
# Used only by the optional outlier block in SECTION 4. Left here so the
# thresholds stay documented alongside everything else.
minRespUntr <- 0.1
minRespTrt  <- 0.1

# ---- Output verbosity ----------------------------------------------------
VERBOSE <- TRUE


# =============================================================================
# SECTION 2: SET UP
# =============================================================================

in_path  <- file.path(DIR, SUBDIR)
out_path <- file.path(DIR, NEW_SUBDIR)
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

# Input files, in the order list.files() returns them (alphabetical). Files
# are named <pair>_runNN.csv, so alphabetical order is run order.
dfnames <- list.files(path = in_path, pattern = "\\.csv$")

if (length(cfg$exclude_files) > 0) {
  skipped <- intersect(dfnames, cfg$exclude_files)
  if (length(skipped) > 0) {
    cat("Excluded by config: ", paste(skipped, collapse = ", "), "\n", sep = "")
    dfnames <- setdiff(dfnames, skipped)
  }
}

if (length(dfnames) == 0) {
  stop("No CSV files found in ", in_path)
}
cat("Found", length(dfnames), "input file(s):\n  ",
    paste(dfnames, collapse = "\n   "), "\n\n")

# Strip the directory and extension from a filename to get a run label.
run_label <- function(fname) gsub(fname, pattern = ".*/|.csv", replacement = "")

# ---- Accumulators --------------------------------------------------------
# QC record of which control probes passed, were optimized out, or failed.
controls <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(controls) <- c("probe", "experiment", "pass_or_fail")

rawdata_list        <- list()   # raw counts, per run
normdata_list       <- list()   # normalized counts, per run
no_control_norm_list <- list()  # corrected response probes, per run
control_norm_list   <- list()   # corrected control probes, per run


# =============================================================================
# SECTION 3: IMPORT RAW DATA
# =============================================================================

for (l in dfnames) {

  if (INPUT_FORMAT == "nsolver") {

    # Read with no header and locate the rows we need by their label in
    # column 1, rather than assuming fixed row numbers. Exports are not
    # consistent: most carry 15 metadata rows (File Name ... Messages), but
    # some are missing the File Name / Description rows entirely, which
    # silently shifted every fixed index in the original script.
    raw_lines <- read.csv(file.path(in_path, l), header = FALSE,
                          colClasses = "character", check.names = FALSE)
    labels <- trimws(raw_lines[[1]])

    # First probe row = first row whose column 1 is a probe class.
    first_probe <- which(labels %in% c("Positive", "Negative",
                                       "Endogenous", "Housekeeping"))[1]
    if (is.na(first_probe)) {
      stop("No probe rows found in ", l,
           " - column 1 should read Positive / Negative / Endogenous.")
    }

    # Column 2 holds the probe name; that becomes the row identifier and is
    # what every downstream grepl() on "_C_", "_R_", "NEG_", "POS_" relies on.
    probe_names <- trimws(raw_lines[[2]][first_probe:nrow(raw_lines)])

    # Columns 1-3 are Class / Name / Accession. Counts start at column 4.
    df <- raw_lines[first_probe:nrow(raw_lines), 4:ncol(raw_lines), drop = FALSE]
    df[] <- lapply(df, function(x) suppressWarnings(as.numeric(x)))

    # Drop trailing all-NA columns produced by ragged CSV export.
    keep <- colSums(is.na(df)) < nrow(df)
    df   <- df[, keep, drop = FALSE]
    rownames(df) <- probe_names

    # ---- Lane names ------------------------------------------------------
    # Normally the "Sample ID" row. Some older exports put bare lane numbers
    # (1, 2, 3, ...) there and the real condition labels in "Description".
    id_row   <- which(labels == "Sample ID")[1]
    desc_row <- which(labels == "Description")[1]
    if (is.na(id_row)) stop("No 'Sample ID' row found in ", l)

    lane_row_values <- function(i) {
      trimws(as.character(raw_lines[i, 4:ncol(raw_lines)]))[keep]
    }

    lane_names <- lane_row_values(id_row)

    if (all(grepl("^[0-9]+$", lane_names[nzchar(lane_names)])) &&
        !is.na(desc_row)) {
      alt <- lane_row_values(desc_row)
      if (any(nzchar(alt))) {
        warning(l, ": the Sample ID row holds bare lane numbers; ",
                "using the Description row for lane names instead.")
        lane_names <- alt
      }
    }

    # Only rename if the labels are unique. If the export carried duplicate
    # sample names, keep the original headers so nothing is silently merged.
    if (anyDuplicated(lane_names) == 0 && all(nzchar(lane_names))) {
      colnames(df) <- lane_names
    } else {
      warning("Duplicate or blank lane labels in ", l,
              " - keeping original column headers.")
    }

  } else if (INPUT_FORMAT == "matrix") {

    # Pre-combined table: probe names in column 1, counts thereafter.
    df_in <- read.csv(file.path(in_path, l))
    rownames(df_in) <- df_in[[1]]
    df <- df_in[, -1]

  } else {
    stop('INPUT_FORMAT must be "nsolver" or "matrix", got: ', INPUT_FORMAT)
  }

  # ---- Drop lanes listed in the config -----------------------------------
  if (length(cfg$exclude_lanes) > 0) {
    hit <- colnames(df) %in% cfg$exclude_lanes
    if (any(hit)) {
      cat("  ", l, ": excluding lane(s) ",
          paste(colnames(df)[hit], collapse = ", "),
          " (cfg$exclude_lanes)\n", sep = "")
      df <- df[, !hit, drop = FALSE]
    }
  }

  # ---- Dead lane check ----------------------------------------------------
  # A lane whose endogenous probes all sit at background is a failed
  # hybridization. It is worth catching here because such a lane will fail
  # every housekeeping probe, and a housekeeping probe that fails in ONE lane
  # is dropped from ALL of them - so one dead lane can wipe out a whole run's
  # normalizers and stop the script further down with a much less obvious
  # error. Add the lane to cfg$exclude_lanes to drop it.
  endo_rows <- !grepl("NEG_|POS_", rownames(df))
  if (any(endo_rows)) {
    lane_max <- suppressWarnings(apply(df[endo_rows, , drop = FALSE], 2,
                                       max, na.rm = TRUE))
    dead <- colnames(df)[is.finite(lane_max) & lane_max <= 10]
    if (length(dead) > 0) {
      warning(l, ": lane(s) ", paste(dead, collapse = ", "),
              " look dead - every endogenous probe is at background ",
              "(max count <= 10). They will fail the housekeeping probes for ",
              "this entire run. Add them to cfg$exclude_lanes if so.")
    }
  }

  rawdata_list[[run_label(l)]] <- df
}

cat("Imported", length(rawdata_list), "run(s).\n\n")


# =============================================================================
# SECTION 4: CORRECT AND NORMALIZE
# =============================================================================

for (i in dfnames) {

  cat("---- ", run_label(i), " ----\n", sep = "")

  df   <- rawdata_list[[run_label(i)]]
  cols <- length(colnames(df))

  # ---- Probe classes -----------------------------------------------------
  # Response probes: everything that is not a spike-in control and not a
  # housekeeping ("_C_") probe. If a panel carries genotypic or baseline
  # probes they land here and must be split out at the visualization stage.
  rspnsPrbs <- row.names(df)[!grepl("NEG_|POS_|_C_", row.names(df))]

  data_out <- data.frame(matrix(ncol = 0, nrow = length(rspnsPrbs)))
  row.names(data_out) <- rspnsPrbs

  # Everything that survives to the normalized output: response probes plus
  # housekeeping probes, but not the ERCC spike-ins.
  txPrbs <- row.names(df)[!grepl("NEG_|POS_", row.names(df))]

  data_norm <- data.frame(matrix(ncol = 0, nrow = length(txPrbs)))
  row.names(data_norm) <- txPrbs

  temp <- df

  # RESET point ------------------------------------------------------------
  # max_cov and removedCtrlProbes are per-file state and are reset here.
  max_cov <- 1
  removedCtrlProbes <- c()

  # ---- 1. Response probe floor -------------------------------------------
  # rm = FALSE: values below minResp become NA but the probe row is kept.
  allRespTrim <- rbind(
    temp[!rownames(temp) %in% rspnsPrbs, ],
    removeFailedProbes(temp[rspnsPrbs, ], minResp, FALSE, untrOnly)
  )

  # ---- 2. Negative control correction (background subtraction) -----------
  # Returns non-control probes only; NEG/POS rows are re-attached below.
  allDataNegCorr <- negCtrlCorr(allRespTrim)

  # ---- 3. Detection floor: 6 x SD of the negative controls ---------------
  # Anything below this is noise. Clamping up (rather than to zero or NA)
  # keeps the probe in play while preventing runaway fold changes.
  minValues <- minValue6xStDev(allRespTrim)

  for (col in 1:ncol(allDataNegCorr)) {
    allDataNegCorr[, col] <- ifelse(
      allDataNegCorr[, col] < minValues[col],
      minValues[col],
      allDataNegCorr[, col]
    )
  }

  # Re-attach the spike-in controls so posCtrlCorr() can find the POS probes.
  allDataNegCorr6xStDevMin <- rbind(
    allRespTrim[grepl("NEG_|POS_", row.names(allRespTrim)), ],
    allDataNegCorr
  )

  # ---- 4. Positive control correction (lane scaling) ---------------------
  allDataCorr <- posCtrlCorr(allDataNegCorr6xStDevMin)

  # ---- Split corrected data by probe class -------------------------------
  allCtrlCorr <- allDataCorr[grepl("_C_", rownames(allDataCorr)), ]
  allRespCorr <- allDataCorr[
    grepl("_R_|_Rb_|_Rdn_|_Ri_|_B_|_Rup_", rownames(allDataCorr)),
  ]

  # ---- 5. Control probe QC -----------------------------------------------
  # rm = TRUE: a housekeeping probe below minCtrl in any sample is dropped.
  allCtrlTrim <- removeFailedProbes(allCtrlCorr, minCtrl)

  # ---- 6. Control probe optimization -------------------------------------
  # Greedily remove the least stable housekeeping probe until every survivor
  # has a coefficient of variation within limitCoV. The final pass removes
  # nothing and simply exits the loop.
  ctrlProbes_toIter <- allCtrlTrim

  # A run with no surviving normalizers cannot be normalized at all. Say so
  # here, naming the likely cause, rather than failing inside the CoV loop
  # with an error that points nowhere near the real problem.
  if (nrow(allCtrlTrim) == 0) {
    worst <- names(which.min(apply(allCtrlCorr, 2, min, na.rm = TRUE)))
    stop("No housekeeping probes survived the minCtrl = ", minCtrl,
         " threshold in ", i, ".\n",
         "  Every one of the ", nrow(allCtrlCorr), " control probes reads ",
         "below it in at least one lane.\n",
         "  Lowest-reading lane: ", worst, " - if that lane is dead, add it ",
         "to cfg$exclude_lanes.\n",
         "  Otherwise lower cfg$min_ctrl for this pair in ",
         "R/00_species_config.R.")
  }
  if (nrow(allCtrlTrim) < 3) {
    warning(i, ": only ", nrow(allCtrlTrim), " housekeeping probe(s) passed ",
            "minCtrl = ", minCtrl, ". Fold changes from this run rest on a ",
            "thin normalization base - check data/control_probe_QC_*.csv.")
  }

  while (max_cov > limitCoV) {
    ctrlProbeCoV_list <- removeCtrlMaxCoV(ctrlProbes_toIter, limitCoV = limitCoV)

    ctrlProbes_toIter <- as.data.frame(ctrlProbeCoV_list[1])

    if (as.character(ctrlProbeCoV_list[2]) != "none") {
      removedCtrlProbes <- c(removedCtrlProbes, as.character(ctrlProbeCoV_list[2]))
    }

    max_cov <- max(unlist(ctrlProbeCoV_list[3]))
  }

  optCtrl <- ctrlProbes_toIter

  if (VERBOSE) {
    cat("  control probes: ", nrow(allCtrlCorr), " -> ", nrow(optCtrl), "\n", sep = "")
    if (length(removedCtrlProbes) > 0) {
      cat("  optimized out:  ", paste(removedCtrlProbes, collapse = ", "), "\n", sep = "")
    }
    cat("  final max CoV:  ", round(max_cov, 4), "\n", sep = "")
  }

  # ---- 7. Normalize to the optimized housekeeping set --------------------
  ctrlNorm    <- colGeoMeans(optCtrl)
  allRespNorm <- sweep(allRespCorr, 2, ctrlNorm, FUN = "/")

  # Normalize the housekeeping probes the same way, so they can be inspected
  # on the same scale (a well-behaved control should sit near 1).
  allCtrlNorm <- sweep(allCtrlTrim, 2, ctrlNorm, FUN = "/")

  # ---- Assemble the normalized table -------------------------------------
  # This merge is doing two jobs at once, and is kept verbatim from the
  # original because downstream code depends on both:
  #   (a) it reindexes onto the full txPrbs list, so probes dropped during
  #       control QC come back as all-NA rows with row name "NA", "NA.1", ...
  #       (02_compute_logfc.R filters these out);
  #   (b) merge() sorts by row name, so the output is alphabetical by probe.
  data_norm <- merge(
    data_norm[txPrbs, ],
    rbind(allCtrlNorm, allRespNorm)[txPrbs, ],
    by = "row.names", all = TRUE
  )
  rownames(data_norm) <- data_norm$Row.names
  data_norm <- data_norm[-1]

  # ---- QC record ---------------------------------------------------------
  failedCtrlProbes <- rownames(allCtrlCorr)[
    !(rownames(allCtrlCorr) %in% rownames(allCtrlTrim))
  ]

  # NOTE: the first column must be named `probe` to match the `controls`
  # accumulator. In the original script it was named `controls`, and because
  # rbind.data.frame() silently drops zero-row data frames before checking
  # names, the accumulator's header was discarded on the first bind rather
  # than raising an error. Fixed so the written QC file is labelled correctly.
  controls_temp <- data.frame(
    probe      = c(rownames(ctrlProbes_toIter), removedCtrlProbes, failedCtrlProbes),
    experiment = c(rep(gsub(i, pattern = "[.]csv", replacement = ""),
                       length(rownames(allCtrlCorr)))),
    pass_or_fail = c(
      rep("pass",    length(rownames(ctrlProbes_toIter))),
      rep("removed", length(removedCtrlProbes)),
      rep("failed",  length(failedCtrlProbes))
    )
  )
  controls <- rbind(controls, controls_temp)

  # ---- Store and write ---------------------------------------------------
  no_control_norm_list[[run_label(i)]] <- allRespCorr
  control_norm_list[[run_label(i)]]    <- allCtrlCorr

  # NEG_G / NEG_H are unused wells on some cartridge layouts.
  data_norm <- data_norm[!grepl("NEG_G|NEG_H", rownames(data_norm)), ]

  write.csv(
    x    = data_norm,
    file = file.path(out_path,
                     paste0(gsub(i, pattern = ".csv", replace = ""),
                            "_normalized.csv"))
  )

  # Keep only real probe rows in the in-memory copy used for the compiled table.
  data_norm <- data_norm[
    grepl("_R_|_Rb_|_Rdn_|_Ri_|_C_|_B_|_Rup_", rownames(data_norm)),
  ]
  normdata_list[[run_label(i)]] <- data_norm
}


# =============================================================================
# SECTION 5: COMPILE ALL RUNS
# =============================================================================
# NOTE: cbind() here will produce duplicate column names if the same sample ID
# appears in more than one run. Resolve by renaming samples upstream.

# Clearing the list names stops cbind() from prefixing columns with the run
# label, so a column stays "Fx041X4" rather than "CalbFluc001.Fx041X4".
#
# Outputs are named per species/drug pair (cfg$compiled_*_file) rather than
# the flat "data/compiled_normdata.csv" of the single-species original, so
# running this script for one pair never clobbers another pair's compiled
# tables.
# ---- Duplicate lane check ---------------------------------------------------
# The same strain is sometimes re-run on a later cartridge. cbind() will
# happily produce two columns with the same name, and everything downstream
# then silently uses whichever comes first. Report every repeat, keep the
# first occurrence (that is, the earliest run), and drop the rest.
lane_files <- unlist(lapply(names(normdata_list), function(run) {
  setNames(rep(run, ncol(normdata_list[[run]])), colnames(normdata_list[[run]]))
}))

dup_lanes <- unique(names(lane_files)[duplicated(names(lane_files))])

if (length(dup_lanes) > 0) {
  dup_report <- do.call(rbind, lapply(dup_lanes, function(lane) {
    runs <- unname(lane_files[names(lane_files) == lane])
    data.frame(lane = lane, runs = paste(runs, collapse = "; "),
               kept = runs[1], dropped = paste(runs[-1], collapse = "; "),
               stringsAsFactors = FALSE)
  }))
  write.csv(dup_report, cfg$duplicate_lanes_file, row.names = FALSE)

  warning(length(dup_lanes), " lane(s) appear in more than one run: ",
          paste(dup_lanes, collapse = ", "),
          ". Keeping the first (earliest) run of each; see ",
          cfg$duplicate_lanes_file,
          ". To use a different run, exclude the other file via ",
          "cfg$exclude_files in R/00_species_config.R.")

  cat("\n!! duplicate lanes across runs (kept the earliest):\n")
  print(dup_report, row.names = FALSE)
  cat("\n")
}

names(rawdata_list) <- NULL
raw_data <- do.call(cbind, rawdata_list)
raw_data <- raw_data[, !duplicated(colnames(raw_data)), drop = FALSE]
write.csv(x = raw_data, file = cfg$compiled_rawdata_file)

names(normdata_list) <- NULL
norm_data <- do.call(cbind, normdata_list)
norm_data <- norm_data[, !duplicated(colnames(norm_data)), drop = FALSE]
write.csv(x = norm_data, file = cfg$compiled_normdata_file)

# Control probe QC across all runs
write.csv(x = controls, file = cfg$control_qc_file, row.names = FALSE)

cat("\n=============================================================\n")
cat("Normalization complete. (", cfg$label, ")\n", sep = "")
cat("  per-run normalized CSVs : ", out_path, "\n", sep = "")
cat("  compiled raw counts     : ", cfg$compiled_rawdata_file, "\n", sep = "")
cat("  compiled normalized     : ", cfg$compiled_normdata_file, "\n", sep = "")
cat("  control probe QC        : ", cfg$control_qc_file, "\n", sep = "")
cat("Next step: R/02_compute_logfc.R (with the same SPECIES_DRUG)\n")
cat("=============================================================\n")


# =============================================================================
# APPENDIX: LEGACY / DISABLED BLOCKS
# =============================================================================
# Kept verbatim from the original script so older analyses remain
# reproducible. Neither block is active.
#
# ---- A. Negative control probe outlier removal --------------------------
# Intended to run once on composite data rather than per file. Removes
# negative probes whose mean or variance is more than 3x the average.
#
# Neg_Probes           <- df[grepl(rownames(df), pattern = "NEG"), ] %>% as.matrix()
# Neg_Probes_Mean      <- Neg_Probes %>% mean(na.rm = TRUE)
# Neg_Probes_Row_Mean  <- Neg_Probes %>% rowMeans()
# rm_Neg_Probes_Mean   <- row.names(Neg_Probes)[Neg_Probes_Row_Mean > 3 * Neg_Probes_Mean]
# Neg_Probes_Row_Vars  <- Neg_Probes %>% matrixStats::rowVars()
# rm_Neg_Probes_Vars   <- row.names(Neg_Probes)[Neg_Probes_Row_Vars > 3 * mean(Neg_Probes_Row_Vars)]
# rm_Neg_Probes        <- unique(c(rm_Neg_Probes_Mean, rm_Neg_Probes_Vars))
# if (length(rm_Neg_Probes) > 0) df[rm_Neg_Probes, ] <- NA
#
# ---- B. Post-normalization outlier handling ------------------------------
# Uses minRespUntr / minRespTrt from SECTION 1. Assumes column 1 is untreated
# and column 2 is treated, so it only works for two-column layouts.
#
# allRespNorm[, 1][allRespNorm[, 1] < 0.05] <- NA
# allRespNorm[, 2][(allRespNorm[, 1] < minRespUntr) & (allRespNorm[, 2] < minRespTrt)] <- NA
#
# ---- C. Baseline-probe variant -------------------------------------------
# For panels analysing baseline ("_B_") expression rather than drug response,
# the original swapped these two lines in:
#
# rspnsPrbs   <- row.names(df)[!grepl("NEG_|POS_|_C_", row.names(df))]
# allRespCorr <- allDataCorr[grepl("_B_", rownames(allDataCorr)), ]
# =============================================================================
