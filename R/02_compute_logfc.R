# =============================================================================
# 02_compute_logfc.R
#
# STEP 2 of the Fungal AST pipeline.
#
# Converts the normalized count table from step 1 into a log2 fold-change
# matrix: for each strain, the drug-treated lane divided by its paired
# untreated lane.
#
#   log2FC[gene, strain] = log2( normalized[gene, <strain>F4]
#                                / normalized[gene, <strain>X4] )
#
# Columns are renamed from internal sample codes (Fx041) to the publication
# identifiers used in the figures (P75016), using the metadata sheet.
#
# INPUT   data/compiled_normdata.csv        from 01_normalize_nanostring.R
#         metadata/albicansFluc_IDs.csv     SampleID -> Identifier mapping
# OUTPUT  data/logfc/logfc_albicansFluc_compiled.csv
#
# NEXT STEP  R/03_heatmap_spr.R
# =============================================================================


# =============================================================================
# SECTION 0: LIBRARIES
# =============================================================================

if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(dplyr)


# =============================================================================
# SECTION 1: USER CONFIG
# =============================================================================

# ---- Paths ---------------------------------------------------------------
NORM_FILE <- "data/compiled_normdata.csv"
IDS_FILE  <- "metadata/albicansFluc_IDs.csv"
OUT_FILE  <- "data/logfc/logfc_albicansFluc_compiled.csv"

# ---- Sample naming convention -------------------------------------------
# Lane columns are named <SampleID><suffix>, e.g. Fx041X4 / Fx041F4.
#   X4 = untreated, 4 hours
#   F4 = fluconazole treated, 4 hours
# Change these for a different drug or timepoint (e.g. "M4" for micafungin).
UNTREATED_SUFFIX <- "X4"
TREATED_SUFFIX   <- "F4"

# ---- Probe filtering -----------------------------------------------------
# Drop housekeeping probes ("_C_") before computing fold changes. They were
# used to normalize, so their fold change is uninformative and near zero by
# construction. Set FALSE to keep them for QC inspection.
DROP_CONTROL_PROBES <- TRUE

# ---- Sample code aliases -------------------------------------------------
# Occasionally a lane is labelled in nSolver with a shortened sample code that
# does not match the SampleID in the metadata sheet. Map those here rather
# than guessing with a regex, so the correction is explicit and reviewable.
#
# In the demo data, run 001 labels one lane "087X4"/"087F4" while the metadata
# sheet calls that strain "Fx087". Without this alias the strain is silently
# dropped from the figure.
SAMPLE_ALIASES <- c(
  "087" = "Fx087"
)

# ---- Zero / missing handling --------------------------------------------
# Normalized values that are NA become this floor before the ratio is taken.
# It must be > 0 or the log2 is undefined.
NA_REPLACEMENT <- 0.1


# =============================================================================
# SECTION 2: LOAD
# =============================================================================

sample_ids <- read.csv(IDS_FILE, stringsAsFactors = FALSE)

stopifnot(all(c("SampleID", "Identifier") %in% colnames(sample_ids)))

# Lookup: internal sample code -> publication identifier
# e.g. "Fx041" -> "P75016"
id_lookup <- setNames(as.character(sample_ids$Identifier), sample_ids$SampleID)

norm_data <- read.csv(NORM_FILE, check.names = FALSE, stringsAsFactors = FALSE)

# Step 1's merge() reindexes onto the full probe list, so probes that were
# dropped during control QC come back as all-NA rows with row name "NA".
# Remove those before doing anything else.
norm_data <- norm_data[!is.na(norm_data[, 1]) &
                         norm_data[, 1] != "NA" &
                         norm_data[, 1] != "", ]

# Keep the first occurrence of any probe that appears more than once.
norm_data <- norm_data[!duplicated(norm_data[, 1]), ]

rownames(norm_data) <- norm_data[, 1]
norm_data <- norm_data[, -1, drop = FALSE]

if (DROP_CONTROL_PROBES) {
  n_before  <- nrow(norm_data)
  norm_data <- norm_data[!grepl("_C_", rownames(norm_data)), , drop = FALSE]
  cat("Dropped", n_before - nrow(norm_data), "housekeeping probes;",
      nrow(norm_data), "response probes remain.\n")
}

norm_data[is.na(norm_data)] <- NA_REPLACEMENT


# =============================================================================
# SECTION 3: PAIR UNTREATED AND TREATED LANES
# =============================================================================

untreated_cols <- grep(paste0(UNTREATED_SUFFIX, "$"), colnames(norm_data), value = TRUE)
treated_cols   <- grep(paste0(TREATED_SUFFIX,   "$"), colnames(norm_data), value = TRUE)

# Strip the suffix to recover the sample code, then keep only codes that have
# BOTH an untreated and a treated lane.
prefixes_untreated <- sub(paste0(UNTREATED_SUFFIX, "$"), "", untreated_cols)
prefixes_treated   <- sub(paste0(TREATED_SUFFIX,   "$"), "", treated_cols)
shared_prefixes    <- intersect(prefixes_untreated, prefixes_treated)

unpaired <- setdiff(union(prefixes_untreated, prefixes_treated), shared_prefixes)
if (length(unpaired) > 0) {
  warning("Samples missing a treated/untreated partner, skipped: ",
          paste(unpaired, collapse = ", "))
}

# Keep only sample codes present in the metadata sheet.
#   ColPrefix  = the code as it appears in the data columns
#   SampleID   = the code as it appears in the metadata sheet (after aliasing)
#   Identifier = the publication name used in the figure
pair_map <- data.frame(
  ColPrefix = shared_prefixes,
  stringsAsFactors = FALSE
)
pair_map$SampleID <- ifelse(
  pair_map$ColPrefix %in% names(SAMPLE_ALIASES),
  SAMPLE_ALIASES[pair_map$ColPrefix],
  pair_map$ColPrefix
)
pair_map$Identifier <- id_lookup[pair_map$SampleID]

unmapped <- pair_map$ColPrefix[is.na(pair_map$Identifier)]
if (length(unmapped) > 0) {
  warning("Samples not found in ", IDS_FILE, ", skipped: ",
          paste(unmapped, collapse = ", "))
  pair_map <- pair_map[!is.na(pair_map$Identifier), ]
}

if (nrow(pair_map) == 0) {
  stop("No sample pairs could be matched to the metadata sheet. ",
       "Check UNTREATED_SUFFIX / TREATED_SUFFIX and the SampleID column.")
}

cat("\nPaired samples:\n")
print(pair_map)


# =============================================================================
# SECTION 4: FOLD CHANGE
# =============================================================================

fold_change_data <- data.frame(row.names = rownames(norm_data))

for (i in seq_len(nrow(pair_map))) {
  # Lane columns use the code as written in the data, not the aliased one.
  sample_code <- pair_map$ColPrefix[i]
  final_id    <- pair_map$Identifier[i]

  untreated_col <- paste0(sample_code, UNTREATED_SUFFIX)
  treated_col   <- paste0(sample_code, TREATED_SUFFIX)

  # Column is named for the publication identifier, not the internal code.
  fold_change_data[[final_id]] <-
    as.numeric(norm_data[[treated_col]]) / as.numeric(norm_data[[untreated_col]])
}

logfc_data <- log2(fold_change_data)

# A zero denominator produces Inf. Treat those as "no measurable change"
# rather than dropping the strain: set them to NA, then to 0.
n_inf <- sum(sapply(logfc_data, is.infinite))
if (n_inf > 0) cat("\nSet", n_inf, "infinite value(s) to 0.\n")

is.na(logfc_data) <- sapply(logfc_data, is.infinite)
logfc_data[is.na(logfc_data)] <- 0


# =============================================================================
# SECTION 5: WRITE
# =============================================================================

dir.create(dirname(OUT_FILE), showWarnings = FALSE, recursive = TRUE)
write.csv(logfc_data, OUT_FILE, row.names = TRUE)

cat("\n=============================================================\n")
cat("log2 fold change complete.\n")
cat("  genes   : ", nrow(logfc_data), "\n", sep = "")
cat("  strains : ", ncol(logfc_data), "\n", sep = "")
cat("  columns : ", paste(colnames(logfc_data), collapse = ", "), "\n", sep = "")
cat("  output  : ", OUT_FILE, "\n", sep = "")
cat("Next step: R/03_heatmap_spr.R\n")
cat("=============================================================\n")
