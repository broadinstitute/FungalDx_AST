# =============================================================================
# 03_heatmap_spr.R
#
# STEP 3 of the Fungal AST pipeline - the publication figure.
#
# Renders a single panel containing, top to bottom:
#   - a susceptibility strip (greyscale: black = susceptible, grey = resistant)
#   - a log2 fold-change heatmap (genes x strains, blue-white-red)
#   - an SPR dot plot, one point per strain
#   - the MIC for each strain, as text
#
# Columns are ordered by MIC, most susceptible on the left.
#
# WHAT SPR IS
#   SPR (Signature Projection Ratio) collapses a strain's whole transcriptional
#   response into one number: how strongly it points along the "susceptible
#   response" direction. The axial vector is the mean log2FC profile of the
#   susceptible strains that went through RNAseq. Each strain's response is
#   projected onto that axis and normalized. A susceptible strain mounts the
#   canonical response and scores near 1; a resistant strain barely responds to
#   the drug and scores near 0. 
#
# INPUT   data/logfc/logfc_albicansFluc_compiled.csv   from 02_compute_logfc.R
#         metadata/albicansFluc_IDs.csv
# OUTPUT  figures/albicansFluc_heatmap_SPR.svg  (and .pdf)
#
# This is the single-run version: one measurement per strain, plain dots, no
# error bars. Averaging replicate runs is a separate analysis.
# =============================================================================


# =============================================================================
# SECTION 0: LIBRARIES
# =============================================================================

if (!requireNamespace("BiocManager",    quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) BiocManager::install("ComplexHeatmap")
if (!requireNamespace("circlize",       quietly = TRUE)) BiocManager::install("circlize")
if (!requireNamespace("dplyr",          quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("svglite",        quietly = TRUE)) install.packages("svglite")

library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(grid)
library(svglite)


# =============================================================================
# SECTION 1: USER CONFIG
# =============================================================================

# ---- Paths ---------------------------------------------------------------
LOGFC_FILE    <- "data/logfc/logfc_albicansFluc_compiled.csv"
METADATA_FILE <- "metadata/albicansFluc_IDs.csv"
OUT_DIR       <- "figures"
OUT_BASENAME  <- "albicansFluc_heatmap_SPR"

# ---- Gene name cleanup ---------------------------------------------------
# Probe names carry a panel prefix and a class tag, e.g.
#   CaFluc3_R_ERG11          -> ERG11
#   CaFluc1_Rup_CDR1         -> CDR1
#   CaFluc1_R_CAALFM_C301540WA -> CAALFM_C301540WA
# Set FALSE to keep full probe names on the row labels.
STRIP_PROBE_PREFIX <- TRUE

# ---- Gene subset ---------------------------------------------------------
# NULL   = use every gene in the matrix, with rows clustered.
# vector = use exactly these genes, in exactly this order, clustering off.
#
# IMPORTANT: the subset is applied BEFORE SPR is computed, so SPR is calculated
# over the selected genes only. That is deliberate and matches the published
# figure - a curated signature should be scored on its own genes. If you want
# SPR over the full panel, leave this NULL.
#
# Example curated set (the ten genes separating S from R most strongly):
#   GENES_TO_KEEP <- c("ERG4", "ERG11", "PMC1", "ERG13", "ERG2",
#                      "ERG5", "ERG25", "ERG3", "ERG6", "CAALFM_C301540WA")
GENES_TO_KEEP <- NULL

CLUSTER_ROWS   <- TRUE     # ignored when GENES_TO_KEEP is set
SHOW_ROW_NAMES <- TRUE

# ---- Column ordering -----------------------------------------------------
# "MIC_asc"  = lowest MIC (most susceptible) on the LEFT
# "MIC_desc" = highest MIC (most resistant) on the LEFT
# Ties within an MIC are broken by descending SPR.
COL_ORDER <- "MIC_asc"

# ---- Susceptibility colours ---------------------------------------------
# Keys must match the values in the metadata Susceptibility column.
SUS_COLORS <- c(
  "S"       = "black",
  "SDD"     = "black",
  "I"       = "grey30",
  "R"       = "grey60",
  "Unknown" = "grey30"
)

# Categories drawn as black in the SPR panel; everything known is grey60.
SUSCEPTIBLE_LEVELS <- c("S", "SDD")

# ---- Figure ---------------------------------------------------------------
FIG_WIDTH  <- 10          # inches; increase for many strains
FIG_HEIGHT <- 8
FIG_TITLE  <- "C. albicans Fluconazole log2 Fold Change - Heatmap + SPR"

POINT_SIZE_MM <- 4
SPR_PANEL_CM  <- 3.8


# =============================================================================
# SECTION 2: LOAD DATA
# =============================================================================

logfc_mat <- as.matrix(read.csv(LOGFC_FILE, row.names = 1, check.names = FALSE))
metadata  <- read.csv(METADATA_FILE, stringsAsFactors = FALSE)

stopifnot(all(c("Identifier", "MIC", "Susceptibility", "RNAseq") %in%
                colnames(metadata)))

cat("logFC matrix: ", nrow(logfc_mat), " genes x ", ncol(logfc_mat),
    " strains\n", sep = "")
cat("Metadata rows: ", nrow(metadata), "\n", sep = "")


# ---- Strip probe prefixes to bare gene names -----------------------------
# Matches "<panel>_<class>_" at the start of the name, where class is one of
# the probe-type tags used across the lab's panels. Anything that does not
# match is left untouched, so running this twice is harmless.
if (STRIP_PROBE_PREFIX) {
  clean <- sub("^[^_]+_(Rup|Rdn|Rb|Ri|R|B|C)_", "", rownames(logfc_mat))

  dups <- clean[duplicated(clean)]
  if (length(dups) > 0) {
    warning("Probe names collide after stripping prefixes: ",
            paste(unique(dups), collapse = ", "),
            ". Suffixes appended to keep them distinct.")
    clean <- make.unique(clean)
  }

  n_changed <- sum(clean != rownames(logfc_mat))
  cat("Stripped probe prefix from ", n_changed, " of ", nrow(logfc_mat),
      " row names.\n", sep = "")
  rownames(logfc_mat) <- clean
}


# =============================================================================
# SECTION 3: GENE SUBSET
# =============================================================================
# Applied before SPR - see the note in SECTION 1.

if (!is.null(GENES_TO_KEEP)) {
  matched <- match(GENES_TO_KEEP, rownames(logfc_mat))

  missing <- GENES_TO_KEEP[is.na(matched)]
  if (length(missing) > 0) {
    warning("Genes not found in the matrix, dropped: ",
            paste(missing, collapse = ", "))
  }

  logfc_mat    <- logfc_mat[matched[!is.na(matched)], , drop = FALSE]
  cluster_rows <- FALSE
} else {
  cluster_rows <- CLUSTER_ROWS
}

cat("Genes in figure: ", nrow(logfc_mat), "\n", sep = "")
if (nrow(logfc_mat) == 0) stop("No genes left after filtering.")


# =============================================================================
# SECTION 4: MATCH COLUMNS TO METADATA
# =============================================================================

col_meta <- data.frame(Identifier = colnames(logfc_mat), stringsAsFactors = FALSE)
col_meta <- left_join(col_meta, metadata, by = "Identifier")

missing_meta <- col_meta$Identifier[is.na(col_meta$Susceptibility)]
if (length(missing_meta) > 0) {
  warning("No metadata for: ", paste(missing_meta, collapse = ", "),
          ". Shown as Unknown.")
}

col_meta$Susceptibility[is.na(col_meta$Susceptibility)] <- "Unknown"
col_meta$MIC[is.na(col_meta$MIC)] <- NA_real_


# ---- Identify the RNAseq derivation strains ------------------------------
# The metadata RNAseq column marks which strains were used to derive the
# signature, in the form "Y S", "Y R", "Y S1", "Y R2", or "N".
#   Y = went through RNAseq
#   S = susceptible derivation strain, R = resistant derivation strain
rnaseq_strains <- metadata %>% filter(grepl("Y", RNAseq))

DERIV_S_IDS <- rnaseq_strains %>% filter(grepl("S", RNAseq)) %>% pull(Identifier)
DERIV_R_IDS <- rnaseq_strains %>% filter(grepl("R", RNAseq)) %>% pull(Identifier)

cat("DerivS identifiers: ", paste(DERIV_S_IDS, collapse = ", "), "\n", sep = "")
cat("DerivR identifiers: ", paste(DERIV_R_IDS, collapse = ", "), "\n", sep = "")

if (length(DERIV_S_IDS) == 0) {
  stop("No derivS strains found - the RNAseq column needs entries like 'Y S'.")
}
if (length(DERIV_R_IDS) == 0) {
  stop("No derivR strains found - the RNAseq column needs entries like 'Y R'.")
}


# =============================================================================
# SECTION 5: SPR CALCULATION
# =============================================================================

#' Signature Projection Ratio.
#'
#' Strains are split into three groups: the susceptible derivation strains,
#' the resistant derivation strains, and everything else ("valid" - the
#' validation set that did not contribute to the signature).
#'
#' The axial vector is the mean log2FC profile across the susceptible
#' derivation strains: the canonical drug response. Each strain's own profile
#' is projected onto it by an elementwise product and column sum, giving a
#' scalar `s` per strain. The score is then
#'
#'     SPR = s * |s| / (sum(axial^2))^2
#'
#' The `s * |s|` keeps the sign while squaring the magnitude, so strains that
#' respond in the opposite direction go negative. The denominator normalizes
#' by the axial vector's own length so a derivation-susceptible strain lands
#' near 1.
#'
#' @param logfc_data Matrix of log2 fold changes, genes x strains.
#' @param derivS Identifiers of the susceptible derivation strains.
#' @param derivR Identifiers of the resistant derivation strains.
#' @return A list with spr_S, spr_R, spr_valid, and spr_all (all three
#'   concatenated, named by strain identifier).
calculate_SPR <- function(logfc_data, derivS, derivR) {

  # ---- Split into the three strain groups --------------------------------
  allRespNormLog_derivS <- logfc_data[,  colnames(logfc_data) %in% derivS, drop = FALSE]
  allRespNormLog_derivR <- logfc_data[,  colnames(logfc_data) %in% derivR, drop = FALSE]
  allRespNormLog_valid  <- logfc_data[, !colnames(logfc_data) %in% c(derivS, derivR),
                                      drop = FALSE]

  # ---- Axial vector: the canonical susceptible response ------------------
  axialVec       <- rowMeans(allRespNormLog_derivS, na.rm = TRUE)
  sumSqrAxialVec <- sum(axialVec^2, na.rm = TRUE)

  # ---- Project each strain onto the axis ---------------------------------
  # sweep multiplies every row by that gene's axial weight; colSums then
  # completes the dot product, one scalar per strain.
  strainDisp_derivS <- sweep(allRespNormLog_derivS, 1, axialVec, FUN = "*")
  strainDisp_derivR <- sweep(allRespNormLog_derivR, 1, axialVec, FUN = "*")
  strainDisp_valid  <- sweep(allRespNormLog_valid,  1, axialVec, FUN = "*")

  sumStrainDisp_derivS <- colSums(strainDisp_derivS, na.rm = TRUE)
  sumStrainDisp_derivR <- colSums(strainDisp_derivR, na.rm = TRUE)
  sumStrainDisp_valid  <- colSums(strainDisp_valid,  na.rm = TRUE)

  # ---- Signed-square, normalized by the axial vector length --------------
  spr_S     <- sumStrainDisp_derivS * abs(sumStrainDisp_derivS) / sumSqrAxialVec^2
  spr_R     <- sumStrainDisp_derivR * abs(sumStrainDisp_derivR) / sumSqrAxialVec^2
  spr_valid <- sumStrainDisp_valid  * abs(sumStrainDisp_valid)  / sumSqrAxialVec^2

  spr_all <- c(spr_S, spr_R, spr_valid)

  return(list(
    spr_S     = spr_S,
    spr_R     = spr_R,
    spr_valid = spr_valid,
    spr_all   = spr_all
  ))
}

# Computed on the unsorted matrix so the names line up naturally.
spr_results <- calculate_SPR(logfc_mat, derivS = DERIV_S_IDS, derivR = DERIV_R_IDS)

cat("\nSPR values (pre-sort):\n")
print(round(spr_results$spr_all, 4))


# =============================================================================
# SECTION 6: ORDER COLUMNS BY MIC, THEN BY SPR
# =============================================================================

sort_df <- data.frame(
  Identifier = col_meta$Identifier,
  MIC        = col_meta$MIC,
  SPR        = spr_results$spr_all[col_meta$Identifier],
  stringsAsFactors = FALSE
)

col_order <- if (COL_ORDER == "MIC_asc") {
  order(sort_df$MIC, -sort_df$SPR, na.last = TRUE)
} else if (COL_ORDER == "MIC_desc") {
  order(-sort_df$MIC, -sort_df$SPR, na.last = TRUE)
} else {
  stop('COL_ORDER must be "MIC_asc" or "MIC_desc", got: ', COL_ORDER)
}

logfc_mat <- logfc_mat[, col_order, drop = FALSE]
col_meta  <- col_meta[col_order, ]

spr_values <- spr_results$spr_all[colnames(logfc_mat)]

cat("\nColumn order (left to right): ",
    paste(col_meta$Identifier, collapse = ", "), "\n", sep = "")
cat("\nSPR values (sorted):\n")
print(round(spr_values, 4))


# =============================================================================
# SECTION 7: BUILD THE FIGURE
# =============================================================================

# ---- Colour scale, symmetric around zero ---------------------------------
# Symmetry matters: it keeps white at exactly 0, so induced and repressed
# genes of equal magnitude read as equally saturated.
max_abs <- max(abs(logfc_mat), na.rm = TRUE)
col_fun <- colorRamp2(c(-max_abs, 0, max_abs), c("blue", "white", "red"))

# ---- Greyscale point fills, matching the susceptibility strip ------------
point_fills <- ifelse(col_meta$Susceptibility %in% SUSCEPTIBLE_LEVELS,
                      "black", "grey60")
point_fills[is.na(col_meta$Susceptibility) |
              col_meta$Susceptibility == "Unknown"] <- "grey30"

# ---- Top annotation: susceptibility strip --------------------------------
ha_top <- HeatmapAnnotation(
  Susceptibility       = col_meta$Susceptibility,
  col                  = list(Susceptibility = SUS_COLORS),
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = 9),
  show_legend          = TRUE
)

# ---- Bottom annotation: SPR dot plot, then the MIC row -------------------
ha_bottom <- HeatmapAnnotation(
  SPR = anno_points(
    spr_values,
    gp         = gpar(col = point_fills, fill = point_fills),
    pch        = 16,
    size       = unit(POINT_SIZE_MM, "mm"),
    axis_param = list(side = "left", gp = gpar(fontsize = 7)),
    height     = unit(SPR_PANEL_CM, "cm")
  ),
  MIC = anno_text(
    as.character(col_meta$MIC),
    gp       = gpar(fontsize = 8, border = "black"),
    height   = unit(0.9, "cm"),
    just     = "centre",
    location = 0.5,
    rot      = 0
  ),
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = 9, fontface = "bold")
)

# ---- Main heatmap --------------------------------------------------------
ht <- Heatmap(
  logfc_mat,
  name              = "log2FC",
  col               = col_fun,
  cluster_rows      = cluster_rows,
  cluster_columns   = FALSE,      # column order is set by MIC, never clustered
  show_row_dend     = FALSE,
  top_annotation    = ha_top,
  bottom_annotation = ha_bottom,
  column_labels     = col_meta$Identifier,
  column_names_rot  = 45,
  column_names_gp   = gpar(fontsize = 8),
  column_names_side = "bottom",
  row_names_gp      = gpar(fontsize = 6),
  show_row_names    = SHOW_ROW_NAMES,
  column_title      = FIG_TITLE,
  column_title_gp   = gpar(fontsize = 13, fontface = "bold"),
  heatmap_legend_param = list(
    title     = "log2FC",
    title_gp  = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)

render_figure <- function() {
  draw(ht, padding = unit(c(8, 8, 8, 8), "mm"), heatmap_legend_side = "right")
}


# =============================================================================
# SECTION 8: SAVE
# =============================================================================

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

out_svg <- file.path(OUT_DIR, paste0(OUT_BASENAME, ".svg"))
out_pdf <- file.path(OUT_DIR, paste0(OUT_BASENAME, ".pdf"))
out_spr <- file.path(OUT_DIR, paste0(OUT_BASENAME, "_SPR_values.csv"))

svglite(out_svg, width = FIG_WIDTH, height = FIG_HEIGHT)
render_figure()
dev.off()

pdf(out_pdf, width = FIG_WIDTH, height = FIG_HEIGHT)
render_figure()
dev.off()

# The underlying numbers, for supplementary tables
write.csv(
  data.frame(
    Identifier     = col_meta$Identifier,
    MIC            = col_meta$MIC,
    Susceptibility = col_meta$Susceptibility,
    SPR            = as.numeric(spr_values),
    stringsAsFactors = FALSE
  ),
  out_spr, row.names = FALSE
)

# On-screen copy
render_figure()

cat("\n=============================================================\n")
cat("Figure saved to:\n  ", out_svg, "\n  ", out_pdf, "\n", sep = "")
cat("SPR table saved to:\n  ", out_spr, "\n", sep = "")
cat("=============================================================\n")
