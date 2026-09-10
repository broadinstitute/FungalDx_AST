# =============================================================================
# 04_gene_pathway_overlap.R
#
# STEP 4 of the Fungal AST pipeline - cross-species/drug-class gene recurrence.
#
# Answers: for a given drug class (e.g. azoles), which genes turn up as
# annotated hits in more than one species/strain-panel dataset, and how many
# times? Example: ERG11 appearing in the logFC tables for C. albicans/Fluc,
# C. glabrata/Fluc, C. parapsilosis/Fluc, AND C. auris/Vori scores 4/4.
#
# Probes with no annotated common gene name (bare systematic locus tags such
# as CAALFM_C301540WA, GVI51_G04301, CJI97_000267T0, or a bare locus number
# like the C. parapsilosis Cf_Rup_101870 probes) are dropped rather than
# guessed at, since they cannot be compared across species without an
# orthology call (e.g. OrthoFinder) - see the note in SECTION 3.
#
# INPUT   data/logfc/logfc_albicansFluc_compiled.csv
#         data/logfc/logfc_glabrataFluc_formatted.csv
#         data/logfc/logfc_parapsilosisFluc_formatted.csv
#         data/logfc/logfc_aurisVori_formatted.csv
# OUTPUT  figures/azole_gene_overlap.svg  (and .pdf)
#         figures/azole_gene_overlap_table.csv
# =============================================================================


# =============================================================================
# SECTION 0: LIBRARIES
# =============================================================================

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr",   quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("svglite", quietly = TRUE)) install.packages("svglite")

library(ggplot2)
library(dplyr)
library(svglite)


# =============================================================================
# SECTION 1: USER CONFIG
# =============================================================================

# ---- Datasets in this drug class -----------------------------------------
# One entry per species/strain-panel logFC table (rows = probes, first column
# = probe name). All four below are azole datasets (fluconazole or
# voriconazole); to build the same figure for a different drug class, point
# this at that class's logFC tables instead.
DRUG_CLASS_LABEL <- "Azoles (Fluconazole / Voriconazole)"

# Paths come from the pair registry, so this stays in step with 01-03 and
# with whatever 02 last wrote. Run 02 for each pair listed here first.
source("R/00_species_config.R")

DATASETS <- vapply(
  c("C. albicans / Fluc"     = "albicansFluc",
    "C. glabrata / Fluc"     = "glabrataFluc",
    "C. parapsilosis / Fluc" = "parapsilosisFluc",
    "C. auris / Vori"        = "aurisVori"),
  function(k) get_species_config(k)$logfc_file,
  character(1)
)

# A pair whose logFC table has not been generated yet is skipped rather than
# stopping the figure - e.g. glabrataFluc, which has a strain sheet but no
# raw exports in the repo yet.
missing_sets <- DATASETS[!file.exists(DATASETS)]
if (length(missing_sets) > 0) {
  warning("Skipping ", length(missing_sets), " dataset(s) with no logFC table ",
          "yet - run 02 for them first: ",
          paste(names(missing_sets), collapse = ", "))
  DATASETS <- DATASETS[file.exists(DATASETS)]
}
if (length(DATASETS) < 2) {
  stop("Need at least two datasets to look at overlap; ", length(DATASETS),
       " available.")
}

# ---- Threshold -------------------------------------------------------------
# Only genes recurring in at least this many datasets are plotted.
MIN_OCCURRENCES <- 2

# ---- Rows that are metadata, not probes -----------------------------------
# Some of these tables carry MIC/susceptibility rows above the probe rows.
NON_GENE_ROWS <- c("MIC", "SIR")

# ---- Output -----------------------------------------------------------------
OUT_DIR      <- "figures"
OUT_BASENAME <- "azole_gene_overlap"
FIG_WIDTH    <- 8
FIG_HEIGHT   <- 5.5


# =============================================================================
# SECTION 2: HELPERS
# =============================================================================

# Same probe-naming convention used in 03_heatmap_spr.R:
#   <panel>_<class>_<gene>   e.g. CaFluc3_R_ERG11 -> ERG11
# class is one of Rup / Rdn / Rb / Ri / R / B / C. Names that don't match this
# shape (e.g. a bare "ERG11" already, or a bare locus tag with no panel
# prefix) are left untouched.
strip_probe_prefix <- function(x) {
  sub("^[^_]+_(Rup|Rdn|Rb|Ri|R|B|C)_", "", x)
}

# A name counts as "annotated" only if it looks like a common gene symbol.
# Anything else is a systematic locus tag with no assigned name, and gets
# dropped - comparing those across species requires an orthology call
# (OrthoFinder), which is future work, not a name match.
#   - still has an underscore after stripping (e.g. CAALFM_C301540WA,
#     GVI51_G04301, CJI97_000267T0) -> locus tag, drop
#   - pure digits (e.g. the C. parapsilosis "107090" locus numbers) -> drop
#   - yeast-style systematic ORF name, e.g. YNL208W -> drop
is_annotated_gene <- function(x) {
  !grepl("_", x) &
    !grepl("^[0-9]+$", x) &
    !grepl("^Y[A-P][LR][0-9]{3}[WC](-[AB])?$", x)
}


# =============================================================================
# SECTION 3: LOAD + CLEAN EACH DATASET
# =============================================================================

species_genes <- list()

for (label in names(DATASETS)) {
  path <- DATASETS[[label]]
  raw  <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)

  probe_names <- trimws(raw[[1]])
  probe_names <- probe_names[!is.na(probe_names) &
                                probe_names != "" &
                                !(probe_names %in% NON_GENE_ROWS)]

  cleaned <- strip_probe_prefix(probe_names)

  kept    <- unique(cleaned[is_annotated_gene(cleaned)])
  dropped <- unique(cleaned[!is_annotated_gene(cleaned)])

  species_genes[[label]] <- kept

  cat(sprintf("%-24s kept %2d annotated genes, dropped %2d unannotated locus tag(s)\n",
              label, length(kept), length(dropped)))
  if (length(dropped) > 0) {
    cat("    dropped: ", paste(dropped, collapse = ", "), "\n", sep = "")
  }
}


# =============================================================================
# SECTION 4: COUNT OCCURRENCES ACROSS DATASETS
# =============================================================================

all_genes <- unique(unlist(species_genes))

occurrence <- data.frame(
  gene  = all_genes,
  count = sapply(all_genes, function(g) sum(sapply(species_genes, function(s) g %in% s))),
  datasets = sapply(all_genes, function(g) {
    paste(names(species_genes)[sapply(species_genes, function(s) g %in% s)], collapse = "; ")
  }),
  stringsAsFactors = FALSE
)

recurrent <- occurrence %>%
  filter(count >= MIN_OCCURRENCES) %>%
  arrange(desc(count), gene)

cat("\nGenes recurring in >=", MIN_OCCURRENCES, "datasets:", nrow(recurrent), "\n")
print(recurrent)

if (nrow(recurrent) == 0) {
  stop("No genes met the MIN_OCCURRENCES threshold - nothing to plot.")
}


# =============================================================================
# SECTION 5: PLOT
# =============================================================================

recurrent$gene <- factor(recurrent$gene, levels = recurrent$gene)

p <- ggplot(recurrent, aes(x = gene, y = count)) +
  geom_col(fill = "grey30", width = 0.65) +
  scale_y_continuous(breaks = seq(1, length(DATASETS)), limits = c(0, length(DATASETS))) +
  labs(
    title    = paste0(DRUG_CLASS_LABEL, " - genes recurring across datasets"),
    subtitle = paste0("Genes with an annotated name appearing in ≥", MIN_OCCURRENCES,
                       " of ", length(DATASETS), " datasets"),
    x = NULL,
    y = "Number of datasets"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

print(p)


# =============================================================================
# SECTION 6: WRITE
# =============================================================================

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

svg_path <- file.path(OUT_DIR, paste0(OUT_BASENAME, ".svg"))
pdf_path <- file.path(OUT_DIR, paste0(OUT_BASENAME, ".pdf"))
csv_path <- file.path(OUT_DIR, paste0(OUT_BASENAME, "_table.csv"))

ggsave(svg_path, p, width = FIG_WIDTH, height = FIG_HEIGHT, device = svglite::svglite)
ggsave(pdf_path, p, width = FIG_WIDTH, height = FIG_HEIGHT)
write.csv(recurrent, csv_path, row.names = FALSE)

cat("\n=============================================================\n")
cat("Gene overlap figure complete.\n")
cat("  genes plotted : ", nrow(recurrent), "\n", sep = "")
cat("  figure        : ", svg_path, " / ", pdf_path, "\n", sep = "")
cat("  table         : ", csv_path, "\n", sep = "")
cat("=============================================================\n")
