#!/usr/bin/env Rscript

suppressMessages(library(docopt))
suppressMessages(library(DESeq2))

#source("/gsap/fungaldx/Cglabrata/Brooks_work/R_scripts/Control_Gene_Scripts/utils.R")
source("/gsap/fungaldx/Cglabrata/Brooks_work/R_scripts/Control_Gene_Scripts/utils_Updated.R")
#source("/broad/IDP-Dx_work/nirmalya/gene_selection/ControlGenes_lib.R")
source ("/gsap/fungaldx/Cglabrata/Brooks_work/R_scripts/Control_Gene_Scripts/ControlGenes_lib_non_hardcoded.R")

'Control genes script over combined dataset

Usage: ControlGenes_auto.R (--from_all|--from_part) --stag <stag> -c <counts> -o <outdir> [--p_cand <padj_cand> --l_cand <l2fc_cand> --p_phase <padj_phase> --l_phase <l2fc_phase> --base_lim <baseMean % limit> --no_res] 

options:
  -c <count> --count <count>
  -o <outdir> --outdir <outdir>
  --stag <stag>
  --p_cand <padj_cand> [default: 0.05]
  --l_cand <l2fc_cand> [default: 0.5]
  --p_phase <padj_phase> [default: 0.05]
  --l_phase <l2fc_phase> [default: 0.5]
  --base_lim <baseMean % limit> [default: 50]' -> doc
# what are the options? Note that stripped versions of the parameters are added to the returned list

opts <- docopt(doc)
str(opts)  


count_file <- opts$count
outdir <- opts$outdir
padj1 <- as.numeric(opts$p_cand)
lfcth1 <- as.numeric(opts$l_cand)
padj2 <- as.numeric(opts$p_phase)
lfcth2 <- as.numeric(opts$l_phase)
from_all <- opts$from_all
from_part <- opts$from_part
stag <- opts$stag
base_lim <- as.numeric(opts$base_lim)
no_res <- opts$no_res

data_usage_str <- NULL

if (from_part) {
    data_usage_str <- "partdata"
} else if (from_all) {
    data_usage_str <- "alldata"
}

dir.create(outdir, recursive = TRUE)
logfile <- paste0(outdir, "/", "logfile_", data_usage_str, ".txt")
print(logfile)
file.create(logfile)

print_log(paste0("count_file: ", count_file))
print_log(paste0("p_cand: ", padj1))
print_log(paste0("l_cand: ", lfcth1))
print_log(paste0("p_phase: ", padj2))
print_log(paste0("l_phase: ", lfcth2))
print_log(paste0("no_res: ", no_res))

print_log(paste0("count_file: " , count_file))
count_tbl <- get_count_tbl(count_file)
print_log(paste0("count_tbl: ", dim(count_tbl)[2]))

sample_ids <- colnames(count_tbl)
sample_groups <- get_sample_groups(sample_ids)
print_log("sample_groups")
print(sample_groups)
#lapply(sample_groups, write, logfile, append=TRUE, ncolumns=1000)
#print_log(sample_groups)
colData <- get_colData(sample_groups)
print(str(colData))
lsamples <- rownames(colData)
countData <- count_tbl[, lsamples] 
print("countData")
print(str(countData))
has_res <- !no_res
deseq_main(countData = round(countData), outdir, stag, sample_groups, padj1, lfcth1, padj2, lfcth2, data_usage_str, base_lim = base_lim, hasRes = has_res)

