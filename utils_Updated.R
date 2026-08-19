get_splitted_parts <- function(lstr, reg_str) {
    # lstr: a string like "S,R" or "4h,8h,16h"
    # reg_str: regex pattern for splitting (e.g. ",")
    part_lst <- strsplit(lstr, reg_str)
    parts <- trimws(part_lst[[1]])
    return(parts)
}

combine_parts <- function(sample_ids, strains, treatments, times) {
    # Ensure times is a character vector of timepoint labels (e.g. "4h", "8h", "16h")
    # This returns all sample_ids that match any combination of strain, treatment, and timepoint
    combinations <- expand.grid(strains, treatments, times, stringsAsFactors = FALSE)
    combined_res_pre <- apply(combinations, 1, paste, collapse = "_")
    combined_res_lst <- lapply(combined_res_pre, function(pre, sample_ids) {
        matched_ones <- sample_ids[startsWith(sample_ids, pre)]
        return(matched_ones)
    }, sample_ids)
    combined_res <- unlist(combined_res_lst, use.names = FALSE)
    return(combined_res)
}

get_tp_parts <- function(sample_ids, time_pos) {
    vals <- lapply(sample_ids, function(x) {strsplit(x, '_')[[1]][time_pos]})
    tp_parts <- sort(unique(unlist(vals)))
    return (tp_parts)
}

get_samples <- function(sample_ids, pos, val) {
    if (val == "S" || val == "R") {
        #print(paste0("THIS IS THE VAL:", val))
      retvals <- lapply(sample_ids, function(x, lpos, lval){
        lval1 <- strsplit(x, '_')[[1]][lpos]
        #print(paste0("lval1:",lval1))
        #print(paste0("lval",lval))
        #print(paste0("This is lval1[1]:", substr(lval1,1,1)))
        if(substr(lval1,1,1) == lval) {
        

          return (TRUE)}
        else
          return (FALSE)
      }, pos, val)
    }
    else {
      retvals <- lapply(sample_ids, function(x, lpos, lval){
        lval1 <- strsplit(x, '_')[[1]][lpos]
        if(lval1 == lval)
          return (TRUE)
        else
          return (FALSE)
      }, pos, val)
      
    }
    
    retvals1 <- unlist(retvals)
    ret_samples <- sample_ids[retvals1]
    return (ret_samples)
}


# PATCHED FUNCTION 1: get_sample_groups()
get_sample_groups <- function(sample_ids) {
    straintype_pos <- 2
    treatment_pos <- 4
    time_pos <- 5

    sus_samples <- get_samples(sample_ids, straintype_pos, 'S')
    res_samples <- get_samples(sample_ids, straintype_pos, 'R')
    treated_samples <- get_samples(sample_ids, treatment_pos, '+')
    untreated_samples <- get_samples(sample_ids, treatment_pos, '-')
    tp_parts <- get_tp_parts(sample_ids, time_pos)

    retval <- list()
    for (j in seq_along(tp_parts)) {
        ltime_point <- tp_parts[j]
        ltime_samples <- get_samples(sample_ids, time_pos, ltime_point)

        treated_str <- paste0("sus_treated_time", ltime_point)
        retval[[treated_str]] <- Reduce(intersect, list(sus_samples, treated_samples, ltime_samples))

        untreated_str <- paste0("sus_untreated_time", ltime_point)
        retval[[untreated_str]] <- Reduce(intersect, list(sus_samples, untreated_samples, ltime_samples))

        res_treated_str <- paste0("res_treated_time", ltime_point)
        retval[[res_treated_str]] <- Reduce(intersect, list(res_samples, treated_samples, ltime_samples))

        res_untreated_str <- paste0("res_untreated_time", ltime_point)
        retval[[res_untreated_str]] <- Reduce(intersect, list(res_samples, untreated_samples, ltime_samples))
    }

    retval$exp_conds <- list("treated_parts" = c("+"), "untreated_parts" = c("-"), "tp_parts" = tp_parts)
    return(retval)
}


get_colData_small <- function(sample_groups, cond2, cond1) {
    lcond2_samples <- sample_groups[[cond2]]
    lcond1_samples <- sample_groups[[cond1]]
    lcond2_rep <- rep(cond2, length(lcond2_samples))
    lcond1_rep <- rep(cond1, length(lcond1_samples))
    condition <- c(lcond2_rep, lcond1_rep)
    all_samples <- c(lcond2_samples, lcond1_samples)

    colData <- data.frame(condition)
    colnames(colData) <- "condition"
    rownames(colData) <- all_samples
    return (colData)
}


# PATCHED FUNCTION 2: get_repeated_vals()
get_repeated_vals <- function(repeat_tag, sample_groups) {
    tp_parts <- sample_groups$exp_conds$tp_parts
    condition <- c()
    samples <- c()
    for (tp in tp_parts) {
        rep_str <- paste0(repeat_tag, tp)
        samples_arr <- unlist(sample_groups[[rep_str]])
        rep_arr <- rep(rep_str, length(samples_arr))
        condition <- c(condition, rep_arr)
        samples <- c(samples, samples_arr)
    }

    return(list("condition" = condition, "samples" = samples))
}


get_colData <- function(sample_groups, hasRes = TRUE) {

    sus_treated_vals <- get_repeated_vals("sus_treated_time", sample_groups)
    sus_untreated_vals <- get_repeated_vals("sus_untreated_time", sample_groups)
    condition <- c(sus_treated_vals[["condition"]], sus_untreated_vals[["condition"]])
    all_samples <- c(sus_treated_vals[["samples"]], sus_untreated_vals[["samples"]])

    if (hasRes) {
        res_treated_vals <- get_repeated_vals("res_treated_time", sample_groups)
        res_untreated_vals <- get_repeated_vals("res_untreated_time", sample_groups)
        condition <- c(condition, res_treated_vals[["condition"]], res_untreated_vals[["condition"]])
        all_samples <- c(all_samples, res_treated_vals[["samples"]], res_untreated_vals[["samples"]])
    }

    colData <- data.frame(condition)
    colnames(colData) <- "condition"
    rownames(colData) <- all_samples
    return (colData)
}


get_count_tbl <- function(count_file) {
    count_tbl1 <- read.csv(count_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
    genes <- rownames(count_tbl1)
    #cds_genes <- grepl("^CDS", genes)
    #count_tbl <- count_tbl1[cds_genes, ]
    return (count_tbl1)
}

get_baseMean_lim <- function(lres, base_lim) {
    lres_sorted <- lres[order(-lres$baseMean), ]
    rowlen <- dim(lres)[1]
    rowlen_lim <- floor(rowlen*(100-base_lim)/100)
    
    if (rowlen_lim > rowlen) {
        rowlen_lim = rowlen
    } else if (rowlen_lim < 1) {
        rowlen_lim = 1
    }
    baseMean_lim <- lres_sorted[rowlen_lim, "baseMean"]
    baseMean_100p <- lres_sorted[1, "baseMean"]
    print(paste0("baseMean_100p: ", baseMean_100p))
    return (baseMean_lim)
}


deseq_condwise <- function(dds, cond2, cond1, lfcth, padjth, altH = "greaterAbs", base_lim = 50) {
    lres <- results(dds, altHypothesis = altH, lfcThreshold = lfcth, contrast=c("condition", cond2, cond1))
    lres_top1 <- subset(lres, padj < padjth)
    baseMean_lim_val <- get_baseMean_lim(lres, base_lim)
    lres_top <- subset(lres_top1, baseMean >= baseMean_lim_val)
    retval <- list("lres_top" = lres_top, "lres" = lres)
    print_log(paste0("baseMean_", base_lim, "p: ", baseMean_lim_val))
    print_log(paste0("cond2: ", cond2))
    print_log(paste0("cond1: ", cond1))
    print_log(paste0("lfcth: ", lfcth))
    print_log(paste0("padjth: ", padjth))
    print_log(paste0("altH: ", altH))
    print_log(paste0("Gene_count: ", dim(lres)[1]))
    print_log(paste0("Gene_count after padj: ", dim(lres_top1)[1]))
    print_log(paste0("Gene_count_top: ", dim(lres_top)[1]))
    print_log("....................................")
    cat("\n")

    return (retval)
}


deseq_condwise_part <- function(countData, sample_groups, cond2, cond1, lfcth, padjth, altH = "greaterAbs", use_beta_prior = FALSE, base_lim = 50) {
    colData <- get_colData_small(sample_groups, cond2, cond1)
    lsamples <- rownames(colData)
    countData <- countData[, lsamples]
    print(colData)
    print_log(paste0("Info: beta prior is set to ", use_beta_prior))
    dds <- NULL
    if (use_beta_prior) {
        dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
    } else {
        dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ 0 + condition)
    }
    dds <- DESeq(dds, betaPrior = use_beta_prior)
    retval <- deseq_condwise(dds, cond2, cond1, lfcth, padjth, altH, base_lim)
    return (retval)
}

print_data_ma_between_conds <- function(out_tbl, outdir, stag, exp_conds, timepos, beta_prior_str, data_usage_str) {
    prefix <- paste(stag, exp_conds$"treated_parts"[1], '_VS_', exp_conds$"untreated_parts"[1], exp_conds$"tp_parts"[timepos], beta_prior_str, data_usage_str, sep = "_")
    print_log(paste0("prefix_str: ", prefix))
    outfile <- paste0(outdir, "/", prefix, ".txt")
    write.table(as.data.frame(out_tbl$"lres_top"), outfile, sep = "\t")
    prefix <- paste(stag, exp_conds$"treated_parts"[1], '_VS_', exp_conds$"untreated_parts"[1], exp_conds$"tp_parts"[timepos], beta_prior_str, data_usage_str, "allgenes", sep = "_")
    outfile <- paste0(outdir, "/", prefix, ".txt")
    write.table(as.data.frame(out_tbl$"lres"), outfile, sep = "\t")
    outfile_MA <- paste0(outdir, "/", prefix, "_MA.pdf")
    pdf(outfile_MA)
    plotMA(out_tbl$"lres", main=prefix, ylim=c(-5,5))
    dev.off()

}

# PATCHED FUNCTION 4: print_data_ma_between_tps()
print_data_ma_between_tps <- function(out_tbl, outdir, stag, exp_conds, cond, tp2, tp1, beta_prior_str, data_usage_str) {
    cond_str <- if (cond == "treated") exp_conds$"treated_parts"[1] else exp_conds$"untreated_parts"[1]
    prefix <- paste(stag, cond_str, tp2, "_VS_", tp1, beta_prior_str, data_usage_str, sep = "_")
    print_log(paste0("prefix_str: ", prefix))

    outfile <- paste0(outdir, "/", prefix, ".txt")
    write.table(as.data.frame(out_tbl$"lres_top"), outfile, sep = "\t")

    prefix_all <- paste(stag, cond_str, tp2, "_VS_", tp1, beta_prior_str, data_usage_str, "allgenes", sep = "_")
    outfile_all <- paste0(outdir, "/", prefix_all, ".txt")
    write.table(as.data.frame(out_tbl$"lres"), outfile_all, sep = "\t")

    outfile_MA <- paste0(outdir, "/", prefix_all, "_MA.pdf")
    pdf(outfile_MA)
    plotMA(out_tbl$"lres", main = prefix, ylim = c(-5, 5))
    dev.off()
}

get_z_score <- function(pval) {
    if (is.na(pval))
        return (0)
    val <- qnorm((pval/2), lower.tail = FALSE)
    return (val)
}

compile_z_scores <- function(deseq_res, time_set, pval_file, zscore_file) {

    timepoints <- length(time_set)
    phase_pval_lst <- list()
    lgenes_sorted <- NULL
    for (j in timepoints:2) {
        select_str <- paste0("sus_untreated_time", j, "_time1")
        select_str_obj <- (deseq_res[[select_str]])$"lres"
        select_z_str <- paste0(select_str, "_z")
        if (j == timepoints) {
            lgenes <- rownames(select_str_obj)
            lgenes_sorted <- sort(lgenes)
        }
        phase_pval_lst[[select_z_str]] <- select_str_obj[lgenes_sorted, "padj"]

    }
    pval_phase_tbl <- data.frame(phase_pval_lst)
    rownames(pval_phase_tbl) <- lgenes_sorted
    zscore_phase_tbl <- apply(pval_phase_tbl, c(1, 2), get_z_score)
    phase_invariant_z <- apply(zscore_phase_tbl, 1, min)

    time_pval_lst <- list()
    time_lfc_lst <- list()
    for (j in time_set) {
        select_str <- paste0("sus_treated_untreated_time", j)
        select_str_obj <- (deseq_res[[select_str]])$"lres"
        select_z_str <- paste0(select_str, "_z")
        time_pval_lst[[select_z_str]] <- select_str_obj[lgenes_sorted, "padj"]
        time_lfc_lst[[select_z_str]] <- select_str_obj[lgenes_sorted, "log2FoldChange"]

    }

pval_time_tbl <- data.frame(time_pval_lst)
    rownames(pval_time_tbl) <- lgenes_sorted
   
    lfc_time_tbl <- data.frame(time_lfc_lst)
    rownames(lfc_time_tbl) <- lgenes_sorted

    zscore_time_tbl_u <- apply(pval_time_tbl, c(1, 2), get_z_score)
    zscore_sign <- apply(lfc_time_tbl, c(1,2),
        function(x) {
            if (is.na(x)) { return (0) }
            if (x >=0) return (1) else return (-1)
        })
    zscore_time_tbl <- zscore_time_tbl_u * zscore_sign
    zscore_time_root_agv <- NULL
    if (length(time_set) == 1)
        zscore_time_root_agv <- zscore_time_tbl
    else
        zscore_time_root_agv <- data.frame(apply(zscore_time_tbl, 1, function(x) {sum(x)/sqrt(length(x))}))
    colnames(zscore_time_root_agv) <- "zscore_time_root_agv"

    zscore_tbl <- data.frame(zscore_phase_tbl, zscore_time_tbl, phase_invariant_z, zscore_time_root_agv)

    write.table(zscore_tbl, zscore_file,  sep = "\t")

    pval_tbl <- data.frame(pval_phase_tbl, pval_time_tbl)
    write.table(pval_tbl, pval_file, sep = "\t")

}

get_lst <- function(lst, end_str) {
    lnames <- names(lst)
    reg_str <- paste0(end_str, "$")
    lvals <- grepl(reg_str, lnames)
    ret_lst <- lst[lvals]
    return(ret_lst)
}


get_sample_groups_TB <- function(sample_ids, sus, treated, untreated, timepoins) {

    reg_str = "\\s*[,|;]\\s*"

    sus_parts <- get_splitted_parts(sus, reg_str)
    treated_parts <- get_splitted_parts(treated, reg_str)
    untreated_parts <- get_splitted_parts(untreated, reg_str)
    tp_parts <- get_splitted_parts(timepoints, reg_str)

    # Create prefixes for differenet categories

    retval <- list()

	for (j in seq_along(tp_parts)) {
        treated_str <- paste0("sus_treated_time", j)
        retval[[treated_str]] <- combine_parts(sample_ids, sus_parts, treated_parts, tp_parts[j])
        untreated_str <- paste0("sus_untreated_time", j)
        retval[[untreated_str]] <- combine_parts(sample_ids, sus_parts, untreated_parts, tp_parts[j])

    }

    exp_conds = list("treated_parts" = treated_parts, "untreated_parts" = untreated_parts,
            "tp_parts" = tp_parts)
    retval$"exp_conds" = exp_conds
    return (retval)
}

