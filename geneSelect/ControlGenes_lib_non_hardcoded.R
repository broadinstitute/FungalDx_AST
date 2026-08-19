#!/usr/bin/env Rscript
# Updated script to use flexible timepoints, always use actual timepoint labels (e.g. "4h") in key construction
suppressMessages(library(docopt))
suppressMessages(library(DESeq2))

deseq_from_part <- function(countData, sample_groups, padj1, lfcth1, padj2, lfcth2, base_lim = 50, hasRes = TRUE) {
    altH = "lessAbs"
    retval <- list()
    tp_parts <- sample_groups$exp_conds$tp_parts
    n_tp <- length(tp_parts)
    for (i in n_tp:2) {
        tp <- tp_parts[i]
        ret_str <- paste0("sus_untreated_time", tp, "_time", tp_parts[1])
        print_log(paste0("comp_str: ", ret_str))
        cond2 <- paste0("sus_untreated_time", tp)
        cond1 <- paste0("sus_untreated_time", tp_parts[1])
        retval[[ret_str]] <- deseq_condwise_part(countData, sample_groups, cond2, cond1, lfcth2, padj2, altH = altH, base_lim = base_lim)
    }
    if (hasRes) {
        for (i in n_tp:2) {
            tp <- tp_parts[i]
            ret_str <- paste0("res_untreated_time", tp, "_time", tp_parts[1])
            print_log(paste0("comp_str: ", ret_str))
            cond2 <- paste0("res_untreated_time", tp)
            cond1 <- paste0("res_untreated_time", tp_parts[1])
            retval[[ret_str]] <- deseq_condwise_part(countData, sample_groups, cond2, cond1, lfcth2, padj2, altH = altH, base_lim = base_lim)
        }
    }
    for (i in n_tp:2) {
        tp <- tp_parts[i]
        ret_str <- paste0("sus_treated_untreated_time", tp)
        print_log(paste0("comp_str: ", ret_str))
        cond2 <- paste0("sus_treated_time", tp)
        cond1 <- paste0("sus_untreated_time", tp)
        retval[[ret_str]] <- deseq_condwise_part(countData, sample_groups, cond2, cond1, lfcth1, padj1, altH = altH, base_lim = base_lim)
    }
    return(retval)
}

deseq_from_all <- function(countData, sample_groups, padj1, lfcth1, padj2, lfcth2, base_lim = 50, hasRes = TRUE) {
    colData <- get_colData(sample_groups)
    dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ 0 + condition)
    print_log("Using beta prior = FALSE")
    dds <- DESeq(dds, betaPrior = FALSE)
    altH = "lessAbs"
    retval <- list()
    tp_parts <- sample_groups$exp_conds$tp_parts
    n_tp <- length(tp_parts)
    for (i in n_tp:2) {
        tp <- tp_parts[i]
        ret_str <- paste0("sus_untreated_time", tp, "_time", tp_parts[1])
        print_log(paste0("comp_str: ", ret_str))
        cond2 <- paste0("sus_untreated_time", tp)
        cond1 <- paste0("sus_untreated_time", tp_parts[1])
        retval[[ret_str]] <- deseq_condwise(dds, cond2, cond1, lfcth2, padj2, altH = altH, base_lim = base_lim)
    }
    if (hasRes) {
        for (i in n_tp:2) {
            tp <- tp_parts[i]
            ret_str <- paste0("res_untreated_time", tp, "_time", tp_parts[1])
            print_log(paste0("comp_str: ", ret_str))
            cond2 <- paste0("res_untreated_time", tp)
            cond1 <- paste0("res_untreated_time", tp_parts[1])
            retval[[ret_str]] <- deseq_condwise(dds, cond2, cond1, lfcth2, padj2, altH = altH, base_lim = base_lim)
        }
    }
    for (i in n_tp:2) {
        tp <- tp_parts[i]
        ret_str <- paste0("sus_treated_untreated_time", tp)
        print_log(paste0("comp_str: ", ret_str))
        cond2 <- paste0("sus_treated_time", tp)
        cond1 <- paste0("sus_untreated_time", tp)
        retval[[ret_str]] <- deseq_condwise(dds, cond2, cond1, lfcth1, padj1, altH = altH, base_lim = base_lim)
    }
    return(retval)
}

print_data_MA <- function(deseq_res, outdir, stag, exp_conds, data_usage_str, hasRes = TRUE) {
    head_len <- 0
    beta_prior_str <- FALSE
    final_gene_lst <- list()
    tp_parts <- exp_conds$tp_parts

    stag_strain <- paste0(stag, "_S")
    for (i in seq(length(tp_parts), 2)) {
        tp <- tp_parts[i]
        lstr <- paste0("sus_untreated_time", tp, "_time", tp_parts[1])
        lstr2 <- paste0(lstr, "_genes")
        print_log(paste0("stag_strain: ", stag_strain))
        print_log(paste0("Printing data for: ", lstr))
        print(head(deseq_res[[lstr]]$"lres", head_len))
        print_data_ma_between_tps(deseq_res[[lstr]], outdir, stag_strain, exp_conds, "untreated", tp, tp_parts[1], beta_prior_str, data_usage_str)
        final_gene_lst[[lstr2]] <- rownames(deseq_res[[lstr]]$"lres_top")
    }

    if (hasRes) {
        stag_strain <- paste0(stag, "_R")
        for (i in seq(length(tp_parts), 2)) {
            tp <- tp_parts[i]
            lstr <- paste0("res_untreated_time", tp, "_time", tp_parts[1])
            lstr2 <- paste0(lstr, "_genes")
            print_log(paste0("stag_strain: ", stag_strain))
            print_log(paste0("Printing data for: ", lstr))
            print(head(deseq_res[[lstr]]$"lres", head_len))
            print_data_ma_between_tps(deseq_res[[lstr]], outdir, stag_strain, exp_conds, "untreated", tp, tp_parts[1], beta_prior_str, data_usage_str)
            final_gene_lst[[lstr2]] <- rownames(deseq_res[[lstr]]$"lres_top")
        }
    }

    stag_strain <- paste0(stag, "_S")
    for (i in seq(length(tp_parts), 2)) {
        tp <- tp_parts[i]
        lstr <- paste0("sus_treated_untreated_time", tp)
        lstr2 <- paste0(lstr, "_genes")
        print_log(paste0("stag_strain: ", stag_strain))
        print_log(paste0("Printing data for: ", lstr))
        print(head(deseq_res[[lstr]]$"lres", head_len))
        print_data_ma_between_conds(deseq_res[[lstr]], outdir, stag_strain, exp_conds, tp, beta_prior_str, data_usage_str)
        final_gene_lst[[lstr2]] <- rownames(deseq_res[[lstr]]$"lres_top")
    }
    return(final_gene_lst)
}

deseq_main <- function(countData, outdir, stag, sample_groups, padj1, lfcth1, padj2, lfcth2, data_usage_str, base_lim = 50, hasRes = TRUE) {
    deseq_res = NULL
    if (data_usage_str == "alldata") {
        print("Gene lists are prepared by fitting all data together.")
        deseq_res <- deseq_from_all(countData, sample_groups, padj1, lfcth1, padj2, lfcth2, base_lim = base_lim, hasRes)
    } else if (data_usage_str == "partdata") {
        print("Gene lists are prepared by fitting data separately to each condition.")
        deseq_res <- deseq_from_part(countData, sample_groups, padj1, lfcth1, padj2, lfcth2, base_lim, hasRes)
    }
    exp_conds <- sample_groups$"exp_conds"
    print("Printing exp_conds:")
    print(exp_conds)
    all_lst <- print_data_MA(deseq_res, outdir, stag, exp_conds, data_usage_str, hasRes)
    print_log("Intersecting all lists: ")
    print(str(all_lst))
    control_lst <- Reduce(intersect, all_lst)
    control_lst_count <- length(control_lst)
    print_log(paste0("Number of control genes: ", control_lst_count))
    write(control_lst, file = logfile, append = TRUE)
    # updated outfile to allow for dynamic timepoints
    tp_parts <- sample_groups$exp_conds$tp_parts
    last_tp <- tail(tp_parts, 1)
    last_key <- paste0("sus_treated_untreated_time", last_tp)
    print_log("Names in deseq_res:")
    print(names(deseq_res))
    print_log(paste("last_key:", last_key))
    print_log(paste("last_tp:", last_tp))
    if (!is.null(deseq_res[[last_key]]) && !is.null(deseq_res[[last_key]]$lres) && length(control_lst) > 0) {
        matched_genes <- intersect(control_lst, rownames(deseq_res[[last_key]]$lres))
        if (length(matched_genes) > 0) {
            out_tbl <- deseq_res[[last_key]]$lres[matched_genes, , drop = FALSE]
            outfile = paste0(outdir, "/", stag, "_control_genes_", data_usage_str, ".txt")
            write.table(out_tbl, outfile, sep = "\t")
            # Get the locus tag for all the genes
            allgenes <- rownames(out_tbl)
            gene_tab <- data.frame(matrix(ncol = 1, nrow = 0))
            for (lgene in allgenes) {
                locus_tag <- sub("^CDS:(\\S+?):.*$", "\\1", lgene)
                gene_tab[locus_tag, 1] <- '|'
            }
            homology_file = paste0(outdir, "/", stag, "_control_genes_homology_", data_usage_str, ".txt")
            write.table(gene_tab, homology_file, sep = "\t", quote = FALSE, col.names = FALSE)
        } else {
            print_log("No matching genes between control_lst and lres rownames. Output file not written.")
        }
    } else {
        print_log("Cannot write output: missing DESeq results or no control genes.")
    }
}

print_log <- function(lstr) {
    print(lstr)
    write(lstr, file = logfile, append = TRUE)
}
