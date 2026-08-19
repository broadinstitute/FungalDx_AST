get_sign <- function(var) {
    if (var < 0) 
        return ("neg")
    else if (var >= 0)
        return ("nonneg")
    else {
        print("problem is getting sign")
    }
}


get_splitted_parts <- function(lstr, reg_str) {
    part_lst <- strsplit(lstr, reg_str)
    parts <- part_lst[[1]]
    return (parts)
}

combine_parts <- function(sample_ids, strains, treatments, times) {
    combined_res_pre <- apply(expand.grid(strains, treatments, times), 1, paste, collapse="_")
    combined_res_lst <- lapply(combined_res_pre,
        function(pre, sample_ids) {
            matched_ones <- sample_ids[startsWith(sample_ids, pre)]
            return (matched_ones)
         }, sample_ids)


    combined_res <- unlist(combined_res_lst)
    return (combined_res)

}

get_tp_parts <- function(sample_ids, time_pos) {
    vals <- lapply(sample_ids, function(x) {strsplit(x, '_')[[1]][time_pos]})
    tp_parts <- sort(unique(unlist(vals)))
    return (tp_parts)
}

get_samples <- function(sample_ids, pos, val) {
  #adding in this IF statement to prevent the need from switching S1,S2,etc.. to all just S
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

get_sample_groups <- function(sample_ids) {

    # A sample: Ebc_S_gent_+_010m_rep1
    # Strain types : S & R
    # treatment types : + & -

    print("sample_ids")
    print(sample_ids)

    straintype_pos <- 2
    treatment_pos <- 4
    time_pos <- 5
    sus_samples <- get_samples(sample_ids, straintype_pos, 'S')
    print("sus_samples: ")
    print(sus_samples)
    res_samples <- get_samples(sample_ids, straintype_pos, 'R')
    print("res_samples: ")
    print(res_samples)
    treated_samples <- get_samples(sample_ids, treatment_pos, '+')
    print("treated_samples: ")
    print(treated_samples)
    untreated_samples <- get_samples(sample_ids, treatment_pos, '-')
    print("untreated_samples")
    print(untreated_samples)
    tp_parts <- get_tp_parts(sample_ids, time_pos)
    print("tp_parts: ")
    print(tp_parts)


    retval <- list()

    tp_len <- length(tp_parts)
     for (j in 1:tp_len) {
        ltime_point <- tp_parts[j]
        ltime_samples <- get_samples(sample_ids, time_pos, ltime_point)

        treated_str <- paste0("sus_treated_time", j)
        lsamples_treated <- list(sus_samples, treated_samples, ltime_samples)
        retval[[treated_str]] <- Reduce(intersect, lsamples_treated)

        untreated_str <- paste0("sus_untreated_time", j)
        lsamples_untreated <- list(sus_samples, untreated_samples, ltime_samples)
        retval[[untreated_str]] <- Reduce(intersect, lsamples_untreated)

        res_treated_str <- paste0("res_treated_time", j)
        lsamples_res_treated <- list(res_samples, treated_samples, ltime_samples)
        retval[[res_treated_str]] <- Reduce(intersect, lsamples_res_treated)


        res_untreated_str <- paste0("res_untreated_time", j)
        lsamples_res_untreated <- list(res_samples, untreated_samples, ltime_samples)
        retval[[res_untreated_str]] <- Reduce(intersect, lsamples_res_untreated)
    }
    
    exp_conds = list("treated_parts" = c("+"), "untreated_parts" = c("-"), "tp_parts" = tp_parts)
    retval$"exp_conds" = exp_conds
    return (retval)
}

print_map_tbls <- function(lres, stag, contrast_val, outdir){

    prefix <- paste0("stag: ", stag, "\ncontrast: ", contrast_val)
    outfile_MA <- paste0(outdir, "/MA_plots/", contrast_val, "_MA.pdf")
    pdf(outfile_MA)
    plotMA(lres$"lres", main=prefix, ylim=c(-5,5))
    dev.off()

    outfile_tbl <- paste0(outdir, "/DESeq_tbls/", contrast_val, "_tbl.tsv")
    write.table(as.data.frame(lres$"lres"), outfile_tbl, sep = "\t")

}


collapse_samples <- function(lcond_samples) {
    parts <- strsplit(lcond_samples, "_")
    lcond_len <- length(parts[[1]])
    lcond_len_2 <- lcond_len -1
    parts[[1]][2] <- substr(parts[[1]][2],1,1)
    retval <- paste(parts[[1]][1:lcond_len_2], collapse = "_")
    return (retval)
}

get_contrast_str <- function(sample_groups, cond2, cond1) {

    print("from get_contrast_str")
    lcond2_samples <- sample_groups[[cond2]]
    lcond2_str <- collapse_samples(lcond2_samples)
    lcond1_samples <- sample_groups[[cond1]]
    lcond1_str <- collapse_samples(lcond1_samples)
    lstr <- paste0(lcond2_str, "_VS_", lcond1_str)
    return (lstr)
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


get_repeated_vals <- function(repeat_tag, sample_groups) {
    tp_len <- length(sample_groups$"exp_conds"$"tp_parts")
    timepoints <- tp_len
    condition <- c()
    samples <- c()
    for (j in 1:timepoints) {
        rep_str <- paste0(repeat_tag, j)
        samples_arr <- unlist(sample_groups[[rep_str]])
        rep_arr <- rep(rep_str, length(samples_arr))
        condition <- c(condition, rep_arr)
        samples <- c(samples, samples_arr)
    }

    retval <- list("condition" = condition, "samples" = samples)
    return (retval)
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
    rowlen_lim <- floor(rowlen*(100.0 - base_lim)/100.0)
    
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
    retval <- list("lres_top" = lres_top, "lres" = lres, "baseMean_lim_val" = baseMean_lim_val)
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

print_logf <- function(lstr) {
    write(lstr, file = logfile, append = TRUE)
}

print_log <- function(lstr) {
    print(lstr)
    write(lstr, file = logfile, append = TRUE)
}


deseq_condwise_part <- function(countData, sample_groups, cond2, cond1, lfcth, padjth, altH = "greaterAbs", use_beta_prior = FALSE, base_lim = 50) {
    colData_s <- get_colData_small(sample_groups, cond2, cond1)
    lsamples <- rownames(colData_s)
    countData_s <- countData[, lsamples]
    print(colData_s)
    print_log(paste0("Info: beta prior is set to ", use_beta_prior))
    dds <- NULL
    if (use_beta_prior) {
        dds <- DESeqDataSetFromMatrix(countData = countData_s, colData = colData_s, design = ~ condition)
    } else {
        dds <- DESeqDataSetFromMatrix(countData = countData_s, colData = colData_s, design = ~ 0 + condition)
    }
    dds <- DESeq(dds, betaPrior = use_beta_prior)
    retval <- deseq_condwise(dds, cond2, cond1, lfcth, padjth, altH, base_lim)
    return (retval)
}


get_lst <- function(lst, end_str) {
    lnames <- names(lst)
    reg_str <- paste0(end_str, "$")
    lvals <- grepl(reg_str, lnames)
    ret_lst <- lst[lvals]
    return(ret_lst)
}

check_for_drop <- function(lgene, res_lst, basemean_lst, treated_term, untreated_tag, tp_len) {
  
    pval_term <- "padj"
    pval_term_raw <- "pvalue"
    logFC_sign_j <- NA
    success_j <- FALSE
    lgene_cond <- paste0(lgene, "__", treated_term)
    pval_arr <- c()
    pval_raw_arr <- c()
    basemean_arr <- c()
    logFC_arr <- c()

    for (k in 1:tp_len) {
        untreated_term <- paste0(untreated_tag, k)
        name_term = paste0(treated_term, "__", untreated_term)
        ldata <- res_lst[[name_term]]
        basemean_lim_val <- basemean_lst[[name_term]]
        basemean_k <- ldata[lgene, "baseMean"]
        PValue_k <- ldata[lgene, pval_term]
        PValue_raw_k <- ldata[lgene, pval_term_raw]
        logFC_k <- ldata[lgene, "log2FoldChange"]

        drop_tag <- paste0("lgene: ", lgene, " treated: ", treated_term, " untreated: ", untreated_term)
        if (is.na(logFC_k)) {
            success_j <- FALSE
            err_str <- paste0("Dropped: NA val of logFC, ", drop_tag)
            print_logf(err_str)
            break
        }

        if (is.na(PValue_k)) {
            success_j <- FALSE
            err_str <- paste0("Dropped: NA val of adjusted PValue, ", drop_tag)
            print_logf(err_str)
            break
        }
        
        if (is.na(basemean_k)) {
            success_j <- FALSE
            err_str <- paste0("Dropped: NA val of basemean, ", drop_tag)
            print_logf(err_str)
            break
        }

        if (is.na(PValue_raw_k)) {
            success_j <- FALSE
            err_str <- paste0("Dropped: NA val of unadjusted p value, ", drop_tag)
            print_logf(err_str)
            break
        }
         
        if (basemean_k < basemean_lim_val) {
            success_j <- FALSE
            err_str <- paste0("Dropped: basemean ", basemean_k, " < basemean_lim_val: ", basemean_lim_val, " , ", drop_tag)
            print_logf(err_str)
            break
        }

        if (PValue_raw_k > pval_lim_raw) {
            err_str <- paste0("Dropped: raw p-value ", PValue_raw_k, " >  raw p-value limit ", pval_lim_raw, " , ", drop_tag)
            print_logf(err_str)
            success_j <- FALSE
            break
        }

        logFC_sign_k <- get_sign(logFC_k)

        if (k == 1) {
            logFC_sign_j <- logFC_sign_k
            success_j <- TRUE
            pval_arr <- c(pval_arr, PValue_k)
            pval_raw_arr <- c(pval_raw_arr, PValue_raw_k)
            basemean_arr <- c(basemean_arr, basemean_k)
            logFC_arr <- c(logFC_arr, logFC_k)

        } else {
            if (logFC_sign_j != logFC_sign_k) {
                errstr <- paste0("Dropped: change of sign og logFC from ", logFC_sign_j, " to ", logFC_sign_k, " , ", drop_tag)
                print_logf(errstr)
                success_j <- FALSE
                break
            } else {
                pval_arr <- c(pval_arr, PValue_k)
                pval_raw_arr <- c(pval_raw_arr, PValue_raw_k)
                basemean_arr <- c(basemean_arr, basemean_k)
                logFC_arr <- c(logFC_arr, logFC_k)
            }
        }
    }

    ret_lst <- list("success_j" = success_j, "pval_arr" = pval_arr, "pval_raw_arr" = pval_raw_arr, "basemean_arr" = basemean_arr, "logFC_arr" = logFC_arr)
    return (ret_lst)

}


combine_arrs <- function(lst1, lst2, lst3, select_str, use_first) {
    arr1 <- lst1[[select_str]]
    if (use_first) {
        return (arr1)
    } else {
        arr2 <- lst2[[select_str]]
        arr3 <- lst3[[select_str]]
        arr_combined <- c(arr1, arr2, arr3)
        return (arr_combined)
    }
}

get_full_str <- function(ret_lst, type_str) {

    pvals_tag <- NA
    pvals_raw_tag <- NA
    logFCs_tag <- NA
    basemeans_tag <- NA

    if (type_str == "")
    {
        pvals_tag <- "pvals"
        pvals_raw_tag <- "pvals_raw"
        logFCs_tag <- "logFCs"
        basemeans_tag <- "basemeans"
    } else {
        pvals_tag <- paste0("pvals_", type_str)
        pvals_raw_tag <- paste0("pvals_raw_", type_str)
        logFCs_tag <- paste0("logFCs_", type_str)
        basemeans_tag <- paste0("basemeans_", type_str)
    }

    pval_arr <- ret_lst$"pval_arr"
    pval_raw_arr <- ret_lst$"pval_raw_arr"
    basemean_arr <- ret_lst$"basemean_arr"
    logFC_arr <- ret_lst$"logFC_arr"

    pval_arr_str <- paste(pval_arr, collapse = ", ")
    pval_raw_arr_str <- paste(pval_raw_arr, collapse = ", ")
    logFC_arr_str <- paste(logFC_arr, collapse = ", ")
    basemean_arr_str <- paste(basemean_arr, collapse = ", ")

    full_str <- paste0(pvals_tag, ": ", pval_arr_str, " ", pvals_raw_tag, ": ", pval_raw_arr_str, " ", logFCs_tag, ": ", logFC_arr_str, " ", basemeans_tag, ": ", basemean_arr_str)
    
}


select_genes <- function(gene_lst, tp_len, res_lst, basemean_lst, pval_lim_raw, treated_start, use_res) {

    count <- 1
    gene_cond_lst <- list()
    gene_cond_raw_lst <- list()

    gene_lst_len <- length(gene_lst)
    total_count <- 0
    last_print <- -1
    for (lgene in gene_lst) {
        for (j in treated_start:tp_len) {
            
            treated_term <- paste0('sus_treated_time', j)

            untreated_tag <- "sus_untreated_time"
            ret_j_su <- check_for_drop(lgene, res_lst, basemean_lst, treated_term, untreated_tag, tp_len)
            success_j_su <- ret_j_su$"success_j"

            ret_j_rt <- NA
            ret_j_ru <- NA
            success_j_rt <- NA
            success_j_ru <- NA 
            if (use_res) {
                untreated_tag <- "res_treated_time"
                ret_j_rt <- check_for_drop(lgene, res_lst, basemean_lst, treated_term, untreated_tag, tp_len)
                untreated_tag <- "res_untreated_time"
                ret_j_ru <- check_for_drop(lgene, res_lst, basemean_lst, treated_term, untreated_tag, tp_len)
                success_j_rt <- ret_j_rt$"success_j"
                success_j_ru <- ret_j_ru$"success_j"
            }

            success_j_com <- NA
            if (use_res) {
                success_j_com <- success_j_su && success_j_rt && success_j_ru
            } else {
                success_j_com <- success_j_su
            }

            if (success_j_com) {
                
                success_str <- paste0("count: ", count, " gene: ", lgene, " cond: ", treated_term)
                print_logf(success_str)
                lgene_cond <- paste0(lgene, "__", treated_term)

                use_first <- NA
                if (use_res) {
                    use_first <- FALSE
                } else {
                    use_first <- TRUE
                }

                pval_arr <- combine_arrs(ret_j_su, ret_j_rt, ret_j_ru, "pval_arr", use_first)
                gene_cond_lst[[lgene_cond]] <- pval_arr
                
                pval_raw_arr <- combine_arrs(ret_j_su, ret_j_rt, ret_j_ru, "pval_raw_arr", use_first)
                gene_cond_raw_lst[[lgene_cond]] <- pval_raw_arr

                # Since a gene at a specific time point is selected, collected 
                # most of the useful information

                # 1. Adjusted p values
                # 2. raw p-values
                # 3. log fold change values
                # 4. basemean values
                all_vals <- NA
                if (!use_res) {
                    full_str_su <- get_full_str(ret_j_su, "")
                    all_vals <- paste0("count: ", count, " gene: ", lgene, " cond: ", treated_term, " ", full_str_su)
                } else {
                    full_str_su <- get_full_str(ret_j_su, "su")
                    full_str_rt <- get_full_str(ret_j_rt, "rt")
                    full_str_ru <- get_full_str(ret_j_ru, "ru")

                    all_vals <- paste0("count: ", count, " gene: ", lgene, " cond: ", treated_term, " ", full_str_su, " ", full_str_rt, " ", full_str_ru)
                }
                write(all_vals, file = pval_logfile, append = TRUE)

                count <- count + 1
                
            }
        }

        total_count <- total_count + 1
        total_comp <- (total_count / gene_lst_len) *100
        total_comp_floor <- floor(total_comp)
        if (total_comp_floor %% 5 == 0) {
            if (total_comp_floor != last_print) {
                out_str <- paste0("Gene selection complete: ", total_comp_floor, "%")
                print_log(out_str)
                last_print <- total_comp_floor
            }
        }
    }

    ret_lst <- list("gene_cond_lst" = gene_cond_lst, "gene_cond_raw_lst" = gene_cond_raw_lst)
    return (ret_lst)
}


get_meta_pval <- function(lpvals, method) { 
   min_val <- .Machine$double.xmin
   lpvals[lpvals == 0] <- min_val
   val <- NA 
   if (method == "fisher") { 
       val <- sumlog(lpvals) 
   } else if (method == "stouffer") { 
       val <- sumz(lpvals) 
   } 
   return (val$p)
}


get_sample_groups_TB <- function(sample_ids, sus, treated, untreated, timepoins) {

    reg_str = "\\s*[,|;]\\s*"

    sus_parts <- get_splitted_parts(sus, reg_str)
    treated_parts <- get_splitted_parts(treated, reg_str)
    untreated_parts <- get_splitted_parts(untreated, reg_str)
    tp_parts <- get_splitted_parts(timepoints, reg_str)

    # Create prefixes for differenet categories

    retval <- list()

    tp_len <- length(tp_parts)
    for (j in 1:tp_len) {
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

