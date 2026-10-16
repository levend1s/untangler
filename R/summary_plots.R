library(dplyr)
library(edgeR)
library(ggplot2)

# ============================================================
# ---- 1. Data loading ----
# ============================================================

load_counts <- function(counts_table) {
  counts <- read.delim(counts_table, row.names = 1, check.names = FALSE)
  counts <- counts[!grepl("_clusterNA$", rownames(counts)), ]
  colnames(counts) <- sub("_read_depth$", "", colnames(counts))
  counts
}

make_group <- function(counts) {
  factor(ifelse(grepl("C", colnames(counts)), "C", "K"), levels = c("C", "K"))
}

load_highlight_clusters <- function(cluster_summary) {
  df <- read.delim(cluster_summary, sep = "\t", header = TRUE, quote = "",
                   fill = TRUE, na.strings = c("", "NA"), stringsAsFactors = FALSE)
  res <- df[!is.na(df$predictor_combo) & trimws(df$predictor_combo) != "[]", ]
  data.frame(ID_cluster = res$index, stringsAsFactors = FALSE)
}

# ============================================================
# ---- 2. Proportions and edgeR diffSplice ----
# ============================================================

compute_proportions <- function(counts, gene, group) {
  gene_totals <- rowsum(counts, gene)
  prop <- as.matrix(counts) / gene_totals[gene, ]
  
  prop_res <- data.frame(
    gene       = gene,
    subcluster = rownames(counts),
    prop_C     = rowMeans(prop[, group == "C", drop = FALSE], na.rm = TRUE),
    prop_K     = rowMeans(prop[, group == "K", drop = FALSE], na.rm = TRUE)
  )
  prop_res$delta_prop <- prop_res$prop_K - prop_res$prop_C
  
  list(gene_totals = gene_totals, prop_res = prop_res)
}

run_diffsplice <- function(counts, gene, group) {
  y <- DGEList(counts = counts, group = group,
               genes = data.frame(GeneID = gene, Subcluster = rownames(counts)))
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  
  design <- model.matrix(~ group)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  ds  <- diffSpliceDGE(fit, coef = 2, geneid = "GeneID", exonid = "Subcluster")
  
  list(
    y = y, design = design, fit = fit, ds = ds,
    gene_hits  = topSpliceDGE(ds, test = "gene",  number = Inf),
    simes_hits = topSpliceDGE(ds, test = "Simes", number = Inf),
    sub_hits   = topSpliceDGE(ds, test = "exon",  number = Inf)
  )
}

# ============================================================
# ---- 3. Plot: volcano of subclusters ----
# ============================================================

plot_subcluster_volcano <- function(sub_hits, fdr_cut = 0.05, lfc_cut = 0.5, n_label = 15) {
  d <- sub_hits
  d$negLogFDR <- -log10(d$FDR)
  d$sig       <- d$FDR < fdr_cut & abs(d$logFC) >= lfc_cut
  d$label     <- paste0(d$GeneID, " ", sub("^.*_(cluster\\d+).*$", "\\1", d$Subcluster))
  
  d$direction <- factor(
    ifelse(d$sig & d$logFC > 0, "Significant, increased",
           ifelse(d$sig & d$logFC < 0, "Significant, decreased",
                  "Not significant")),
    levels = c("Not significant", "Significant, decreased", "Significant, increased"))
  d <- d[order(d$direction), ]     # draw significant points on top
  
  cat("Subclusters tested:", nrow(d), "\n")
  cat("Significant, increased:", sum(d$direction == "Significant, increased"), "\n")
  cat("Significant, decreased:", sum(d$direction == "Significant, decreased"), "\n")
  
  p <- ggplot(d, aes(x = logFC, y = negLogFDR, colour = direction)) +
    geom_point(aes(size = direction, alpha = direction)) +
    scale_colour_manual(values = c("Not significant" = "grey75",
                                   "Significant, decreased" = "steelblue",
                                   "Significant, increased" = "firebrick"),
                        name = NULL) +
    scale_size_manual(values = c(1, 1.8, 1.8), guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.9, 0.9), guide = "none") +
    geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed") +
    labs(x = "log2 fold change in subcluster proportion (K vs C, relative to rest of gene)",
         y = "-log10(FDR)",
         title = "Change in subcluster usage within genes: K vs C") +
    theme_classic()
  
  if (lfc_cut > 0) p <- p + geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed")
  
  lab <- head(d[d$sig, ][order(d$FDR[d$sig]), ], n_label)
  if (nrow(lab) > 0) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(data = lab, aes(label = label),
                                        size = 3, max.overlaps = Inf, show.legend = FALSE)
    } else {
      p <- p + geom_text(data = lab, aes(label = label), size = 3, vjust = -0.5, show.legend = FALSE)
    }
  } else {
    message("No subclusters passed the thresholds; nothing to label. Try lfc_cut <- 0 or a looser fdr_cut.")
  }
  p
}

# ============================================================
# ---- 3b. Gene-level differential expression volcano (edgeR) ----
# ============================================================

load_edger_de <- function(de_file) {
  d <- read.delim(de_file, check.names = FALSE, stringsAsFactors = FALSE)
  
  # Gene IDs: first column if it's non-numeric text, otherwise rownames
  id_col <- intersect(c("GeneID", "gene", "Gene", "gene_id", "ID", "Geneid"), names(d))
  if (length(id_col) > 0) {
    d$GeneID <- d[[id_col[1]]]
  } else if (!is.numeric(d[[1]])) {
    d$GeneID <- d[[1]]
  } else {
    d$GeneID <- rownames(d)
  }
  
  pick <- function(candidates) {
    hit <- intersect(candidates, names(d))
    if (length(hit) == 0) stop("None of these columns found: ", paste(candidates, collapse = ", "),
                               "\nAvailable: ", paste(names(d), collapse = ", "))
    hit[1]
  }
  d$logFC   <- d[[pick(c("logFC", "log2FoldChange", "log2FC", "logFC (KS / control)", "logFC (control / KS)"))]]
  d$PValue  <- d[[pick(c("PValue", "P.Value", "pvalue", "p_value"))]]
  d$FDR     <- d[[pick(c("FDR", "adj.P.Val", "padj", "qvalue"))]]
  d[!is.na(d$FDR) & !is.na(d$logFC), ]
}

plot_de_volcano <- function(de, fdr_cut = 0.05, lfc_cut = 1, n_label = 15,
                            highlight = NULL, highlight_name = "Differential splicing hit") {
  d <- de
  d$negLogFDR <- -log10(d$FDR)
  d$sig <- d$FDR < fdr_cut & abs(d$logFC) >= lfc_cut
  d$direction <- factor(
    ifelse(d$sig & d$logFC > 0, "Significant, up in K",
           ifelse(d$sig & d$logFC < 0, "Significant, down in K", "Not significant")),
    levels = c("Not significant", "Significant, down in K", "Significant, up in K"))
  d <- d[order(d$direction), ]   # significant points drawn on top
  
  cat("Genes tested:", nrow(d), "\n")
  cat("Significant, up in K:  ", sum(d$direction == "Significant, up in K"), "\n")
  cat("Significant, down in K:", sum(d$direction == "Significant, down in K"), "\n")
  
  p <- ggplot(d, aes(x = logFC, y = negLogFDR, colour = direction)) +
    geom_point(aes(size = direction, alpha = direction)) +
    scale_colour_manual(values = c("Not significant" = "grey75",
                                   "Significant, down in K" = "steelblue",
                                   "Significant, up in K" = "firebrick"),
                        name = NULL) +
    scale_size_manual(values = c(1, 1.8, 1.8), guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.9, 0.9), guide = "none") +
    geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed") +
    labs(x = "log2 fold change in gene expression (K vs C)",
         y = "-log10(FDR)",
         title = "Differential gene expression: K vs C") +
    theme_classic()
  
  if (lfc_cut > 0) p <- p + geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed")
  
  # Optional: ring genes of interest (e.g. significant diffSplice genes)
  if (!is.null(highlight)) {
    hl <- d[d$GeneID %in% highlight, ]
    if (nrow(hl) > 0) {
      p <- p + geom_point(data = hl, aes(shape = highlight_name),
                          colour = "black", size = 2.6, stroke = 0.8, fill = NA, alpha = 1) +
        scale_shape_manual(values = c(21), name = NULL)
    }
    cat("Highlighted genes found in DE table:", nrow(hl), "of", length(highlight), "\n")
  }
  
  lab <- d[d$sig, ]
  lab <- head(lab[order(lab$FDR), ], n_label)
  if (nrow(lab) > 0) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(data = lab, aes(label = GeneID), size = 3,
                                        max.overlaps = Inf, show.legend = FALSE)
    } else {
      p <- p + geom_text(data = lab, aes(label = GeneID), size = 3, vjust = -0.5, show.legend = FALSE)
    }
  } else {
    message("No genes passed the thresholds; try lfc_cut = 0 or a looser fdr_cut.")
  }
  p
}

# ============================================================
# ---- 4. GFF loading + chromosome map ----
# ============================================================

load_gff_genes <- function(annotation_file) {
  lines <- readLines(annotation_file)
  fa <- grep("^##FASTA", lines)
  if (length(fa)) lines <- lines[seq_len(fa[1] - 1)]
  
  sr <- strsplit(grep("^##sequence-region", lines, value = TRUE), "\\s+")
  chr_len <- do.call(rbind, lapply(sr, function(x)
    data.frame(seqid = x[2], len = as.numeric(x[4]))))
  
  gff <- read.delim(text = lines[!grepl("^#", lines)], header = FALSE, sep = "\t",
                    quote = "", stringsAsFactors = FALSE,
                    col.names = c("seqid", "source", "type", "start", "end",
                                  "score", "strand", "phase", "attr"))
  genes <- gff[grepl("gene$", gff$type), ]
  genes$GeneID <- sub(".*ID=([^;]+).*", "\\1", genes$attr)
  genes$pos <- (genes$start + genes$end) / 2
  
  nuc <- grepl("^Pf3D7_\\d{2}_v3$", genes$seqid)
  genes <- genes[nuc, ]
  chr_len <- chr_len[grepl("^Pf3D7_\\d{2}_v3$", chr_len$seqid), ]
  
  chr_num <- function(s) as.integer(sub("^Pf3D7_(\\d{2})_v3$", "\\1", s))
  chr_len$label <- factor(paste("Chr", chr_num(chr_len$seqid)),
                          levels = rev(paste("Chr", sort(chr_num(chr_len$seqid)))))
  genes$label <- factor(paste("Chr", chr_num(genes$seqid)), levels = levels(chr_len$label))
  
  list(genes = genes, chr_len = chr_len)
}

plot_chromosome_map <- function(gene_hits, genes, chr_len, fdr_cut = 0.05) {
  gh <- merge(gene_hits, genes[, c("GeneID", "pos", "label")], by = "GeneID")
  gh$sig <- gh$FDR < fdr_cut
  gh$neglogFDR <- -log10(gh$FDR)
  
  cat("Matched", nrow(gh), "of", nrow(gene_hits), "tested genes\n")
  
  ggplot() +
    geom_segment(data = chr_len, aes(x = 0, xend = len, y = label, yend = label),
                 linewidth = 7, colour = "grey88", lineend = "round") +
    geom_point(data = gh[!gh$sig, ],
               aes(x = pos, y = label), shape = 124, size = 3, colour = "grey55") +
    geom_point(data = gh[gh$sig, ],
               aes(x = pos, y = label, colour = neglogFDR), shape = 124, size = 6) +
    scale_colour_viridis_c(option = "inferno", direction = -1, name = "-log10(FDR)") +
    scale_x_continuous(labels = function(x) paste0(x / 1e6, " Mb"),
                       expand = expansion(mult = c(0.01, 0.02))) +
    labs(x = "Position along chromosome", y = NULL,
         title = "Genes with significant change in subcluster proportions") +
    theme_classic() +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())
}

# ============================================================
# ---- 5. Per-gene proportion shift summary ----
# ============================================================

compute_gene_shift_summary <- function(prop_res, gene_totals, group, gene_hits, simes_hits,
                                       min_reads = 20) {
  reads_C <- rowSums(gene_totals[, group == "C", drop = FALSE])
  reads_K <- rowSums(gene_totals[, group == "K", drop = FALSE])
  
  gene_sum <- prop_res %>%
    dplyr::filter(!is.na(prop_C), !is.na(prop_K)) %>%
    dplyr::group_by(gene) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::mutate(pC = prop_C / sum(prop_C),
                  pK = prop_K / sum(prop_K),
                  d  = pK - pC) %>%
    dplyr::summarise(
      n_subclusters  = dplyr::n(),
      shift          = 0.5 * sum(abs(d)),
      top_subcluster = subcluster[which.max(abs(d))],
      top_change     = d[which.max(abs(d))],
      .groups = "drop") %>%
    dplyr::filter(!is.na(shift))
  
  gene_sum$reads_C <- reads_C[gene_sum$gene]
  gene_sum$reads_K <- reads_K[gene_sum$gene]
  gene_sum <- gene_sum[pmin(gene_sum$reads_C, gene_sum$reads_K) >= min_reads, ]
  
  gene_sum$FDR       <- gene_hits$FDR[match(gene_sum$gene, gene_hits$GeneID)]
  gene_sum$FDR_simes <- simes_hits$FDR[match(gene_sum$gene, simes_hits$GeneID)]
  gene_sum$sig       <- !is.na(gene_sum$FDR) & gene_sum$FDR < 0.05
  gene_sum$negLogFDR <- -log10(gene_sum$FDR)
  gene_sum$gene.F    <- gene_hits$gene.F[match(gene_sum$gene,  gene_hits$GeneID)]
  gene_sum$P.Value   <- gene_hits$P.Value[match(gene_sum$gene, gene_hits$GeneID)]
  
  gene_sum[order(-gene_sum$shift), ]
}

# ============================================================
# ---- 6. Stacked proportion bar plots ----
# ============================================================

plot_gene_props <- function(genes_vec, title, gene_sum, prop_res) {
  info <- gene_sum[match(genes_vec, gene_sum$gene), ]
  panel_lab <- paste0(info$gene,
                      "\nshift = ", sprintf("%.2f", info$shift),
                      "\nF = ", ifelse(is.na(info$gene.F), "NA", sprintf("%.2f", info$gene.F)),
                      ", FDR = ", ifelse(is.na(info$FDR), "NA", formatC(info$FDR, format = "e", digits = 1))
  )
  names(panel_lab) <- info$gene
  
  long <- prop_res %>%
    dplyr::filter(gene %in% genes_vec) %>%
    dplyr::group_by(gene) %>%
    dplyr::mutate(C = prop_C / sum(prop_C), K = prop_K / sum(prop_K)) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_longer(c(C, K), names_to = "condition", values_to = "proportion")
  long$cluster <- sub("^.*_(cluster\\d+).*$", "\\1", long$subcluster)
  long$panel   <- factor(panel_lab[long$gene], levels = panel_lab[genes_vec])
  
  ggplot(long, aes(x = condition, y = proportion, fill = cluster)) +
    geom_col(colour = "white", linewidth = 0.2) +
    facet_wrap(~ panel, nrow = 2) +
    labs(x = NULL, y = "Share of gene's reads", fill = "Subcluster", title = title) +
    theme_classic() +
    theme(strip.text = element_text(size = 8))
}

#' Top N genes by proportion shift (descending), or a custom gene list
#'
#' @param genes  Optional character vector of gene IDs. If supplied, these
#'               genes are plotted (sorted by shift, descending) instead of
#'               the automatic top-N selection, and `n` is ignored unless
#'               the list is longer than it -- in which case it's truncated
#'               to the top `n` *within the supplied list*.
plot_top_shifted <- function(gene_sum, prop_res, n = 12, genes = NULL) {
  if (!is.null(genes)) {
    missing <- setdiff(genes, gene_sum$gene)
    if (length(missing) > 0) {
      message("Not found in gene_sum (skipped): ", paste(missing, collapse = ", "))
    }
    sel <- gene_sum[gene_sum$gene %in% genes, ]
    sel <- sel[order(-sel$shift), ]
    top_genes <- head(sel$gene, n)
    title <- "Selected genes: subcluster usage in C vs K"
  } else {
    top_genes <- head(gene_sum$gene, n)   # gene_sum is pre-sorted by shift, descending
    title <- "Most shifted genes: subcluster usage in C vs K"
  }
  plot_gene_props(top_genes, title, gene_sum, prop_res)
}

#' Bottom N genes by proportion shift among well-covered genes, or a custom gene list
#'
#' @param genes  Optional character vector of gene IDs. If supplied, these
#'               genes are plotted (sorted by shift, ascending) instead of
#'               the automatic least-shifted selection; `min_reads_stable`
#'               is ignored in that case since you've chosen the genes
#'               explicitly.
plot_least_shifted <- function(gene_sum, prop_res, n = 12, min_reads_stable = 100, genes = NULL) {
  if (!is.null(genes)) {
    missing <- setdiff(genes, gene_sum$gene)
    if (length(missing) > 0) {
      message("Not found in gene_sum (skipped): ", paste(missing, collapse = ", "))
    }
    sel <- gene_sum[gene_sum$gene %in% genes, ]
    least_genes <- head(sel$gene[order(sel$shift)], n)
    title <- "Selected genes: subcluster usage in C vs K"
  } else {
    stable <- gene_sum[pmin(gene_sum$reads_C, gene_sum$reads_K) >= min_reads_stable, ]
    least_genes <- head(stable$gene[order(stable$shift)], n)
    title <- "Least shifted genes: subcluster usage in C vs K"
  }
  print(least_genes)
  plot_gene_props(least_genes, title, gene_sum, prop_res)
}

#' Plot an arbitrary list of genes, with no ranking/filtering applied at all
#'
#' Use this instead of plot_top_shifted()/plot_least_shifted() when the genes
#' are chosen for some external reason (e.g. a predicted-m6A gene set) rather
#' than by their proportion shift, and you don't want them reordered or capped.
#'
#' @param genes_vec  character vector of gene IDs, in the order you want them
#'                   plotted (left to right, wrapping at `nrow`)
#' @param title      plot title; defaults to a generic label if not supplied
#' @param sort_by_shift  if TRUE, reorders genes_vec by shift (descending)
#'                       before plotting instead of using the given order
plot_selected_genes <- function(genes_vec, gene_sum, prop_res, title = NULL,
                                sort_by_shift = FALSE) {
  missing <- setdiff(genes_vec, gene_sum$gene)
  if (length(missing) > 0) {
    message("Not found in gene_sum (skipped): ", paste(missing, collapse = ", "))
    genes_vec <- intersect(genes_vec, gene_sum$gene)
  }
  if (length(genes_vec) == 0) {
    stop("None of the requested genes were found in gene_sum.")
  }
  
  if (sort_by_shift) {
    ord <- gene_sum$gene[gene_sum$gene %in% genes_vec]
    ord <- ord[order(-gene_sum$shift[match(ord, gene_sum$gene)])]
    genes_vec <- ord
  }
  
  if (is.null(title)) title <- "Selected genes: subcluster usage in C vs K"
  
  plot_gene_props(genes_vec, title, gene_sum, prop_res)
}

# ============================================================
# ---- 7. Single-gene inspection ----
# ============================================================

inspect_gene <- function(gid, counts, gene, prop_res, gene_hits, simes_hits, sub_hits,
                         gene_sum, ds, gene_info = NULL) {
  rows <- which(gene == gid)
  if (length(rows) == 0) { message(gid, " not found in the counts table"); return(invisible(NULL)) }
  
  cat("\n=====", gid, "=====\n")
  if (!is.null(gene_info)) {
    d <- gene_info[match(gid, gene_info$GeneID), ]
    cat("Symbol:", d$symbol, "\nDescription:", d$description, "\n")
  }
  
  cat("\n-- Raw counts per sample --\n")
  print(counts[rows, , drop = FALSE])
  
  cat("\n-- Subcluster proportions (mean per condition) --\n")
  print(prop_res[prop_res$gene == gid, c("subcluster", "prop_C", "prop_K", "delta_prop")],
        row.names = FALSE, digits = 3)
  
  cat("\n-- Gene-level tests (edgeR diffSplice) --\n")
  gh <- gene_hits[gene_hits$GeneID == gid, ]
  sh <- simes_hits[simes_hits$GeneID == gid, ]
  if (nrow(gh) == 0) {
    cat("Not tested (filtered by filterByExpr or <2 subclusters left)\n")
  } else {
    print(gh[, c("GeneID", "NExons", "gene.F", "P.Value", "FDR")], row.names = FALSE)
    cat("Simes FDR:", signif(sh$FDR, 3), "\n")
  }
  
  cat("\n-- Subcluster-level tests (vs rest of gene) --\n")
  sb <- sub_hits[sub_hits$GeneID == gid, ]
  if (nrow(sb) > 0) {
    keep_cols <- intersect(c("Subcluster", "logFC", "exon.F", "t", "P.Value", "FDR"), names(sb))
    print(sb[, keep_cols], row.names = FALSE)
  }
  
  if (gid %in% gene_sum$gene) {
    cat("\n-- Proportion shift --\n")
    print(gene_sum[gene_sum$gene == gid, c("shift", "top_subcluster", "top_change", "reads_C", "reads_K")],
          row.names = FALSE)
  }
  
  if (gid %in% ds$genes$GeneID) plotSpliceDGE(ds, geneid = gid, genecol = "GeneID")
  if (gid %in% gene_sum$gene) print(plot_gene_props(gid, gid, gene_sum, prop_res))
  
  invisible(NULL)
}

# ============================================================
# ---- 8. Shift vs CPM ----
# ============================================================
plot_shift_vs_cpm <- function(gene_sum, y, gene, fdr_cut = 0.05) {
  # Average logCPM per gene, using all samples (from the already-normalised y)
  logcpm <- cpm(y, log = TRUE, prior.count = 1)
  gene_all <- gene[match(rownames(y$counts), rownames(y$counts))]  # y is already subset/filtered
  gene_for_y <- sub("_cluster.*$", "", rownames(y$counts))
  gene_cpm <- tapply(rowMeans(logcpm), gene_for_y, sum)   # sum subcluster CPM back to gene level
  
  gene_sum$logCPM <- gene_cpm[gene_sum$gene]
  
  d <- gene_sum[!is.na(gene_sum$logCPM) & !is.na(gene_sum$FDR), ]
  d$direction <- factor(
    ifelse(d$FDR < fdr_cut & d$top_change > 0, "Significant, increased",
           ifelse(d$FDR < fdr_cut & d$top_change < 0, "Significant, decreased",
                  "Not significant")),
    levels = c("Not significant", "Significant, decreased", "Significant, increased"))
  d <- d[order(d$direction), ]   # draw significant points on top
  
  cat("Significant, increased (top subcluster gained share):", sum(d$direction == "Significant, increased"), "\n")
  cat("Significant, decreased (top subcluster lost share):",  sum(d$direction == "Significant, decreased"), "\n")
  
  ggplot(d, aes(x = logCPM, y = shift, colour = direction)) +
    geom_point(aes(size = direction, alpha = direction)) +
    scale_colour_manual(values = c("Not significant" = "grey75",
                                   "Significant, decreased" = "steelblue",
                                   "Significant, increased" = "firebrick"),
                        name = NULL) +
    scale_size_manual(values = c(1, 2, 2), guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.9, 0.9), guide = "none") +
    labs(x = "Gene abundance (log2 CPM)", y = "Proportion shift",
         title = "Proportion shift vs abundance",
         subtitle = "Colour = direction of change in the top-shifting subcluster (FDR < 0.05)") +
    theme_classic()
}

# ============================================================
# ---- 8. Diagnostic plots ----
# ============================================================

plot_depth_bias <- function(gene_sum) {
  plot(log10(pmin(gene_sum$reads_C, gene_sum$reads_K)), gene_sum$shift,
       xlab = "log10 min reads (C or K)", ylab = "Proportion shift", pch = 20, col = "grey40")
}

plot_pvalue_histogram <- function(gene_sum) {
  d <- gene_sum[!is.na(gene_sum$P.Value), ]
  ggplot(d, aes(x = P.Value)) +
    geom_histogram(breaks = seq(0, 1, by = 0.02), colour = "white", fill = "steelblue") +
    geom_hline(yintercept = nrow(d) * 0.02, linetype = "dashed", colour = "red") +
    labs(x = "p-value", y = "Number of genes",
         title = "Distribution of gene-level p-values") +
    theme_classic()
}

plot_shift_vs_fdr <- function(gene_sum, fdr_cut = 0.05, shift_cut = 0.2,
                              large_shift_only = FALSE, show_nonsig = TRUE,
                              n_label = 5) {
  gene_plot <- gene_sum[!is.na(gene_sum$FDR), ]
  gene_plot$negLogFDR <- -log10(gene_plot$FDR)
  
  gene_plot$sig <- gene_plot$FDR < fdr_cut &
    (!large_shift_only | gene_plot$shift >= shift_cut)
  
  gene_plot$status <- factor(
    ifelse(!gene_plot$sig, "Not significant",
           ifelse(gene_plot$shift >= shift_cut, "Significant & large shift",
                  "Significant, small shift")),
    levels = c("Not significant", "Significant, small shift", "Significant & large shift"))
  
  cat(sprintf("Genes tested: %d; significant: %d\n",
              nrow(gene_plot), sum(gene_plot$sig)))
  print(table(gene_plot$status))
  
  if (!show_nonsig) gene_plot <- gene_plot[gene_plot$sig, ]
  if (sum(gene_plot$sig) == 0) message("No genes pass the thresholds; try a looser fdr_cut.")
  
  gene_plot <- gene_plot[order(gene_plot$status), ]   # grey drawn first, significant on top
  
  p <- ggplot(gene_plot, aes(x = shift, y = negLogFDR, colour = status)) +
    geom_point(aes(size = status, alpha = status)) +
    scale_colour_manual(values = c("Not significant" = "grey75",
                                   "Significant, small shift" = "orange",
                                   "Significant & large shift" = "firebrick"),
                        name = NULL, drop = FALSE) +
    scale_size_manual(values = c(1, 1.8, 1.8), guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.85, 0.85), guide = "none") +
    geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed") +
    geom_vline(xintercept = shift_cut, linetype = "dashed") +
    labs(x = "Proportion shift", y = "-log10(FDR), gene-level test",
         title = "All genes: proportion shift vs significance",
         subtitle = sprintf("Significant = FDR < %.2f%s", fdr_cut,
                            if (large_shift_only) sprintf(" and shift >= %.2f", shift_cut) else "")) +
    theme_classic()
  
  lab <- gene_plot %>%
    dplyr::filter(sig) %>%
    dplyr::group_by(status) %>%
    dplyr::slice_max(negLogFDR + shift, n = n_label) %>%
    dplyr::ungroup()
  
  if (nrow(lab) > 0) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 3,
                                        max.overlaps = Inf, show.legend = FALSE)
    } else {
      p <- p + geom_text(data = lab, aes(label = gene), size = 3, vjust = -0.7,
                         show.legend = FALSE)
    }
  }
  p
}
# ============================================================
# ---- 8b. Shift vs edgeR differential expression ----
# ============================================================

plot_shift_vs_de <- function(gene_sum, de, fdr_cut = 0.05, lfc_cut = 1,
                             splice_fdr_cut = fdr_cut, use_abs_lfc = FALSE,
                             sig_only = TRUE, label_all = TRUE, n_label = 15) {
  d <- merge(gene_sum, de[, c("GeneID", "logFC", "PValue", "FDR")],
             by.x = "gene", by.y = "GeneID", suffixes = c("", "_de"))
  d <- d[!is.na(d$shift) & !is.na(d$logFC) & !is.na(d$FDR) & !is.na(d$FDR_de), ]
  
  d$de_sig     <- d$FDR_de < fdr_cut & abs(d$logFC) >= lfc_cut
  d$splice_sig <- d$FDR < splice_fdr_cut
  d$status <- factor(
    ifelse(d$splice_sig & d$de_sig, "Both",
           ifelse(d$splice_sig, "Differential splicing only",
                  ifelse(d$de_sig, "Differential expression only", "Neither"))),
    levels = c("Neither", "Differential expression only",
               "Differential splicing only", "Both"))
  
  cat(sprintf("Genes with both results: %d\n", nrow(d)))
  print(table(d$status))
  
  if (sig_only) d <- d[d$status == "Both", ]
  if (nrow(d) == 0) {
    message("No genes to plot; try lfc_cut = 0 or looser FDR cutoffs.")
    return(invisible(NULL))
  }
  d <- d[order(d$status), ]   # grey first, significant on top
  
  d$x <- if (use_abs_lfc) abs(d$logFC) else d$logFC
  d$negLogFDR_splice <- -log10(d$FDR)
  
  subtitle <- sprintf("Splicing FDR<%.2f; DE FDR<%.2f & |logFC|>=%.1f", splice_fdr_cut, fdr_cut, lfc_cut)
  if (nrow(d) >= 5) {
    ct <- suppressWarnings(cor.test(d$x, d$shift, method = "spearman"))
    cat(sprintf("Spearman rho (plotted genes) = %.3f (p = %.2g)\n", ct$estimate, ct$p.value))
    subtitle <- paste0(subtitle, sprintf("; Spearman rho = %.2f", ct$estimate))
  }
  
  p <- ggplot(d, aes(x = x, y = shift, colour = status)) +
    geom_point(aes(size = status, alpha = status)) +
    scale_colour_manual(values = c("Neither" = "grey75",
                                   "Differential expression only" = "steelblue",
                                   "Differential splicing only" = "orange",
                                   "Both" = "firebrick"),
                        name = NULL, drop = FALSE) +
    scale_size_manual(values = c(1, 1.8, 1.8, 2.5), guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.8, 0.8, 0.95), guide = "none") +
    labs(x = if (use_abs_lfc) "|log2 fold change| in gene expression (K vs C)"
         else "log2 fold change in gene expression (K vs C)",
         y = "Proportion shift",
         title = if (sig_only) "Genes significant for both subcluster shift and differential expression"
         else "Subcluster proportion shift vs differential expression (all genes)",
         subtitle = subtitle) +
    theme_classic()
  
  if (!use_abs_lfc) p <- p + geom_vline(xintercept = 0, colour = "grey60")
  if (lfc_cut > 0) {
    p <- p + geom_vline(xintercept = if (use_abs_lfc) lfc_cut else c(-lfc_cut, lfc_cut),
                        linetype = "dashed")
  }
  
  # Label only genes significant for both, unless label_all = FALSE
  lab <- d[d$status == "Both", ]
  if (!label_all) lab <- head(lab[order(lab$FDR), ], n_label)
  if (nrow(lab) > 0) {
    lab$label <- sprintf("%s\nDE p=%.1e\nsplice p=%.1e", lab$gene, lab$PValue, lab$P.Value)
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(data = lab, aes(label = label), size = 2.6,
                                        max.overlaps = Inf, show.legend = FALSE,
                                        lineheight = 0.9, colour = "black")
    } else {
      p <- p + geom_text(data = lab, aes(label = label), size = 2.6, vjust = -0.7,
                         show.legend = FALSE, colour = "black")
    }
  }
  p
}




cluster_features <- function(df,
                             features = c("logCPM", "logFC", "shift", "n_subclusters",
                                          "n_canonical_m6a", "wam_delta", "utr_delta"),
                             id_col = "GeneID",
                             method = c("kmeans", "hclust"),
                             k = NULL, k_range = 2:8, seed = 1,
                             transforms = list(
                               shift           = sqrt,
                               n_subclusters   = log2,
                               n_canonical_m6a = log1p,
                               utr_delta       = function(x) sign(x) * log1p(abs(x)))) {
  method <- match.arg(method)
  d <- as.data.frame(df)
  
  ok <- complete.cases(d[, features]) & apply(d[, features], 1, function(r) all(is.finite(r)))
  cat(sprintf("Genes used: %d of %d\n", sum(ok), nrow(d)))
  d <- d[ok, ]
  
  # Transform skewed features, then z-score
  X <- as.matrix(d[, features])
  for (f in intersect(names(transforms), features)) X[, f] <- transforms[[f]](X[, f])
  Xs <- scale(X)
  rownames(Xs) <- d[[id_col]]
  
  dm <- dist(Xs)
  hc <- if (method == "hclust") hclust(dm, method = "ward.D2") else NULL
  
  # Cluster assignment for a given k
  assign_k <- function(kk) {
    if (method == "kmeans") { set.seed(seed); kmeans(Xs, centers = kk, nstart = 50)$cluster }
    else cutree(hc, k = kk)
  }
  
  # Choose k by average silhouette if not supplied
  sil_df <- NULL
  if (is.null(k)) {
    sil_df <- do.call(rbind, lapply(k_range, function(kk)
      data.frame(k = kk, avg_sil = mean(cluster::silhouette(assign_k(kk), dm)[, 3]))))
    k <- sil_df$k[which.max(sil_df$avg_sil)]
    cat("Silhouette-selected k =", k, "\n"); print(sil_df, row.names = FALSE)
  }
  
  cl <- assign_k(k)
  d$cluster <- factor(match(cl, order(table(cl), decreasing = TRUE)))   # 1 = largest
  
  # PCA for display
  pca <- prcomp(Xs)
  d$PC1 <- pca$x[, 1]; d$PC2 <- pca$x[, 2]
  ve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  
  # Summary of raw values per cluster
  summ <- d %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(n = dplyr::n(),
                     dplyr::across(dplyr::all_of(features), median, .names = "med_{.col}"),
                     .groups = "drop")
  print(as.data.frame(summ), digits = 3)
  
  p_pca <- ggplot(d, aes(PC1, PC2, colour = cluster)) +
    geom_point(size = 1.3, alpha = 0.75) +
    labs(x = sprintf("PC1 (%.0f%%)", ve[1]), y = sprintf("PC2 (%.0f%%)", ve[2]),
         title = sprintf("Gene clusters (%s, k = %d)", method, k)) +
    theme_classic()
  
  prof <- as.data.frame(Xs) %>%
    dplyr::mutate(cluster = d$cluster) %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), mean), .groups = "drop") %>%
    tidyr::pivot_longer(-cluster, names_to = "feature", values_to = "z")
  prof$feature <- factor(prof$feature, levels = features)
  p_profile <- ggplot(prof, aes(feature, cluster, fill = z)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", z)), size = 3.3) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                         midpoint = 0, name = "Mean z") +
    labs(x = NULL, y = "Cluster", title = "Cluster profiles (mean z-score)") +
    theme_classic() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  out <- list(data = d, summary = summ, scaled = Xs, hclust = hc, silhouette = sil_df,
              p_pca = p_pca, p_profile = p_profile)
  if (!is.null(sil_df))
    out$p_silhouette <- ggplot(sil_df, aes(k, avg_sil)) + geom_line() + geom_point() +
    labs(x = "k", y = "Average silhouette width", title = "Choice of k") + theme_classic()
  out
}

umap_cluster <- function(df,
                         features = c("logCPM", "logFC", "shift", "n_subclusters",
                                      "n_canonical_m6a", "wam_delta", "utr_delta"),
                         id_col = "GeneID",
                         k = NULL, k_range = 2:10,
                         n_neighbors = 30, min_dist = 0.2, seed = 1,
                         transforms = list(
                           shift           = sqrt,
                           n_subclusters   = log2,
                           n_canonical_m6a = log1p,
                           utr_delta       = function(x) sign(x) * log1p(abs(x)))) {
  stopifnot(requireNamespace("uwot", quietly = TRUE),
            requireNamespace("cluster", quietly = TRUE))
  
  d <- as.data.frame(df)
  ok <- complete.cases(d[, features]) & apply(d[, features], 1, function(r) all(is.finite(r)))
  cat(sprintf("Genes used: %d of %d\n", sum(ok), nrow(d)))
  d <- d[ok, ]
  
  # Transform skewed features, then z-score
  X <- as.matrix(d[, features])
  for (f in intersect(names(transforms), features)) X[, f] <- transforms[[f]](X[, f])
  Xs <- scale(X)
  
  # UMAP
  set.seed(seed)
  emb <- uwot::umap(Xs, n_neighbors = n_neighbors, min_dist = min_dist)
  d$UMAP1 <- emb[, 1]; d$UMAP2 <- emb[, 2]
  
  # k-means on the embedding; choose k by silhouette if not given
  if (is.null(k)) {
    dm <- dist(emb)
    sil <- sapply(k_range, function(kk) {
      set.seed(seed)
      mean(cluster::silhouette(kmeans(emb, kk, nstart = 25)$cluster, dm)[, 3])
    })
    k <- k_range[which.max(sil)]
    cat("Silhouette-selected k =", k, "\n")
  }
  set.seed(seed)
  cl <- kmeans(emb, centers = k, nstart = 50)$cluster
  d$cluster <- factor(match(cl, order(table(cl), decreasing = TRUE)))   # 1 = largest
  
  # Cluster summary (raw values)
  summ <- d %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(n = dplyr::n(),
                     dplyr::across(dplyr::all_of(features), median, .names = "med_{.col}"),
                     .groups = "drop")
  print(as.data.frame(summ), digits = 3)
  
  p_umap <- ggplot(d, aes(UMAP1, UMAP2, colour = cluster)) +
    geom_point(size = 1.3, alpha = 0.75) +
    labs(title = sprintf("UMAP clusters (k = %d)", k)) +
    theme_classic()
  
  # UMAP coloured by each feature (z-scores capped at +/-3)
  long <- cbind(d[, c("UMAP1", "UMAP2")], as.data.frame(Xs)) %>%
    tidyr::pivot_longer(-c(UMAP1, UMAP2), names_to = "feature", values_to = "z")
  long$feature <- factor(long$feature, levels = features)
  long$z <- pmax(pmin(long$z, 3), -3)
  p_features <- ggplot(long, aes(UMAP1, UMAP2, colour = z)) +
    geom_point(size = 0.7) +
    facet_wrap(~ feature, nrow = 2) +
    scale_colour_gradient2(low = "steelblue", mid = "grey90", high = "firebrick",
                           midpoint = 0, name = "z-score") +
    theme_classic()
  
  list(data = d, summary = summ, scaled = Xs, embedding = emb,
       p_umap = p_umap, p_features = p_features)
}


# ============================================================
# ---- 11. Equivalence test for "no change in subcluster usage" ----
# ============================================================

# ============================================================
# ---- 11. Equivalence test (fast version) ----
# ============================================================

shift_from_counts <- function(Mc, Mk) {
  pc <- rowMeans(sweep(Mc, 2, pmax(colSums(Mc), 1), "/"))
  pk <- rowMeans(sweep(Mk, 2, pmax(colSums(Mk), 1), "/"))
  pc <- pc / sum(pc); pk <- pk / sum(pk)
  0.5 * sum(abs(pk - pc))
}

# B bootstrap shifts for one gene, vectorised across the B replicates
boot_shifts <- function(Mc, Mk, B) {
  boot_props <- function(M) {
    n_s <- ncol(M)
    tot <- colSums(M)
    idx <- matrix(sample.int(n_s, B * n_s, replace = TRUE), nrow = B)  # which replicate fills each slot
    acc <- matrix(0, nrow(M), B)
    for (i in seq_len(n_s)) {
      for (j in unique(idx[, i])) {
        if (tot[j] == 0) next
        bs <- which(idx[, i] == j)
        draws <- rmultinom(length(bs), tot[j], M[, j] / tot[j])   # read-level noise
        acc[, bs] <- acc[, bs] + draws / tot[j]
      }
    }
    acc / n_s
  }
  pc <- boot_props(Mc); pk <- boot_props(Mk)
  pc <- sweep(pc, 2, colSums(pc), "/")
  pk <- sweep(pk, 2, colSums(pk), "/")
  s <- 0.5 * colSums(abs(pk - pc))
  s[is.na(s)] <- 1     # a bootstrap sample with no reads in a group: treat as not equivalent
  s
}

test_shift_equivalence <- function(counts, gene, group, genes = NULL,
                                   delta = 0.1, B = 500, min_reads = 20,
                                   alpha = 0.05, seed = 1,
                                   n_cores = max(1, parallel::detectCores() - 1),
                                   B_screen = 100, screen_cut = 0.25) {
  counts <- as.matrix(counts)
  idx_by_gene <- split(seq_along(gene), gene)          # one pass, instead of gene == g per gene
  genes <- if (is.null(genes)) names(idx_by_gene) else intersect(genes, names(idx_by_gene))
  iC <- which(group == "C"); iK <- which(group == "K")
  
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed)             # reproducible across cores
  
  one_gene <- function(g) {
    M <- counts[idx_by_gene[[g]], , drop = FALSE]
    if (nrow(M) < 2) return(NULL)
    Mc <- M[, iC, drop = FALSE]; Mk <- M[, iK, drop = FALSE]
    if (min(sum(Mc), sum(Mk)) < min_reads) return(NULL)
    
    obs  <- shift_from_counts(Mc, Mk)
    boot <- boot_shifts(Mc, Mk, min(B_screen, B))      # cheap screen first
    # Only genes with a plausible chance of passing get the full B
    if (B > B_screen && mean(boot >= delta) <= screen_cut)
      boot <- c(boot, boot_shifts(Mc, Mk, B - B_screen))
    data.frame(gene = g, n_subclusters = nrow(M), shift_obs = obs,
               shift_upper = unname(quantile(boot, 1 - alpha)),
               p_equiv = (1 + sum(boot >= delta)) / (length(boot) + 1))
  }
  
  res <- parallel::mclapply(genes, one_gene, mc.cores = n_cores, mc.set.seed = TRUE)
  out <- do.call(rbind, res)
  out$FDR_equiv  <- p.adjust(out$p_equiv, method = "BH")
  out$equivalent <- out$FDR_equiv < alpha
  cat(sprintf("Genes tested: %d; confidently unchanged (shift < %.2f, FDR < %.2f): %d\n",
              nrow(out), delta, alpha, sum(out$equivalent)))
  out[order(out$p_equiv), ]
}


# ============================================================
# ---- RUN ----
# ============================================================

counts_table    <- "/Users/jlevendis/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/untangler_test/cluster_transcripts_28_introns/counts_table.tsv"
cluster_summary <- "/Users/jlevendis/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/untangler_test/cluster_transcripts_28_introns/cluster_summary.tsv"


# counts_table    <- "/Users/jlevendis/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/untangler_test/cluster_transcripts_28_m6A/counts_table.tsv"
# cluster_summary <- "/Users/jlevendis/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/untangler_test/cluster_transcripts_28_m6A/cluster_summary.tsv"

# cluster_summary <- "/Users/jlevendis/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/untangler_test/cluster_transcripts_28hpi/cluster_summary.tsv"
# counts_table    <- "/Users/jlevendis/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/untangler_test/cluster_transcripts_28hpi/counts_table.tsv"

counts_table <- "/Users/jlevendis/untangler/cluster_transcripts_PF3D7_1249300_28_m6A_90p/counts_table.tsv"
cluster_summary <- "/Users/jlevendis/untangler/cluster_transcripts_PF3D7_1249300_28_m6A_90p/cluster_summary.tsv"

counts_table <- "/Users/jlevendis/untangler/cluster_transcripts_PF3D7_1249300_28_90p/counts_table.tsv"
cluster_summary <- "/Users/jlevendis/untangler/cluster_transcripts_PF3D7_1249300_28_90p/cluster_summary.tsv"

annotation_file <- "~/Downloads/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff"

edge_de_28_chimeric <- "~/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/untangler_test/edgeR28hpichimeric.tsv"
edge_de_28_chimeric_100nt_cds <- "~/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/untangler_test/edgeR28hpichimeric100ntcds.tsv"
edge_de_28 <- "~/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/untangler_test/edgeR28hpi.tsv"

wam_28_file <- "~/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/untangler_test/wam28hpi.tsv"
utrs_28_file <- "~/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/untangler_test/utrs28hpi.tsv"

wam_28 <- read.delim(wam_28_file, row.names = 2, check.names = FALSE)
utrs_28 <- read.delim(utrs_28_file, row.names = 2, check.names = FALSE)

utrs_28$gene <- sub("\\.\\d+$", "", rownames(utrs_28))
wam_28$gene <- sub("\\.\\d+$", "", rownames(wam_28))

wam_28 <- wam_28[!duplicated(wam_28$gene), ]
rownames(wam_28) <- wam_28$gene

utrs_28 <- utrs_28[!duplicated(utrs_28$gene), ]
rownames(utrs_28) <- utrs_28$gene

count_list_items <- function(x) {
  x <- trimws(gsub("^\\[|\\]$", "", x))              # drop the brackets
  n <- lengths(strsplit(x, ",", fixed = TRUE))        # count comma-separated items
  n[is.na(x) | x == ""] <- 0L                         # "[]" or empty -> 0
  n[is.na(x)] <- NA_integer_                          # real NA stays NA
  n
}

wam_28$n_canonical_m6a <- count_list_items(wam_28$canonical_mods)
wam_28$wam_delta <- wam_28$canonical_wam_g2 - wam_28$canonical_wam_g1
utrs_28$utr_delta <- utrs_28$mean_poly_a_length_g2 - utrs_28$mean_utr_length_g1

counts  <- load_counts(counts_table)
gene    <- sub("_cluster.*$", "", rownames(counts))
group   <- make_group(counts)
new_df  <- load_highlight_clusters(cluster_summary)

pr <- compute_proportions(counts, gene, group)
gene_totals <- pr$gene_totals
prop_res    <- pr$prop_res

de <- run_diffsplice(counts, gene, group)
y <- de$y; fit <- de$fit; ds <- de$ds
gene_hits <- de$gene_hits; simes_hits <- de$simes_hits; sub_hits <- de$sub_hits

gff <- load_gff_genes(annotation_file)
plot_chromosome_map(gene_hits, gff$genes, gff$chr_len)

gene_sum <- compute_gene_shift_summary(prop_res, gene_totals, group, gene_hits, simes_hits)
gene_table <- gene_sum[, c("gene", "n_subclusters", "reads_C", "reads_K",
                           "shift", "gene.F", "P.Value", "FDR",
                           "top_subcluster", "top_change")]

eq <- test_shift_equivalence(counts, gene, group, genes = gene_sum$gene, delta = 0.2)
head(eq)

# Add to gene_sum
gene_sum$FDR_equiv  <- eq$FDR_equiv[match(gene_sum$gene, eq$gene)]
gene_sum$equivalent <- eq$equivalent[match(gene_sum$gene, eq$gene)]

# Compare with the original test
table(splice_sig = gene_sum$FDR < 0.05, equivalent = gene_sum$equivalent, useNA = "ifany")




chimeric_counts <- load_edger_de(edge_de_28_chimeric)

chimeric_counts$chimeric_logFC <- -chimeric_counts$logFC
chimeric_counts$chimeric_logCPM <- chimeric_counts$logCPM

# gene_sum <- gene_sum[!gene_sum$gene %in% rownames(chimeric_counts), ]

# diagnostic plots
plotBCV(y)
plot_depth_bias(gene_sum)
plot_pvalue_histogram(gene_sum)
plot_shift_histogram(gene_sum)

# cool plots
plot_top_shifted(gene_sum, prop_res)
plot_least_shifted(gene_sum, prop_res)
plot_shift_vs_fdr(gene_sum) # the more I think about this, the more irrelevant the FDR metric is. The null hypothesis is that there isn't a change?
plot_subcluster_volcano(sub_hits)
plot_shift_vs_cpm(gene_sum, y, gene)


# edgeR differential expression volcano
de_28 <- load_edger_de(edge_de_28)
# plot_de_volcano(de_28)
de_28 <- de_28[!de_28$GeneID %in% chimeric_counts$GeneID, ]
plot_de_volcano(de_28)
plot_de_volcano(chimeric_counts)

# Optional: ring genes that also show significant differential splicing
# ds_sig_genes <- gene_hits$GeneID[gene_hits$FDR < 0.05]
ds_sig_genes <- gene_sum$gene[gene_sum$shift >= 0.3]
de_volcano <- de_28#[de_28$GeneID %in% ds_sig_genes, ]
plot_de_volcano(de_28, highlight = single_genes)

ds_sig_genes <- gene_sum$gene[gene_sum$shift < 0.1 ]
de_volcano <- de_28#[de_28$GeneID %in% ds_sig_genes, ]
plot_de_volcano(chimeric_counts, highlight = ds_sig_genes)

# genes with a single cluster are excluded from shift analysis
single_genes <- names(which(table(gene) == 1 & rowSums(gene_totals[, group == "C", drop = FALSE])[names(table(gene))] >= 20 & rowSums(gene_totals[, group == "K", drop = FALSE])[names(table(gene))] >= 20))


merged <- merge(
  de_28[, c("GeneID", "logCPM", "logFC"), drop = FALSE],
  gene_sum[, c("gene", "shift", "n_subclusters")],
  by.x = "GeneID", by.y = "gene"
)

merged <- merge(
  merged,
  wam_28[, c("gene", "n_canonical_m6a", "wam_delta")],
  by.x = "GeneID", by.y = "gene"
)

merged <- merge(
  merged,
  utrs_28[, c("gene", "utr_delta")],
  by.x = "GeneID", by.y = "gene"
)

# merged <- merge(
#   merged,
#   chimeric_counts[, c("GeneID", "chimeric_logCPM", "chimeric_logFC")],
#   by.x = "GeneID", by.y = "GeneID"
# )


# plot_shift_vs_de(gene_sum, de_28, sig_only = FALSE, label_all = FALSE, n_label = 0)
# plot_shift_vs_de(gene_sum, de_28, lfc_cut = 0)               # DE by FDR only
# plot_shift_vs_de(gene_sum, de_28, label_all = FALSE, n_label = 10)  # if crowded

# TODO: DE abundance change for 0 shift genes with only 1 m6A 


# res <- umap_cluster(merged)
# res$p_umap
# res$p_features
# res$data$GeneID[res$data$cluster == 1]
# 
# # umap_cluster(merged, k = 5)
# # umap_cluster(merged, features = c("logCPM", "logFC", "shift", "n_subclusters", "n_canonical_m6a"))
# umap_cluster(merged)


# 
# res <- cluster_features(merged)
# res$p_pca
# res$p_profile
# res$p_silhouette
# 
# cluster_features(merged, k = 4)
# cluster_features(merged, method = "hclust", k = 4)
# plot(cluster_features(merged, method = "hclust", k = 4)$hclust, labels = FALSE)  # dendrogram
# 
# res$data$GeneID[res$data$cluster == 1]    # genes in a cluster


inspect_gene("PF3D7_0603200", counts, gene, prop_res, gene_hits, simes_hits, sub_hits, gene_sum, ds)


m6A_intron_predicted_genes <- c("PF3D7_0708400", "PF3D7_0832200", "PF3D7_0929000", "PF3D7_1249300", "PF3D7_0309000", "PF3D7_0715500", "PF3D7_1029300", "PF3D7_1326300", "PF3D7_0526200", 
                                "PF3D7_0820600", "PF3D7_0924700", "PF3D7_1138350", "PF3D7_1231700", "PF3D7_1361500")

# Option B: use the dedicated function for an arbitrary gene list, unranked
plot_selected_genes(m6A_intron_predicted_genes, gene_sum, prop_res, title = "Predicted m6A intron genes: subcluster usage in C vs K")
