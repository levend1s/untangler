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
                      "\nF = ", ifelse(is.na(info$gene.F), "NA", sprintf("%.2f", info$gene.F))
                      # ", FDR = ", ifelse(is.na(info$FDR), "NA", formatC(info$FDR, format = "e", digits = 1))
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

plot_top_shifted <- function(gene_sum, prop_res, n = 12) {
  top_genes <- head(gene_sum$gene, n)   # gene_sum is pre-sorted by shift, descending
  plot_gene_props(top_genes, "Most shifted genes: subcluster usage in C vs K", gene_sum, prop_res)
}

plot_least_shifted <- function(gene_sum, prop_res, n = 12, min_reads_stable = 100) {
  stable <- gene_sum[pmin(gene_sum$reads_C, gene_sum$reads_K) >= min_reads_stable, ]
  least_genes <- head(stable$gene[order(stable$shift)], n)
  plot_gene_props(least_genes, "Least shifted genes: subcluster usage in C vs K", gene_sum, prop_res)
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

plot_shift_histogram <- function(gene_sum) {
  d <- gene_sum[!is.na(gene_sum$shift), ]
  ggplot(d, aes(x = shift)) +
    geom_histogram(breaks = seq(0, 1, by = 0.02), colour = "white", fill = "steelblue") +
    geom_hline(yintercept = nrow(d) * 0.02, linetype = "dashed", colour = "red") +
    labs(x = "p-value", y = "Number of genes",
         title = "Distribution of shift values") +
    theme_classic()
}

plot_shift_vs_fdr <- function(gene_sum, fdr_cut = 0.05, shift_cut = 0.2) {
  gene_plot <- gene_sum[!is.na(gene_sum$FDR), ]
  gene_plot$negLogFDR <- -log10(gene_plot$FDR)
  gene_plot$status <- factor(
    ifelse(gene_plot$FDR < fdr_cut & gene_plot$shift >= shift_cut, "Significant & large shift",
           ifelse(gene_plot$FDR < fdr_cut,                                "Significant, small shift",
                  ifelse(gene_plot$shift >= shift_cut,                           "Large shift, not significant",
                         "Not significant"))),
    levels = c("Not significant", "Large shift, not significant",
               "Significant, small shift", "Significant & large shift"))
  
  print(table(gene_plot$status))
  
  p <- ggplot(gene_plot, aes(x = shift, y = negLogFDR, colour = status)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_colour_manual(values = c("Not significant" = "grey75",
                                   "Large shift, not significant" = "steelblue",
                                   "Significant, small shift" = "orange",
                                   "Significant & large shift" = "firebrick"),
                        name = NULL) +
    geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed") +
    geom_vline(xintercept = shift_cut, linetype = "dashed") +
    labs(x = "Proportion shift", y = "-log10(FDR), gene-level test",
         title = "All genes: proportion shift vs significance") +
    theme_classic()
  
  lab <- gene_plot %>%
    dplyr::group_by(status) %>%
    dplyr::slice_max(negLogFDR + shift, n = 5) %>%
    dplyr::ungroup()
  
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 3,
                                      max.overlaps = Inf, show.legend = FALSE)
  } else {
    p <- p + geom_text(data = lab, aes(label = gene), size = 3, vjust = -0.7, show.legend = FALSE)
  }
  p
}



# ============================================================
# ---- RUN ----
# ============================================================

# counts_table    <- "/Users/jlevendis/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/untangler_test/cluster_transcripts_28_introns/counts_table.tsv"
# cluster_summary <- "/Users/jlevendis/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/untangler_test/cluster_transcripts_28_introns/cluster_summary.tsv"


counts_table    <- "/Users/jlevendis/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/untangler_test/cluster_transcripts_28_m6A/counts_table.tsv"
cluster_summary <- "/Users/jlevendis/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/untangler_test/cluster_transcripts_28_m6A/cluster_summary.tsv"

annotation_file <- "~/Downloads/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff"

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

# diagnostic plots
plotBCV(y)
plot_depth_bias(gene_sum)
plot_pvalue_histogram(gene_sum)
plot_shift_histogram(gene_sum)

# cool plots
plot_top_shifted(gene_sum, prop_res)
plot_least_shifted(gene_sum, prop_res)
plot_shift_vs_fdr(gene_sum)
plot_subcluster_volcano(sub_hits)
plot_shift_vs_cpm(gene_sum, y, gene)

inspect_gene("PF3D7_0914400", counts, gene, prop_res, gene_hits, simes_hits, sub_hits, gene_sum, ds)
