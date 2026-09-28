library(dplyr)
library(edgeR)
library(ggplot2)

# counts_table    <- "/Users/joshualevendis/rqc/cluster_transcripts_28_introns/counts_table.tsv"
# cluster_summary <- "/Users/joshualevendis/rqc/cluster_transcripts_28_introns/cluster_summary.tsv"  # check this path

counts_table    <- "/Users/joshualevendis/rqc/cluster_transcripts_28_m6A/counts_table.tsv"
cluster_summary <- "/Users/joshualevendis/rqc/cluster_transcripts_28_m6A/cluster_summary.tsv"  # check this path


# ---- Load counts ----
counts <- read.delim(counts_table, row.names = 1, check.names = FALSE)
counts <- counts[!grepl("_clusterNA$", rownames(counts)), ]
colnames(counts) <- sub("_read_depth$", "", colnames(counts))

gene <- sub("_cluster.*$", "", rownames(counts))
group <- factor(ifelse(grepl("C", colnames(counts)), "C", "K"), levels = c("C", "K"))

# ---- Clusters to highlight (those with predictor combos) ----
df <- read.delim(cluster_summary, sep = "\t", header = TRUE, quote = "",
                 fill = TRUE, na.strings = c("", "NA"), stringsAsFactors = FALSE)
res <- df[!is.na(df$predictor_combo) & trimws(df$predictor_combo) != "[]", ]
new_df <- data.frame(ID_cluster = res$index, stringsAsFactors = FALSE)   # index = gene_clusterN_features

# ---- Descriptive proportion change ----
gene_totals <- rowsum(counts, gene)
prop <- as.matrix(counts) / gene_totals[gene, ]

prop_res <- data.frame(
  gene       = gene,
  subcluster = rownames(counts),
  prop_C     = rowMeans(prop[, group == "C", drop = FALSE], na.rm = TRUE),
  prop_K     = rowMeans(prop[, group == "K", drop = FALSE], na.rm = TRUE)
)
prop_res$delta_prop <- prop_res$prop_K - prop_res$prop_C

# ---- edgeR: differential proportion within genes ----
y <- DGEList(counts = counts, group = group,
             genes = data.frame(GeneID = gene, Subcluster = rownames(counts)))
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

design <- model.matrix(~ group)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

ds <- diffSpliceDGE(fit, coef = 2, geneid = "GeneID", exonid = "Subcluster")

gene_hits <- topSpliceDGE(ds, test = "gene",  number = Inf)
sig_genes <- gene_hits[gene_hits$FDR < 0.05, ]

simes_hits <- topSpliceDGE(ds, test = "Simes", number = Inf)
sig_simes  <- simes_hits[simes_hits$FDR < 0.05, ]

sub_hits <- topSpliceDGE(ds, test = "exon", number = Inf)
sig_subs <- sub_hits[sub_hits$FDR < 0.05, ]
sig_subs <- merge(sig_subs, prop_res, by.x = "Subcluster", by.y = "subcluster", all.x = TRUE)

# ---- Plot a top gene (set gene_to_plot manually if you want a specific one) ----
# gene_to_plot <- if (nrow(sig_genes) > 0) sig_genes$GeneID[1] else gene_hits$GeneID[1]
# # gene_to_plot <- "PF3D7_0603200"
# if (gene_to_plot %in% ds$genes$GeneID) {
#   plotSpliceDGE(ds, geneid = gene_to_plot, genecol = "GeneID")
# } else {
#   message(gene_to_plot, " was not tested (filtered out or single subcluster)")
# }

# ---- Volcano plot of subclusters ----
fdr_cut <- 0.05
lfc_cut <- 0.5      # set to 0 to call significance on FDR alone
n_label <- 15

short_key <- function(x) sub("^(.*_cluster\\d+).*$", "\\1", x)

d <- sub_hits
d$negLogFDR <- -log10(d$FDR)
d$in_new_df <- short_key(d$Subcluster) %in% short_key(new_df$ID_cluster)
d$sig       <- d$FDR < fdr_cut & abs(d$logFC) >= lfc_cut
d$label     <- paste0(d$GeneID, " ", sub("^.*_(cluster\\d+).*$", "\\1", d$Subcluster))

d$status <- factor(
  ifelse(d$sig & d$in_new_df, "Significant, has predictor combo",
         ifelse(d$sig,               "Significant",
                ifelse(d$in_new_df,         "Predictor combo, not significant",
                       "Not significant"))),
  levels = c("Not significant", "Predictor combo, not significant",
             "Significant", "Significant, has predictor combo"))
d <- d[order(d$status), ]     # draw significant points on top

# Diagnostics
cat("Subclusters tested:", nrow(d), "\n")
cat("In new_df:", sum(d$in_new_df), "of", nrow(new_df), "\n")
cat("Significant (FDR <", fdr_cut, "& |logFC| >=", lfc_cut, "):", sum(d$sig), "\n")
cat("Significant AND in new_df:", sum(d$sig & d$in_new_df), "\n")

p <- ggplot(d, aes(x = logFC, y = negLogFDR, colour = status)) +
  geom_point(aes(size = status, alpha = status)) +
  scale_colour_manual(values = c("Not significant" = "grey75",
                                 "Predictor combo, not significant" = "#F4A6A0",
                                 "Significant" = "steelblue",
                                 "Significant, has predictor combo" = "firebrick"),
                      name = NULL) +
  scale_size_manual(values = c(1, 1.8, 1.8, 2.5), guide = "none") +
  scale_alpha_manual(values = c(0.5, 0.9, 0.9, 1), guide = "none") +
  geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed") +
  { if (lfc_cut > 0) geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") } +
  labs(x = "log2 fold change in subcluster proportion (K vs C, relative to rest of gene)",
       y = "-log10(FDR)",
       title = "Change in subcluster usage within genes: K vs C") +
  theme_classic()

# Label the most significant subclusters (whether or not they have a predictor combo)
lab <- head(d[d$sig, ][order(d$FDR[d$sig]), ], n_label)
if (nrow(lab) > 0) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(data = lab, aes(label = label),
                                      size = 3, max.overlaps = Inf, show.legend = FALSE)
  } else {
    p <- p + geom_text(data = lab, aes(label = label), size = 3, vjust = -0.5,
                       show.legend = FALSE)
  }
} else {
  message("No subclusters passed the thresholds; nothing to label. Try lfc_cut <- 0 or a looser fdr_cut.")
}
p

# ---- Read GFF (stop at ##FASTA if the file has sequences appended) ----
lines <- readLines(annotation_file)
fa <- grep("^##FASTA", lines)
if (length(fa)) lines <- lines[seq_len(fa[1] - 1)]

# Chromosome lengths from ##sequence-region headers
sr <- strsplit(grep("^##sequence-region", lines, value = TRUE), "\\s+")
chr_len <- do.call(rbind, lapply(sr, function(x)
  data.frame(seqid = x[2], len = as.numeric(x[4]))))

# Gene features
gff <- read.delim(text = lines[!grepl("^#", lines)], header = FALSE, sep = "\t",
                  quote = "", stringsAsFactors = FALSE,
                  col.names = c("seqid", "source", "type", "start", "end",
                                "score", "strand", "phase", "attr"))
genes <- gff[grepl("gene$", gff$type), ]
genes$GeneID <- sub(".*ID=([^;]+).*", "\\1", genes$attr)
genes$pos <- (genes$start + genes$end) / 2

# Nuclear chromosomes only, e.g. Pf3D7_01_v3 -> "Chr 1"
nuc <- grepl("^Pf3D7_\\d{2}_v3$", genes$seqid)
genes <- genes[nuc, ]
chr_len <- chr_len[grepl("^Pf3D7_\\d{2}_v3$", chr_len$seqid), ]

chr_num <- function(s) as.integer(sub("^Pf3D7_(\\d{2})_v3$", "\\1", s))
chr_len$label <- factor(paste("Chr", chr_num(chr_len$seqid)),
                        levels = rev(paste("Chr", sort(chr_num(chr_len$seqid)))))
genes$label <- factor(paste("Chr", chr_num(genes$seqid)), levels = levels(chr_len$label))

# ---- Join with edgeR gene-level results ----
gene_hits <- topSpliceDGE(ds, test = "gene", number = Inf)
gene_hits <- merge(gene_hits, genes[, c("GeneID", "pos", "label")], by = "GeneID")
gene_hits$sig <- gene_hits$FDR < 0.05
gene_hits$neglogFDR <- -log10(gene_hits$FDR)

# Sanity check: how many tested genes were matched to the GFF?
cat("Matched", nrow(gene_hits), "of", nrow(topSpliceDGE(ds, test = "gene", number = Inf)),
    "tested genes\n")

# ---- Chromosome map ----
ggplot() +
  geom_segment(data = chr_len, aes(x = 0, xend = len, y = label, yend = label),
               linewidth = 7, colour = "grey88", lineend = "round") +
  geom_point(data = gene_hits[!gene_hits$sig, ],
             aes(x = pos, y = label), shape = 124, size = 3, colour = "grey55") +
  geom_point(data = gene_hits[gene_hits$sig, ],
             aes(x = pos, y = label, colour = neglogFDR), shape = 124, size = 6) +
  scale_colour_viridis_c(option = "inferno", direction = -1, name = "-log10(FDR)") +
  scale_x_continuous(labels = function(x) paste0(x / 1e6, " Mb"),
                     expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = "Position along chromosome", y = NULL,
       title = "Genes with significant change in subcluster proportions") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())







# proportional subcluster change by gene
min_reads <- 20     # minimum gene reads per condition, to avoid noise-driven shifts

# ---- Per-gene reads in each condition ----
reads_C <- rowSums(gene_totals[, group == "C", drop = FALSE])
reads_K <- rowSums(gene_totals[, group == "K", drop = FALSE])

# ---- Per-gene proportion shift ----
gene_sum <- prop_res %>%
  dplyr::filter(!is.na(prop_C), !is.na(prop_K)) %>%
  dplyr::group_by(gene) %>%
  dplyr::filter(dplyr::n() > 1) %>%                         # need >= 2 subclusters
  dplyr::mutate(pC = prop_C / sum(prop_C),                  # renormalise so each condition sums to 1
                pK = prop_K / sum(prop_K),
                d  = pK - pC) %>%
  dplyr::summarise(
    n_subclusters  = dplyr::n(),
    shift          = 0.5 * sum(abs(d)),                     # fraction of reads that switched subcluster
    top_subcluster = subcluster[which.max(abs(d))],
    top_change     = d[which.max(abs(d))],                  # signed change of the biggest mover
    .groups = "drop") %>%
  dplyr::filter(!is.na(shift))

gene_sum$reads_C <- reads_C[gene_sum$gene]
gene_sum$reads_K <- reads_K[gene_sum$gene]
gene_sum <- gene_sum[pmin(gene_sum$reads_C, gene_sum$reads_K) >= min_reads, ]

# ---- Add edgeR gene-level significance ----
gene_sum$FDR       <- gene_hits$FDR[match(gene_sum$gene, gene_hits$GeneID)]
gene_sum$FDR_simes <- simes_hits$FDR[match(gene_sum$gene, simes_hits$GeneID)]
gene_sum$sig       <- !is.na(gene_sum$FDR) & gene_sum$FDR < 0.05
gene_sum$negLogFDR <- -log10(gene_sum$FDR)

gene_sum <- gene_sum[order(-gene_sum$shift), ]

# ---- Add the test statistic (edgeR gene-level splice test) to gene_sum ----
# gene.F = F statistic, tests whether any subcluster's usage differs between K and C
gene_sum$gene.F  <- gene_hits$gene.F[match(gene_sum$gene,  gene_hits$GeneID)]
gene_sum$P.Value <- gene_hits$P.Value[match(gene_sum$gene, gene_hits$GeneID)]
# gene_sum$FDR was already added in the previous step

gene_table <- gene_sum[, c("gene", "n_subclusters", "reads_C", "reads_K",
                           "shift", "gene.F", "P.Value", "FDR",
                           "top_subcluster", "top_change")]

# ---- Helper: stacked proportion bars, panel titles carry shift + test statistic ----
plot_gene_props <- function(genes_vec, title) {
  info <- gene_sum[match(genes_vec, gene_sum$gene), ]
  panel_lab <- paste0(info$gene,
                      "\nshift = ", sprintf("%.2f", info$shift),
                      "\nF = ", ifelse(is.na(info$gene.F), "NA", sprintf("%.2f", info$gene.F)),
                      ", FDR = ", ifelse(is.na(info$FDR), "NA", formatC(info$FDR, format = "e", digits = 1)))
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
    theme_minimal() +
    theme(strip.text = element_text(size = 8))
}

# ---- Most shifted genes ----
top_genes <- head(gene_sum$gene, 24)                 # gene_sum is already sorted by shift, descending
plot_gene_props(top_genes, "Most shifted genes: subcluster usage in C vs K")

# ---- Least shifted genes ----
# Require deeper coverage here: a gene with 20 reads and no shift is weak evidence of stability
min_reads_stable <- 100
stable <- gene_sum[pmin(gene_sum$reads_C, gene_sum$reads_K) >= min_reads_stable, ]
least_genes <- head(stable$gene[order(stable$shift)], 24)
plot_gene_props(least_genes, "Least shifted genes: subcluster usage in C vs K")

# Tables behind the plots
gene_table[gene_table$gene %in% top_genes, ]
gene_table[gene_table$gene %in% least_genes, ]







inspect_gene <- function(gid) {
  rows <- which(gene == gid)
  if (length(rows) == 0) { message(gid, " not found in the counts table"); return(invisible(NULL)) }
  
  cat("\n=====", gid, "=====\n")
  if (exists("gene_info")) {
    d <- gene_info[match(gid, gene_info$GeneID), ]
    cat("Symbol:", d$symbol, "\nDescription:", d$description, "\n")
  }
  
  # 1. Raw counts per sample
  cat("\n-- Raw counts per sample --\n")
  print(counts[rows, , drop = FALSE])
  
  # 2. Proportions and change
  cat("\n-- Subcluster proportions (mean per condition) --\n")
  print(prop_res[prop_res$gene == gid, c("subcluster", "prop_C", "prop_K", "delta_prop")],
        row.names = FALSE, digits = 3)
  
  # 3. Gene-level tests
  cat("\n-- Gene-level tests (edgeR diffSplice) --\n")
  gh <- gene_hits[gene_hits$GeneID == gid, ]
  sh <- simes_hits[simes_hits$GeneID == gid, ]
  if (nrow(gh) == 0) cat("Not tested (filtered by filterByExpr or <2 subclusters left)\n") else {
    print(gh[, c("GeneID", "NExons", "gene.F", "P.Value", "FDR")], row.names = FALSE)
    cat("Simes FDR:", signif(sh$FDR, 3), "\n")
  }
  
  # 4. Subcluster-level tests
  cat("\n-- Subcluster-level tests (vs rest of gene) --\n")
  sb <- sub_hits[sub_hits$GeneID == gid, ]
  if (nrow(sb) > 0) {
    keep_cols <- intersect(c("Subcluster", "logFC", "exon.F", "t", "P.Value", "FDR"), names(sb))
    print(sb[, keep_cols], row.names = FALSE)
  }
  # 5. Summary shift, if computed
  if (exists("gene_sum") && gid %in% gene_sum$gene) {
    cat("\n-- Proportion shift --\n")
    print(gene_sum[gene_sum$gene == gid, c("shift", "top_subcluster", "top_change", "reads_C", "reads_K")],
          row.names = FALSE)
  }
  
  # 6. Plots
  if (gid %in% ds$genes$GeneID) plotSpliceDGE(ds, geneid = gid, genecol = "GeneID")
  if (exists("plot_gene_props") && gid %in% gene_sum$gene) print(plot_gene_props(gid, gid))
  
  invisible(NULL)
}

inspect_gene("PF3D7_1216200")
