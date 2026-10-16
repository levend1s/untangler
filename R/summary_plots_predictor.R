library(dplyr)
library(edgeR)
library(pheatmap)
library(ggplot2)

# ==========================================================================
# INPUT PATHS
# ==========================================================================

# the input here should be a cluster summary with some predictor information
annotation_file <- "~/Downloads/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff"

cluster_summary <- "/Users/jlevendis/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/untangler_test/cluster_transcripts_28hpi/cluster_summary.tsv"
counts_table    <- "/Users/jlevendis/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/untangler_test/cluster_transcripts_28hpi/counts_table.tsv"

cluster_summary <- "/Users/jlevendis/untangler/cluster_transcripts_PF3D7_1249300_28_90p/cluster_summary.tsv"
counts_table <- "/Users/jlevendis/untangler/cluster_transcripts_PF3D7_1249300_28_90p/counts_table.tsv"

# cluster_summary <- "/Users/jlevendis/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/untangler_test/cluster_transcripts_32hpi/cluster_summary.tsv"
# counts_table    <- "/Users/jlevendis/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/untangler_test/cluster_transcripts_32hpi/counts_table.tsv"

# cluster_summary <- "/Users/jlevendis/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/untangler_test/cluster_transcripts_36hpi/cluster_summary.tsv"
# counts_table    <- "/Users/jlevendis/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/untangler_test/cluster_transcripts_36hpi/counts_table.tsv"

# ==========================================================================
# PLOT FUNCTIONS
# (defined up front; called at the very end of the script once all data
#  prep, normalisation, and testing is done)
# ==========================================================================

#' Cluster abundance heatmap (logCPM, row-scaled)
#'
#' @param mat            matrix of logCPM values (rows = clusters, cols = samples)
#' @param row_labels     character vector of row labels (same order as rows of mat)
#' @param ann_col        data.frame of column annotations (condition, optionally hpi)
#' @param gaps           numeric vector of column gap positions, or NULL
#' @param title          plot title
plot_cluster_heatmap <- function(mat, row_labels, ann_col, gaps, title) {
  pheatmap(
    mat,
    scale = "row",
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    labels_row = row_labels,
    annotation_col = ann_col,
    gaps_col = gaps,
    color = colorRampPalette(c("navy", "white", "firebrick"))(100),
    show_rownames = nrow(mat) <= 60,
    fontsize_row = 7,
    main = title
  )
}

#' Volcano plot of cluster differential abundance
#'
#' @param vol             data.frame with logFC, negLogFDR, in_new_df, sig, label
#' @param contrast_label  string describing the contrast, used in the plot title
#' @param lab             data.frame subset of vol to label (top significant clusters)
plot_volcano <- function(vol, contrast_label, lab) {
  p <- ggplot(vol, aes(x = logFC, y = negLogFDR)) +
    geom_point(data = vol[!vol$in_new_df, ], colour = "grey70", size = 1, alpha = 0.6) +
    geom_point(data = vol[vol$in_new_df, ],  colour = "firebrick", size = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    labs(x = "log2 fold change in cluster abundance",
         y = "-log10(FDR)",
         title = paste("Cluster abundance:", contrast_label),
         subtitle = "Red = clusters in new_df; grey = all other clusters") +
    theme_classic()
  
  if (nrow(lab) > 0) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(data = lab, aes(label = label), size = 3,
                                        max.overlaps = Inf)
    } else {
      p <- p + geom_text(data = lab, aes(label = label), size = 3, vjust = -0.5)
    }
  }
  p
}

#' MD plot (mean abundance vs log fold change) of cluster differential abundance
#'
#' @param vol             data.frame with logCPM, logFC, in_new_df, label
#' @param contrast_label  string describing the contrast, used in the plot title
#' @param lab             data.frame subset of vol to label (top significant clusters,
#'                         same set used in plot_volcano for consistency)
plot_md <- function(vol, contrast_label, lab) {
  p_md <- ggplot(vol, aes(x = logCPM, y = logFC)) +
    geom_point(data = vol[!vol$in_new_df, ], colour = "grey70", size = 1, alpha = 0.6) +
    geom_point(data = vol[vol$in_new_df, ],  colour = "firebrick", size = 2) +
    geom_hline(yintercept = 0, colour = "black") +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed") +
    labs(x = "Average abundance (log2 CPM)",
         y = "log2 fold change in cluster abundance",
         title = paste("MD plot:", contrast_label),
         subtitle = "Red = clusters in new_df; grey = all other clusters") +
    theme_classic()
  
  if (nrow(lab) > 0) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p_md <- p_md + ggrepel::geom_text_repel(data = lab, aes(label = label), size = 3,
                                              max.overlaps = Inf)
    } else {
      p_md <- p_md + geom_text(data = lab, aes(label = label), size = 3, vjust = -0.5)
    }
  }
  p_md
}


# ==========================================================================
# HELPERS
# ==========================================================================

get_attr <- function(attr, key) {
  out <- rep(NA_character_, length(attr))
  hit <- grepl(paste0("(^|;)", key, "="), attr)
  out[hit] <- sub(paste0("^.*(^|;)", key, "=([^;]*).*$"), "\\2", attr[hit])
  out
}

lookup <- function(ids, col) gene_info[[col]][match(ids, gene_info$GeneID)]

clean <- function(x, n = 35) {
  x[grepl("unknown function", x, ignore.case = TRUE)] <- "unannotated"
  x[is.na(x) | x == ""] <- "none"
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 3), "..."), x)
}


# ==========================================================================
# LOAD CLUSTER SUMMARY
# ==========================================================================

df <- read.delim(cluster_summary, sep = "\t", header = TRUE, quote = "",
                 fill = TRUE, na.strings = c("", "NA"), stringsAsFactors = FALSE)

res <- df[!is.na(df$predictor_combo) & trimws(df$predictor_combo) != "[]", ]

res$gene <- sub("_cluster.*$", "", res$index)

new_df <- res %>%
  dplyr::select(gene, cluster, predictor_combo, target)


# ==========================================================================
# READ GFF (stop at ##FASTA if sequences are appended)
# ==========================================================================

lines <- readLines(annotation_file)
fa <- grep("^##FASTA", lines)
if (length(fa)) lines <- lines[seq_len(fa[1] - 1)]

gff <- read.delim(text = lines[!grepl("^#", lines)], header = FALSE, sep = "\t",
                  quote = "", stringsAsFactors = FALSE,
                  col.names = c("seqid", "source", "type", "start", "end",
                                "score", "strand", "phase", "attr"))
genes <- gff[grepl("gene$", gff$type), ]
genes$GeneID <- sub(".*ID=([^;]+).*", "\\1", genes$attr)

genes$description <- vapply(unname(get_attr(genes$attr, "description")),
                            function(x) if (is.na(x)) NA_character_ else URLdecode(x),
                            character(1))
genes$symbol <- get_attr(genes$attr, "Name")

gene_info <- genes[, c("GeneID", "symbol", "description")]
gene_info <- gene_info[!duplicated(gene_info$GeneID), ]


# ==========================================================================
# FIND LEFT/RIGHT NEIGHBOURS ALONG EACH CHROMOSOME
# ==========================================================================

genes_o <- genes[order(genes$seqid, genes$start), ]
genes_o$left_id  <- ave(genes_o$GeneID, genes_o$seqid, FUN = function(x) c(NA, x[-length(x)]))
genes_o$right_id <- ave(genes_o$GeneID, genes_o$seqid, FUN = function(x) c(x[-1], NA))

# 5' / 3' (upstream / downstream) depends on strand
genes_o$id_5p <- ifelse(genes_o$strand == "-", genes_o$right_id, genes_o$left_id)
genes_o$id_3p <- ifelse(genes_o$strand == "-", genes_o$left_id,  genes_o$right_id)

neighbours <- data.frame(
  GeneID  = genes_o$GeneID,
  id_5p   = genes_o$id_5p,
  name_5p = lookup(genes_o$id_5p, "symbol"),
  desc_5p = lookup(genes_o$id_5p, "description"),
  id_3p   = genes_o$id_3p,
  name_3p = lookup(genes_o$id_3p, "symbol"),
  desc_3p = lookup(genes_o$id_3p, "description"),
  stringsAsFactors = FALSE
)
neighbours <- neighbours[!duplicated(neighbours$GeneID), ]


# ==========================================================================
# ADD GENE INFO AND NEIGHBOURS TO new_df (keyed on `gene`)
# ==========================================================================

drop_cols <- c("symbol", "description", "id_5p", "name_5p", "desc_5p",
               "id_3p", "name_3p", "desc_3p")
new_df <- new_df[, !(names(new_df) %in% drop_cols)]

i_info <- match(new_df$gene, gene_info$GeneID)
new_df$symbol      <- gene_info$symbol[i_info]
new_df$description <- gene_info$description[i_info]

i_nb <- match(new_df$gene, neighbours$GeneID)
for (col in c("id_5p", "name_5p", "desc_5p", "id_3p", "name_3p", "desc_3p")) {
  new_df[[col]] <- neighbours[[col]][i_nb]
}

new_df <- new_df %>%
  dplyr::select(gene, cluster, predictor_combo, target,
                symbol, description,
                id_5p, name_5p, desc_5p,
                id_3p, name_3p, desc_3p)


# ==========================================================================
# LOAD COUNTS AND PARSE SAMPLE NAMES (e.g. 36K1 -> 36 hpi, K, rep 1)
# ==========================================================================

counts <- read.delim(counts_table, row.names = 1, check.names = FALSE)
colnames(counts) <- sub("_read_depth$", "", colnames(counts))   # in case of old-style names

meta <- data.frame(sample = colnames(counts), stringsAsFactors = FALSE)
meta$hpi       <- as.integer(sub("^(\\d+)([A-Za-z])(\\d+)$", "\\1", meta$sample))
meta$condition <- sub("^(\\d+)([A-Za-z])(\\d+)$", "\\2", meta$sample)
meta$rep       <- as.integer(sub("^(\\d+)([A-Za-z])(\\d+)$", "\\3", meta$sample))

meta <- meta[order(meta$hpi, meta$condition, meta$rep), ]
counts <- counts[, meta$sample]
meta$hpi       <- factor(meta$hpi)
meta$condition <- factor(meta$condition, levels = c("C", "K"))
rownames(meta) <- meta$sample

tps    <- levels(meta$hpi)
n_tp   <- length(tps)
multi  <- n_tp > 1
cat("Timepoints detected:", paste(tps, collapse = ", "), "hpi\n")


# ==========================================================================
# MATCH new_df CLUSTERS TO THE COUNTS TABLE
# ==========================================================================

new_df$ID_cluster <- paste0(new_df$gene, "_", new_df$cluster)
cat("Matched", sum(new_df$ID_cluster %in% rownames(counts)), "of", nrow(new_df), "clusters\n")

# ---- Normalise using ALL clusters ----
y_all  <- calcNormFactors(DGEList(counts = counts))
logcpm <- cpm(y_all, log = TRUE, prior.count = 1)

# ---- Log2 fold changes of mean logCPM ----
grp <- paste0(meta$hpi, meta$condition)
mean_lcpm <- sapply(unique(grp), function(g) rowMeans(logcpm[, grp == g, drop = FALSE]))

lfc <- data.frame(ID_cluster = rownames(mean_lcpm))
if (multi) {
  # each later timepoint vs the first, within each condition
  for (cond in levels(meta$condition)) {
    for (tp in tps[-1]) {
      a <- paste0(tp, cond); b <- paste0(tps[1], cond)
      if (all(c(a, b) %in% colnames(mean_lcpm)))
        lfc[[paste0(cond, "_", tp, "v", tps[1])]] <- mean_lcpm[, a] - mean_lcpm[, b]
    }
  }
} else {
  # single timepoint: K vs C
  a <- paste0(tps, "K"); b <- paste0(tps, "C")
  if (all(c(a, b) %in% colnames(mean_lcpm)))
    lfc[[paste0("K_vs_C_", tps)]] <- mean_lcpm[, a] - mean_lcpm[, b]
}


# ==========================================================================
# edgeR TEST
# ==========================================================================

if (multi) {
  design <- model.matrix(~ condition + hpi, data = meta)
  test_coef <- grep("^hpi", colnames(design))          # joint test across timepoints
  cat("Testing: change across timepoints (adjusting for condition)\n")
} else {
  design <- model.matrix(~ condition, data = meta)
  test_coef <- 2                                       # K vs C
  cat("Testing: K vs C at", tps, "hpi\n")
}

keep <- filterByExpr(y_all, design)
y <- calcNormFactors(y_all[keep, , keep.lib.sizes = FALSE])
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = test_coef)

test_res <- topTags(qlf, n = Inf)$table
test_res$ID_cluster <- rownames(test_res)

# ---- Add results to new_df ----
new_df <- new_df[, !(names(new_df) %in% c(setdiff(names(lfc), "ID_cluster"), "test_FDR"))]
new_df <- merge(new_df, lfc, by = "ID_cluster", all.x = TRUE, sort = FALSE)
new_df$test_FDR <- test_res$FDR[match(new_df$ID_cluster, test_res$ID_cluster)]


# ==========================================================================
# PREP HEATMAP INPUTS
# ==========================================================================

ids <- intersect(new_df$ID_cluster, rownames(counts))
ids <- ids[rowSums(counts[ids, ]) >= 10]
heatmap_mat <- logcpm[ids, , drop = FALSE]

# Labels: "cluster1 | gene id | description"
i <- match(ids, new_df$ID_cluster)
cluster_short <- sub("^(cluster\\d+).*$", "\\1", new_df$cluster[i])
heatmap_row_labels <- paste0(cluster_short, " | ", new_df$gene[i], " | ",
                             clean(new_df$description[i]))

# Only annotate/gap by timepoint when there is more than one
heatmap_ann_col <- data.frame(condition = meta$condition, row.names = meta$sample)
if (multi) heatmap_ann_col$hpi <- meta$hpi
heatmap_gaps <- if (multi) head(cumsum(table(meta$hpi)), -1) else NULL

heatmap_title <- paste0("Cluster abundance, ", paste(tps, collapse = "/"),
                        " hpi (log2 CPM)")


# ==========================================================================
# PREP VOLCANO / MD PLOT INPUTS
# ==========================================================================

# ---- Choose the contrast for the volcano ----
# "K_vs_C": treatment effect (adjusted for timepoint if there are several)
# "time"  : last timepoint vs first (only meaningful with >1 timepoint)
volcano_contrast <- "K_vs_C"

if (volcano_contrast == "time" && multi) {
  coef_use   <- tail(grep("^hpi", colnames(design)), 1)
  contrast_label <- paste0(tail(tps, 1), " vs ", tps[1], " hpi")
} else {
  coef_use   <- grep("^condition", colnames(design))    # conditionK
  contrast_label <- if (multi) "K vs C (adjusted for timepoint)" else paste0("K vs C at ", tps, " hpi")
}

# ---- Differential abundance across ALL clusters ----
vol <- topTags(glmQLFTest(fit, coef = coef_use), n = Inf)$table
vol$ID_cluster <- rownames(vol)
vol$in_new_df  <- vol$ID_cluster %in% new_df$ID_cluster
vol$sig        <- vol$FDR < 0.05 & abs(vol$logFC) >= 1
vol$negLogFDR  <- -log10(vol$FDR)

vol$label <- paste0(sub("_cluster.*$", "", vol$ID_cluster), " ",
                    sub("^.*_(cluster\\d+).*$", "\\1", vol$ID_cluster))

cat(sum(vol$in_new_df), "of", nrow(new_df), "new_df clusters were tested (rest filtered out as low count)\n")
cat(sum(vol$sig), "clusters significant overall;", sum(vol$sig & vol$in_new_df), "of them in new_df\n")

# Top significant new_df clusters, used to label both the volcano and MD plots
volcano_lab <- head(vol[vol$in_new_df & vol$sig, ][order(vol$FDR[vol$in_new_df & vol$sig]), ], 15)

# ==========================================================================
# PLOT
# ==========================================================================

plot_cluster_heatmap(
  mat = heatmap_mat,
  row_labels = heatmap_row_labels,
  ann_col = heatmap_ann_col,
  gaps = heatmap_gaps,
  title = heatmap_title
)

p_volcano <- plot_volcano(vol, contrast_label, volcano_lab)
print(p_volcano)

p_md <- plot_md(vol, contrast_label, volcano_lab)
print(p_md)


# ==========================================================================
# SUMMARY TABLES
# ==========================================================================

# Significant clusters
new_df %>% dplyr::filter(test_FDR < 0.05) %>%
  dplyr::select(gene, cluster, dplyr::starts_with(c("C_", "K_")), test_FDR)

# Table of significant new_df clusters (volcano contrast)
vol_sig <- merge(vol[vol$in_new_df & vol$sig, c("ID_cluster", "logFC", "FDR")],
                 new_df[, c("ID_cluster", "gene", "cluster", "description")],
                 by = "ID_cluster")
vol_sig[order(vol_sig$FDR), ]