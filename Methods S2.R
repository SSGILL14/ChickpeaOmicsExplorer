##############################################################################
# Methods S2: WGCNA of chickpea tissue expression
##############################################################################

# ── 1. Install / load packages ──────────────────────────────────────────────
required_pkgs <- c("WGCNA", "readxl", "flashClust", "dynamicTreeCut",
                   "ggplot2", "reshape2", "pheatmap", "RColorBrewer",
                   "gridExtra", "plotly", "htmlwidgets")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "WGCNA") {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("WGCNA")
    } else {
      install.packages(pkg)
    }
  }
}

library(WGCNA)
library(readxl)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(plotly)
library(htmlwidgets)

allowWGCNAThreads()

# ── 2. Read and prepare data ────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
preflight_only <- "--preflight" %in% args
positional_args <- args[!grepl("^--", args)]

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE))
} else {
  getwd()
}
input_candidates <- c(file.path(script_dir, "Data S1.xlsx"),
                      file.path(getwd(), "Supplementary Files", "Data S1.xlsx"),
                      file.path(getwd(), "Data S1.xlsx"),
                      "Data S1.xlsx")
default_input <- input_candidates[file.exists(input_candidates)][1]
if (is.na(default_input)) default_input <- input_candidates[1]

file_path <- if (length(positional_args) >= 1) positional_args[1] else default_input
output_dir <- if (length(positional_args) >= 2) positional_args[2] else
  file.path(dirname(normalizePath(file_path, mustWork = FALSE)), "WGCNA outputs")

if (!file.exists(file_path)) stop("Input workbook not found: ", file_path)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- function(name) file.path(output_dir, name)

# Data S1 has a descriptive heading in row 1 and column names in row 2.
dat <- read_excel(file_path, sheet = "Sheet1", skip = 1)

tissue_cols <- c("GS","S","R","Rtip","RH","YL","ML","Brac","SAM",
                 "FB1","FB2","FB3","FB4","FL1","FL2","FL3","FL4","FL5",
                 "Cal","Cor","And","Gyn","Pedi","PodSh","SdCt",
                 "Emb","Endo","5DAP","10DAP","20DAP","30DAP","Nod")

required_cols <- c("Transcript id", "mRNA", "lncRNA", tissue_cols)
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0)
  stop("Missing required Data S1 columns: ", paste(missing_cols, collapse = ", "))

dat$`Transcript id` <- trimws(as.character(dat$`Transcript id`))
if (any(is.na(dat$`Transcript id`) | dat$`Transcript id` == ""))
  stop("Blank transcript identifiers were found in Data S1.")

expr <- as.data.frame(dat[, tissue_cols])
expr <- as.data.frame(lapply(expr, function(x) suppressWarnings(as.numeric(x))),
                      check.names = FALSE)
if (anyNA(expr)) stop("Missing or non-numeric expression values were found in Data S1.")
if (any(as.matrix(expr) < 0)) stop("Negative FPKM values were found in Data S1.")

# Some transcripts occur more than once in Data S1. They must agree exactly
# before one copy is retained, otherwise the analysis stops.
duplicate_ids <- unique(dat$`Transcript id`[duplicated(dat$`Transcript id`)])
if (length(duplicate_ids) > 0) {
  discordant <- vapply(duplicate_ids, function(id) {
    rows <- which(dat$`Transcript id` == id)
    length(unique(apply(expr[rows, , drop = FALSE], 1, paste, collapse = "\r"))) > 1
  }, logical(1))
  if (any(discordant))
    stop("Duplicate transcript IDs with discordant expression values: ",
         paste(duplicate_ids[discordant], collapse = ", "))
  discordant_annot <- vapply(duplicate_ids, function(id) {
    rows <- which(dat$`Transcript id` == id)
    length(unique(paste(dat$mRNA[rows], dat$lncRNA[rows], sep = "\r"))) > 1
  }, logical(1))
  if (any(discordant_annot))
    stop("Duplicate transcript IDs with discordant mRNA/lncRNA annotation: ",
         paste(duplicate_ids[discordant_annot], collapse = ", "))
}

keep_unique <- !duplicated(dat$`Transcript id`)
expr <- expr[keep_unique, , drop = FALSE]
dat_unique <- dat[keep_unique, , drop = FALSE]
rownames(expr) <- dat_unique$`Transcript id`

orig_annot <- data.frame(
  Transcript = dat_unique$`Transcript id`,
  mRNA = dat_unique$mRNA,
  lncRNA = dat_unique$lncRNA,
  stringsAsFactors = FALSE
)

cat("Raw matrix:", nrow(expr), "transcripts x", ncol(expr), "tissues\n")
cat("Input:", normalizePath(file_path), "\n")
cat("Duplicate rows collapsed after concordance check:",
    nrow(dat) - nrow(dat_unique), "\n")

# ── 3. Filter genes (paper's criteria) ──────────────────────────────────────
# Step 1: Keep genes with FPKM >= 0.1 in at least one tissue
keep_expressed <- apply(expr, 1, function(x) max(x) >= 0.1)
expr_filt <- expr[keep_expressed, ]
cat("After FPKM >= 0.1 filter:", nrow(expr_filt), "transcripts\n")

# Step 2: Log2 transform
expr_log <- log2(expr_filt + 1)

# Step 3: Variance > 1.5 filter on log2 values
gene_var <- apply(expr_log, 1, var)
keep_var <- gene_var > 1.5
expr_final <- expr_log[keep_var, ]
cat("After variance > 1.5 filter:", nrow(expr_final), "transcripts\n")
cat("(Paper reported 13,176 transcripts at this step)\n\n")

filter_summary <- data.frame(
  Stage = c("Unique transcript profiles", "FPKM >= 0.1 in >=1 tissue",
            "log2(FPKM + 1) variance > 1.5"),
  Transcripts = c(nrow(expr), nrow(expr_filt), nrow(expr_final))
)
write.csv(filter_summary, out_file("WGCNA_filter_summary.csv"), row.names = FALSE)
write.csv(data.frame(Transcript = rownames(expr_final)),
          out_file("WGCNA_retained_transcripts.csv"), row.names = FALSE)

if (preflight_only) {
  cat("Preflight completed; network construction was not run.\n")
  quit(save = "no", status = 0)
}

# ── 4. Transpose for WGCNA ──────────────────────────────────────────────────
datExpr <- as.data.frame(t(expr_final))
geneNames <- colnames(datExpr)

gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0)
    cat("Removing", sum(!gsg$goodGenes), "problematic genes\n")
  if (sum(!gsg$goodSamples) > 0)
    cat("Removing", sum(!gsg$goodSamples), "problematic samples\n")
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  geneNames <- colnames(datExpr)
}

cat("Final WGCNA input:", nrow(datExpr), "samples x", length(geneNames), "genes\n\n")

# ── 5. Sample clustering ────────────────────────────────────────────────────
sampleTree <- hclust(dist(datExpr), method = "average")

pdf(out_file("WGCNA_01_sample_clustering.pdf"), width = 14, height = 6)
par(mar = c(4, 5, 2, 1))
plot(sampleTree, main = "Sample Clustering to Detect Outliers",
     sub = "", xlab = "", cex = 0.8)
dev.off()
cat("Sample dendrogram saved.\n")

# ── 6. Pick soft-thresholding power ─────────────────────────────────────────
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers,
                         networkType = "unsigned", verbose = 3)

pdf(out_file("WGCNA_02_soft_threshold.pdf"), width = 12, height = 5)
par(mfrow = c(1, 2))

plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (signed R-sq)",
     type = "n", main = "Scale Independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red", cex = 0.9)
abline(h = 0.86, col = "blue", lty = 2)
abline(h = 0.80, col = "green", lty = 2)
legend("bottomright", legend = c("R2=0.86 (paper)", "R2=0.80"),
       col = c("blue","green"), lty = 2, cex = 0.7)

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n", main = "Mean Connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5],
     labels = powers, col = "red", cex = 0.9)

dev.off()
cat("Soft threshold plot saved.\n")

cat("\n-- Soft Threshold Results --\n")
print(sft$fitIndices[, c(1, 2, 5)])

softPower <- 9
write.csv(sft$fitIndices, out_file("WGCNA_soft_threshold_diagnostics.csv"),
          row.names = FALSE)
power_row <- match(softPower, sft$fitIndices[, 1])
if (is.na(power_row)) stop("Selected soft-thresholding power was not evaluated.")
observed_fit <- -sign(sft$fitIndices[power_row, 3]) * sft$fitIndices[power_row, 2]
cat("\nUsing soft-thresholding power beta =", softPower, "\n")
cat("Observed signed R2 at this power:", round(observed_fit, 3), "\n\n")

# ── 7. Build adjacency and TOM ──────────────────────────────────────────────
cat("Calculating adjacency matrix (unsigned)... this takes ~10-20 min\n")
adjacency <- adjacency(datExpr, power = softPower, type = "unsigned")

cat("Calculating TOM (unsigned)...\n")
TOM <- TOMsimilarity(adjacency, TOMType = "unsigned")
dissTOM <- 1 - TOM

rm(adjacency)
gc()

# ── 8. Hierarchical clustering and module detection ─────────────────────────
cat("Clustering genes...\n")
geneTree <- hclust(as.dist(dissTOM), method = "average")

dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = 30)

dynamicColors <- labels2colors(dynamicMods)
names(dynamicColors) <- geneNames

nPreMerge <- sum(unique(dynamicMods) != 0)
cat("Dynamic tree cut found", nPreMerge, "modules (before merging)\n")
cat("Module sizes:\n")
print(table(dynamicMods))

pdf(out_file("WGCNA_03_gene_dendrogram_pre_merge.pdf"), width = 14, height = 6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Modules",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene Dendrogram and Module Colors (before merging)")
dev.off()
cat("Pre-merge dendrogram saved.\n")

# ── 9. Merge similar modules ────────────────────────────────────────────────
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes

MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

pdf(out_file("WGCNA_04_module_eigengene_clustering.pdf"), width = 10, height = 6)
plot(METree, main = "Clustering of Module Eigengenes",
     xlab = "", sub = "", cex = 0.8)
abline(h = 0.2, col = "red", lty = 2)
text(1, 0.22, "Merge cutHeight = 0.2", col = "red", cex = 0.8)
dev.off()

mergedMods <- mergeCloseModules(datExpr, dynamicColors,
                                cutHeight = 0.2, verbose = 3)
mergedColors <- mergedMods$colors
mergedMEs <- mergedMods$newMEs
names(mergedColors) <- geneNames

nModules <- length(unique(mergedColors)) - ifelse("grey" %in% mergedColors, 1, 0)
cat("\nAfter merging:", nModules, "assigned modules\n")

pdf(out_file("WGCNA_05_gene_dendrogram_merged.pdf"), width = 14, height = 7)
plotDendroAndColors(geneTree,
                    cbind(dynamicColors, mergedColors),
                    c("Before Merge", "After Merge"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene Dendrogram with Module Colors")
dev.off()
cat("Merged dendrogram saved.\n\n")

# ── 10. Rename modules M1, M2, ... (sorted by size, descending) ────────────
moduleTable <- sort(table(mergedColors), decreasing = TRUE)
moduleTable <- moduleTable[names(moduleTable) != "grey"]

colorToModule <- setNames(
  paste0("M", seq_along(moduleTable)),
  names(moduleTable)
)

moduleLabels <- ifelse(mergedColors == "grey", "Unassigned",
                       colorToModule[mergedColors])
names(moduleLabels) <- geneNames

module_key <- data.frame(
  ModuleColor = names(colorToModule),
  ModuleName = unname(colorToModule),
  Transcripts = as.integer(moduleTable),
  stringsAsFactors = FALSE
)
write.csv(module_key, out_file("WGCNA_module_key.csv"), row.names = FALSE)

cat("-- Module Sizes (M1 = largest) --\n")
for (i in seq_along(moduleTable)) {
  cat(sprintf("  %s (%s): %d transcripts\n",
              colorToModule[names(moduleTable)[i]],
              names(moduleTable)[i],
              moduleTable[i]))
}
cat(sprintf("  Unassigned (grey): %d transcripts\n\n", sum(mergedColors == "grey")))

# ── 11. Module-Tissue correlation heatmap ───────────────────────────────────
MEs_ordered <- orderMEs(mergedMEs)

# Rename columns to M1, M2, ...
me_colors <- gsub("^ME", "", colnames(MEs_ordered))
new_names <- ifelse(me_colors == "grey", "Unassigned",
                    colorToModule[me_colors])
colnames(MEs_ordered) <- new_names

nSamples <- nrow(datExpr)

# Tissue identity matrix
tissues <- rownames(datExpr)
traitMatrix <- matrix(0, nrow = nSamples, ncol = nSamples)
diag(traitMatrix) <- 1
rownames(traitMatrix) <- tissues
colnames(traitMatrix) <- tissues

# Correlation and p-values
moduleTraitCor <- cor(MEs_ordered, traitMatrix, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
moduleTraitFDR <- matrix(NA_real_, nrow = nrow(moduleTraitPvalue),
                         ncol = ncol(moduleTraitPvalue),
                         dimnames = dimnames(moduleTraitPvalue))
assigned_rows <- rownames(moduleTraitPvalue) != "Unassigned"
moduleTraitFDR[assigned_rows, ] <- matrix(
  p.adjust(as.vector(moduleTraitPvalue[assigned_rows, , drop = FALSE]),
           method = "BH"),
  nrow = sum(assigned_rows),
  dimnames = list(rownames(moduleTraitPvalue)[assigned_rows],
                  colnames(moduleTraitPvalue))
)

# ── Biological tissue order (NO grouping — tissues in natural order) ─────────
tissue_order <- c("GS","S","R","Rtip","RH","YL","ML","Brac","SAM",
                  "FB1","FB2","FB3","FB4","FL1","FL2","FL3","FL4","FL5",
                  "Cal","Cor","And","Gyn","Pedi","PodSh","SdCt",
                  "Emb","Endo","5DAP","10DAP","20DAP","30DAP","Nod")

# Order modules M1 → Mn numerically (M1 at top for publication convention)
all_mnames <- unique(colorToModule)
all_mnums  <- as.numeric(gsub("M", "", all_mnames))
module_order <- all_mnames[order(all_mnums)]   # M1, M2, M3 ... Mn top-to-bottom
if ("Unassigned" %in% rownames(moduleTraitCor)) {
  module_order <- c(module_order, "Unassigned")
}
module_order <- module_order[module_order %in% rownames(moduleTraitCor)]

# Subset and order the matrices
cor_ordered  <- moduleTraitCor[module_order,  tissue_order]
pval_ordered <- moduleTraitPvalue[module_order, tissue_order]
fdr_ordered  <- moduleTraitFDR[module_order, tissue_order]

write.csv(cor_ordered, out_file("WGCNA_module_tissue_correlations.csv"))
write.csv(pval_ordered, out_file("WGCNA_module_tissue_raw_pvalues.csv"))
write.csv(fdr_ordered, out_file("WGCNA_module_tissue_BH_FDR.csv"))

# ──────────────────────────────────────────────────────────────────────────────
# PUBLICATION HEATMAP — Section A: WGCNA labeledHeatmap style
#   • Both r and p displayed inside every cell
#   • blue-white-red colour scale (-1 → 0 → +1)
# ──────────────────────────────────────────────────────────────────────────────

# Cell labels: correlation on line 1, BH-adjusted p-value in parentheses on line 2.
# Four decimal places prevent high adjusted P values such as 0.9853 from being
# displayed as 1, while very small values remain legible in scientific notation.
# The unassigned row is descriptive and therefore displays n/a.
format_adjusted_p <- function(x) {
  out <- rep("n/a", length(x))
  keep <- !is.na(x)
  out[keep] <- ifelse(
    x[keep] < 0.0001,
    formatC(x[keep], format = "e", digits = 1),
    formatC(x[keep], format = "f", digits = 4)
  )
  dim(out) <- dim(x)
  dimnames(out) <- dimnames(x)
  out
}
fdr_display <- format_adjusted_p(fdr_ordered)
textMatrix <- paste0(
  signif(cor_ordered,  2), "\n(",
  fdr_display, ")"
)
dim(textMatrix) <- dim(cor_ordered)

# Row labels include gene count
row_labels_full <- paste0(module_order, " (n=",
                          sapply(module_order, function(m) {
                            if (m == "Unassigned") sum(mergedColors == "grey")
                            else sum(moduleLabels == m)
                          }), ")")

pdf(out_file("WGCNA_06_module_tissue_heatmap_publication.pdf"),
    width = 22, height = max(10, nrow(cor_ordered) * 0.45 + 3))

par(mar = c(9, 12, 4, 3))   # bottom, left, top, right

labeledHeatmap(
  Matrix        = cor_ordered,
  xLabels       = tissue_order,
  yLabels       = row_labels_full,
  ySymbols      = row_labels_full,   # suppress colour squares on y-axis
  colorLabels   = FALSE,
  colors        = blueWhiteRed(100),
  textMatrix    = textMatrix,
  setStdMargins = FALSE,
  cex.text      = 0.40,              # cell text size
  cex.lab       = 0.75,              # axis label size
  zlim          = c(-1, 1),
  main          = paste0(
    "Module-Tissue Correlation\n",
    "Upper value: Pearson r  |  Lower value: BH-adjusted p  |  Colour: r scale (-1 to +1)"
  )
)

dev.off()
cat("Publication labeledHeatmap saved: WGCNA_06_module_tissue_heatmap_publication.pdf\n")

# ──────────────────────────────────────────────────────────────────────────────
# PUBLICATION HEATMAP — Section B: heatmap style (journal-ready)
#   • Significance stars use BH-adjusted p-values
#   • r value shown in each cell
#   • NO column/row clustering  |  NO annotation sidebars
#   • White grid lines between cells for clarity
# ──────────────────────────────────────────────────────────────────────────────

# Build cell annotation: "r\n*" where stars encode significance
sig_stars <- ifelse(is.na(fdr_ordered), "",
             ifelse(fdr_ordered < 0.001, "***",
             ifelse(fdr_ordered < 0.01,  "**",
             ifelse(fdr_ordered < 0.05,  "*", ""))))

# Format: r value + significance star (two lines)
cell_labels <- matrix(
  paste0(sprintf("%.2f", cor_ordered), "\n", sig_stars),
  nrow = nrow(cor_ordered),
  dimnames = list(module_order, tissue_order)
)

# Colour palette: deep blue → white → deep red (publication-grade)
pub_colors <- colorRampPalette(
  c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0",
    "white",
    "#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F")
)(101)

# Dynamic figure dimensions
n_rows <- nrow(cor_ordered)
n_cols <- ncol(cor_ordered)
fig_w  <- max(16, n_cols * 0.55 + 4)
fig_h  <- max(10, n_rows * 0.45 + 3)

pdf(out_file("WGCNA_06b_module_tissue_heatmap_pheatmap_publication.pdf"),
    width = fig_w, height = fig_h)

hm <- pheatmap(
  mat              = cor_ordered,
  color            = pub_colors,
  breaks           = seq(-1, 1, length.out = 102),

  # ── No clustering — tissues shown in biological order ──────────────────
  cluster_rows     = FALSE,
  cluster_cols     = FALSE,

  # ── No annotation sidebars (no GT/RT/FDS/FP/ST grouping) ───────────────
  annotation_col   = NA,
  annotation_row   = NA,

  # ── Labels ──────────────────────────────────────────────────────────────
  labels_row       = row_labels_full,
  labels_col       = tissue_order,

  # ── Cell numbers: r value + significance stars ───────────────────────
  display_numbers  = cell_labels,
  number_color     = "black",
  fontsize_number  = 6.5,

  # ── Typography ──────────────────────────────────────────────────────────
  fontsize_row     = 8,
  fontsize_col     = 8,
  angle_col        = 45,

  # ── Grid / borders ──────────────────────────────────────────────────────
  border_color     = "white",    # white cell borders → clean grid look

  # ── Legend ──────────────────────────────────────────────────────────────
  legend           = TRUE,
  legend_breaks    = seq(-1, 1, by = 0.25),
  legend_labels    = as.character(seq(-1, 1, by = 0.25)),

  # ── Title ───────────────────────────────────────────────────────────────
  main             = paste0(
    "Module-Tissue Correlation (Pearson r)\n",
    "Cell: r value | Stars: BH-adjusted p (*** <0.001, ** <0.01, * <0.05)"
  ),

  # ── Silent rendering ─────────────────────────────────────────────────────
  silent           = TRUE
)
grid::grid.newpage()
grid::grid.draw(hm$gtable)

dev.off()
cat("Publication pheatmap saved: WGCNA_06b_module_tissue_heatmap_pheatmap_publication.pdf\n")

# ──────────────────────────────────────────────────────────────────────────────
# PUBLICATION HEATMAP — Section C: ggplot2 version (maximum customisability)
#   • WGCNA blueWhiteRed colour scheme (matching _06 labeledHeatmap)
#   • Two values per cell: r on line 1, (p-value) on line 2  — matching _06
#   • Bold, large axis labels for print/publication visibility
#   • No grey background, clean white grid borders
# ──────────────────────────────────────────────────────────────────────────────

# Melt correlation and BH-adjusted p-value matrices to long format
cor_long  <- reshape2::melt(cor_ordered,  varnames = c("Module", "Tissue"),
                             value.name = "r")
pval_long <- reshape2::melt(fdr_ordered, varnames = c("Module", "Tissue"),
                             value.name = "p_adj")
heat_df   <- merge(cor_long, pval_long, by = c("Module", "Tissue"))

# Factor levels — tissues in biological order, modules M1 → Mn top-to-bottom
heat_df$Tissue  <- factor(heat_df$Tissue,  levels = tissue_order)
heat_df$Module  <- factor(heat_df$Module,  levels = rev(module_order))  # rev → M1 at top

# Module size label (y-axis row labels)
mod_size_map <- setNames(
  sapply(module_order, function(m) {
    if (m == "Unassigned") sum(mergedColors == "grey")
    else sum(moduleLabels == m)
  }),
  module_order
)
heat_df$ModLabel <- paste0(heat_df$Module, " (n=", mod_size_map[as.character(heat_df$Module)], ")")
heat_df$ModLabel <- factor(heat_df$ModLabel,
                            levels = rev(paste0(module_order, " (n=", mod_size_map[module_order], ")")))

# ── Cell text: r value on line 1, (p-value) on line 2  — matches _06 style ──
heat_df$label <- paste0(
  signif(heat_df$r, 2), "\n(",
  as.vector(format_adjusted_p(matrix(heat_df$p_adj, ncol = 1))), ")"
)

# ── WGCNA blueWhiteRed colour palette (matches labeledHeatmap _06) ───────────
# blueWhiteRed(100) goes: deep blue → white → deep red
bwr_colors <- colorRampPalette(
  c("#00008B",   # deep blue  (blueWhiteRed low end)
    "#4444FF",
    "#8888FF",
    "#BBBBFF",
    "white",     # midpoint
    "#FFBBBB",
    "#FF8888",
    "#FF4444",
    "#8B0000")   # deep red   (blueWhiteRed high end)
)(100)

# ── Figure dimensions ─────────────────────────────────────────────────────────
gg_width  <- max(16, n_cols * 0.60 + 5)
gg_height <- max(10, n_rows * 0.50 + 4)

# ── Text colour logic: white text on deep-coloured cells, black elsewhere ─────
cell_text_color <- ifelse(abs(heat_df$r) > 0.60, "white", "black")

# ── Cell font size: scales with grid density ──────────────────────────────────
cell_fontsize <- max(2.0, min(3.2, 80 / max(n_rows, n_cols)))

gg_pub <- ggplot(heat_df, aes(x = Tissue, y = ModLabel, fill = r)) +

  # Heatmap tiles with white cell borders (clean grid)
  geom_tile(color = "white", linewidth = 0.5) +

  # r value on line 1, (p) on line 2 — same layout as _06 labeledHeatmap
  geom_text(aes(label = label),
            size       = cell_fontsize,
            lineheight = 0.82,
            fontface   = "plain",
            color      = cell_text_color) +

  # ── WGCNA blueWhiteRed colour scale ──────────────────────────────────────
  scale_fill_gradientn(
    colours = bwr_colors,
    limits  = c(-1, 1),
    name    = "Pearson\nr",
    breaks  = seq(-1, 1, 0.25),
    guide   = guide_colorbar(
      barwidth       = 1.2,
      barheight      = 12,
      title.position = "top",
      title.hjust    = 0.5,
      title.theme    = element_text(size = 11, face = "bold"),
      label.theme    = element_text(size = 10)
    )
  ) +

  # ── Publication theme: bold, large labels ─────────────────────────────────
  theme_classic(base_size = 12) +
  theme(
    # X-axis (tissue labels) — bold, angled, larger
    axis.text.x  = element_text(angle    = 45,
                                hjust    = 1,
                                vjust    = 1,
                                size     = 12,
                                face     = "bold",
                                colour   = "black"),

    # Y-axis (module labels) — bold, larger
    axis.text.y  = element_text(size   = 12,
                                face   = "bold",
                                colour = "black"),

    # Remove axis titles (tissue/module labels are self-explanatory)
    axis.title   = element_blank(),

    # Thin tick marks
    axis.ticks   = element_line(colour = "black", linewidth = 0.4),
    axis.line    = element_blank(),

    # Legend
    legend.position  = "right",
    legend.title     = element_text(size = 11, face = "bold"),
    legend.text      = element_text(size = 10),

    # Plot title and subtitle
    plot.title    = element_text(size   = 14,
                                 face   = "bold",
                                 hjust  = 0.5,
                                 margin = margin(b = 4)),
    plot.subtitle = element_text(size   = 10,
                                 colour = "grey35",
                                 hjust  = 0.5,
                                 margin = margin(b = 8)),

    # Generous margins to avoid clipping bold labels
    plot.margin  = margin(12, 16, 12, 12)
  ) +

  labs(
    title    = "Module-Tissue Correlation (Pearson r)",
    subtitle = "Upper value: Pearson r  |  Lower value: BH-adjusted p  |  Colour: r scale (-1 to +1)"
  ) +

  # Remove white space around tile grid
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0))

# ── Save: high-resolution PDF + 600 dpi PNG ───────────────────────────────────
# cairo_pdf can fail silently when the Cairo library is unavailable (for
# example, macOS without XQuartz), so the standard pdf device is the fallback.
gg_pdf_path <- out_file("WGCNA_06c_module_tissue_heatmap_ggplot_publication.pdf")
if (file.exists(gg_pdf_path)) file.remove(gg_pdf_path)
tryCatch(
  ggsave(gg_pdf_path, plot = gg_pub, width = gg_width, height = gg_height,
         units = "in", device = cairo_pdf),
  error = function(e) NULL
)
if (!file.exists(gg_pdf_path))
  ggsave(gg_pdf_path, plot = gg_pub, width = gg_width, height = gg_height,
         units = "in", device = pdf)

ggsave(out_file("WGCNA_06c_module_tissue_heatmap_ggplot_publication.png"),
       plot   = gg_pub,
       width  = gg_width,
       height = gg_height,
       units  = "in",
       dpi    = 600)

cat("Publication ggplot2 heatmap saved (PDF + 600 dpi PNG).\n")

# ── 12. Module composition bar plot (Fig 4a) ────────────────────────────────
annot <- data.frame(
  Transcript = geneNames,
  Module = moduleLabels,
  stringsAsFactors = FALSE
)

annot <- merge(annot, orig_annot, by = "Transcript", all.x = TRUE)
annot$Type <- ifelse(!is.na(annot$lncRNA) & annot$lncRNA != "",
                     "lncRNA", "PCG")

module_composition <- as.data.frame(table(annot$Module, annot$Type))
colnames(module_composition) <- c("Module", "Type", "Count")
module_composition <- module_composition[module_composition$Module != "Unassigned", ]

module_totals <- aggregate(Count ~ Module, data = module_composition, sum)
colnames(module_totals)[2] <- "Total"
module_composition <- merge(module_composition, module_totals)
module_composition$Fraction <- module_composition$Count / module_composition$Total * 100

m_levels <- paste0("M", 1:nModules)
module_composition$Module <- factor(module_composition$Module, levels = m_levels)

pdf(out_file("WGCNA_07_module_composition.pdf"), width = 14, height = 6)
p <- ggplot(module_composition, aes(x = Module, y = Fraction, fill = Type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("PCG" = "#4169E1", "lncRNA" = "#FF6347")) +
  labs(title = "Module Composition: PCGs vs lncRNAs (Fig 4a style)",
       subtitle = paste("Total modules:", nModules,
                        "| Total transcripts in modules:",
                        sum(mergedColors != "grey")),
       x = "Module", y = "Percentage of Transcripts (%)",
       fill = "Type") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data = module_totals[module_totals$Module != "Unassigned", ],
            aes(x = Module, y = -3, label = Total, fill = NULL),
            size = 2.5, inherit.aes = FALSE)
print(p)
dev.off()
cat("Module composition plot saved.\n")

# ── 13. Module eigengene bar plots ──────────────────────────────────────────
ME_df <- as.data.frame(MEs_ordered)
ME_df$Tissue <- factor(rownames(ME_df), levels = tissue_order)

ME_long <- melt(ME_df, id.vars = "Tissue",
                variable.name = "Module", value.name = "Eigengene")
ME_long <- ME_long[ME_long$Module != "Unassigned", ]

me_levels <- paste0("M", sort(as.numeric(gsub("M", "", unique(ME_long$Module)))))
ME_long$Module <- factor(ME_long$Module, levels = me_levels)

pdf(out_file("WGCNA_08_module_eigengene_profiles.pdf"), width = 20, height = 16)
eigengene_plot <- ggplot(ME_long, aes(x = Tissue, y = Eigengene)) +
  geom_bar(stat = "identity", aes(fill = Eigengene > 0), show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "firebrick3", "FALSE" = "steelblue")) +
  facet_wrap(~ Module, scales = "free_y", ncol = 4) +
  labs(title = "Module Eigengene Expression Across 32 Tissues",
       x = "", y = "Module Eigengene Value") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4),
        strip.text = element_text(face = "bold", size = 8),
        panel.grid.minor = element_blank())
print(eigengene_plot)
dev.off()
cat("Module eigengene profiles saved.\n")

# ── 14. Hub genes per module ────────────────────────────────────────────────
kME <- signedKME(datExpr, mergedMEs)

cat("\n-- Top 10 Hub Genes per Module --\n")
hub_genes_list <- list()

for (mod_color in names(colorToModule)) {
  mod_name <- colorToModule[mod_color]
  me_col   <- paste0("kME", mod_color)
  
  if (!(me_col %in% colnames(kME))) next
  
  mod_genes <- geneNames[mergedColors == mod_color]
  if (length(mod_genes) == 0) next
  mod_genes <- mod_genes[mod_genes %in% rownames(kME)]
  if (length(mod_genes) == 0) next
  
  mod_kME <- kME[mod_genes, me_col, drop = FALSE]
  mod_kME <- mod_kME[order(-mod_kME[, 1]), , drop = FALSE]
  top10   <- head(mod_kME, 10)
  
  hub_genes_list[[mod_name]] <- data.frame(
    Module = mod_name,
    Gene   = rownames(top10),
    kME    = round(top10[, 1], 4),
    row.names = NULL
  )
  
  cat(sprintf("\n%s (%s, n=%d) -- Top %d hub genes:\n",
              mod_name, mod_color, length(mod_genes), nrow(top10)))
  print(top10)
}

if (length(hub_genes_list) > 0) {
  hub_genes_df <- do.call(rbind, hub_genes_list)
  write.csv(hub_genes_df, out_file("WGCNA_hub_genes.csv"), row.names = FALSE)
  cat("\nHub genes saved to WGCNA_hub_genes.csv\n")
}

# ── 15. Save complete module assignments ────────────────────────────────────
module_assignments <- data.frame(
  Transcript   = geneNames,
  ModuleColor  = mergedColors,
  ModuleName   = moduleLabels,
  stringsAsFactors = FALSE
)
module_assignments <- merge(module_assignments, orig_annot,
                            by = "Transcript", all.x = TRUE)
module_assignments$Type <- ifelse(
  !is.na(module_assignments$lncRNA) & module_assignments$lncRNA != "",
  "lncRNA", "PCG"
)
write.csv(module_assignments, out_file("WGCNA_module_assignments.csv"),
          row.names = FALSE)
cat("Module assignments saved.\n")

# ── 16. Interactive heatmap (plotly) ────────────────────────────────────────
hover_text <- matrix(
  paste0("Module: ", rep(module_order, ncol(cor_ordered)),
         "\nTissue: ", rep(tissue_order, each = nrow(cor_ordered)),
         "\nr = ", round(cor_ordered, 3),
         "\nraw p = ", signif(pval_ordered, 3),
         "\nBH-adjusted p = ", signif(fdr_ordered, 3)),
  nrow = nrow(cor_ordered)
)

fig <- plot_ly(
  x          = tissue_order,
  y          = module_order,
  z          = cor_ordered,
  type       = "heatmap",
  colorscale = list(c(0, "#053061"), c(0.5, "white"), c(1, "#67001F")),
  zmin = -1, zmax = 1,
  text      = hover_text,
  hoverinfo = "text"
) %>%
  layout(
    title  = "Module-Tissue Correlation (hover for r, raw p and BH-adjusted p)",
    xaxis  = list(title = "", tickangle = 45, tickfont = list(size = 9)),
    yaxis  = list(title = "", tickfont = list(size = 9)),
    margin = list(b = 100, l = 100)
  )

# A self-contained HTML file requires pandoc; otherwise the supporting
# JavaScript is written to a companion folder next to the HTML file.
has_pandoc <- tryCatch(rmarkdown::pandoc_available(),
                       error = function(e) nzchar(Sys.which("pandoc")))
saveWidget(fig, out_file("WGCNA_interactive_heatmap.html"),
           selfcontained = has_pandoc)
cat("Interactive heatmap saved",
    if (has_pandoc) "(self-contained).\n" else "(with companion library folder).\n")

# ── 17. Reproducibility records ─────────────────────────────────────────────
analysis_parameters <- data.frame(
  Parameter = c("input_file", "input_md5", "tissues", "fpkm_filter",
                "variance_filter", "network_type", "soft_power",
                "observed_signed_R2", "deepSplit", "pamRespectsDendro",
                "minimum_module_size", "merge_cut_height", "p_adjustment"),
  Value = c(normalizePath(file_path), unname(tools::md5sum(file_path)),
            length(tissue_cols), "FPKM >= 0.1 in at least one tissue",
            "sample variance of log2(FPKM + 1) > 1.5", "unsigned",
            softPower, signif(observed_fit, 6), 2, FALSE, 30, 0.2,
            "Benjamini-Hochberg across assigned-module x tissue tests"),
  stringsAsFactors = FALSE
)
write.csv(analysis_parameters, out_file("WGCNA_analysis_parameters.csv"),
          row.names = FALSE)
writeLines(capture.output(sessionInfo()), out_file("WGCNA_sessionInfo.txt"))
saveRDS(list(datExpr = datExpr, geneTree = geneTree,
             dynamicColors = dynamicColors, mergedColors = mergedColors,
             moduleLabels = moduleLabels, MEs = MEs_ordered,
             moduleTraitCor = moduleTraitCor,
             moduleTraitPvalue = moduleTraitPvalue,
             moduleTraitFDR = moduleTraitFDR, kME = kME),
        out_file("WGCNA_analysis_objects.rds"))

# ── 18. Summary ─────────────────────────────────────────────────────────────
cat("\n")
cat("===============================================================\n")
cat("              WGCNA ANALYSIS SUMMARY (v3 — Publication)        \n")
cat("===============================================================\n")
cat(sprintf("  Input transcripts (after filtering):  %d\n", length(geneNames)))
cat(sprintf("  Network type:                         unsigned\n"))
cat(sprintf("  Soft-threshold power (beta):          %d\n", softPower))
cat(sprintf("  R2 at beta=%d:                         %.3f\n",
            softPower, observed_fit))
cat(sprintf("  Modules before merging:               %d\n", nPreMerge))
cat(sprintf("  Modules after merging:                %d\n", nModules))
cat(sprintf("  Assigned modules:                     %d\n", nModules))
cat(sprintf("  Unassigned genes (grey):              %d\n", sum(mergedColors == "grey")))
cat(sprintf("  Merge cut height:                     0.2\n"))
cat("===============================================================\n")
cat("\nOutput directory:", normalizePath(output_dir), "\n")
cat("  WGCNA_06_module_tissue_heatmap_publication.pdf       <- labeledHeatmap (r + p in cells)\n")
cat("  WGCNA_06b_module_tissue_heatmap_pheatmap_publication.pdf  <- pheatmap (r + sig stars)\n")
cat("  WGCNA_06c_module_tissue_heatmap_ggplot_publication.pdf    <- ggplot2 PDF (high quality)\n")
cat("  WGCNA_06c_module_tissue_heatmap_ggplot_publication.png    <- ggplot2 PNG @ 600 dpi\n")
cat("\nNOTE: All three heatmaps show tissues in their original biological order\n")
cat("  (GS → S → R → ... → Nod) with no tissue-group clustering or annotation.\n")
cat("\nDone!\n")



##############################################################################
# 19. Robustness of module assignments and hub-candidate rankings
##############################################################################
# This section extends the reference analysis above. No separate R script is
# required. All 32 profiles are assessed by leave-one-profile-out analysis.
#
# Procedures:
#   (a) Verify the reference expression matrix and reproduce its hub rankings.
#   (b) Recompute eigengenes and signed kME within fixed reference modules.
#   (c) Reconstruct all 32 omitted-profile networks and compare partitions by
#       adjusted Rand index; match modules by maximum Jaccard overlap and
#       assess hub retention within the reconstructed modules.
#   (d) Vary power, deepSplit, minimum module size and merge height individually.
#   (e) Retain the additional 20DAP exclusion with filtering reapplied as a
#       sensitivity analysis; it does not replace the 32-profile reference.
#
# The top-ten criterion uses signed kME, although adjacency/TOM are unsigned.
# Eigengene signs are aligned before reference-versus-omission hub comparisons.
# These are internal sensitivity tests, not independent biological validation.
#
# Outputs: a 'robustness' subdirectory of the reference output directory.
# Usage: Rscript "Methods S2.R" "Data S1.xlsx" "WGCNA outputs"
# To run only Sections 1-18, add --skip-robustness.
# Requires the reference WGCNA 1.74 / R 4.4.x environment for exact reproduction.
# Network reconstructions are sequential and may take considerable time.

run_coe_robustness <- function(input_file, reference_dir,
                               robustness_dir = file.path(reference_dir, "robustness"),
                               stages = c("raw_audit", "hubs", "baseline", "powers",
                                          "no20DAP_refiltered", "all_loto"),
                               summarize = TRUE) {
  reference_dir <- normalizePath(reference_dir, mustWork = TRUE)
  robustness_input_file <- normalizePath(input_file, mustWork = TRUE)
  dir.create(robustness_dir, recursive = TRUE, showWarnings = FALSE)
  robustness_dir <- normalizePath(robustness_dir, mustWork = TRUE)
  if (identical(robustness_dir, reference_dir))
    stop("Robustness outputs must use a subdirectory, not the reference directory.")
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  robustness_script_path <- if (length(file_arg))
    normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE) else
    file.path(dirname(robustness_input_file), "Methods S2.R")
  required <- c("WGCNA", "readxl", "ggplot2", "gridExtra", "jsonlite", "fastcluster")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing))
    stop("Install required robustness packages before running: ",
         paste(missing, collapse = ", "))
  reference_files <- file.path(reference_dir,
                              c("WGCNA_analysis_objects.rds", "WGCNA_hub_genes.csv"))
  if (!all(file.exists(reference_files)))
    stop("Run the reference WGCNA analysis before the robustness section.")
  # Isolated function environments keep sensitivity outputs and temporary
  # variables separate from the original reference analysis.
  run_stage <- function(stage) {
    out <- robustness_dir
    src <- dirname(robustness_input_file)
    base_dir <- reference_dir
    suppressPackageStartupMessages(library(WGCNA))
    disableWGCNAThreads()
    options(stringsAsFactors=FALSE)
    set.seed(20261019)
    writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
    stamp <- function(...) { cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n"); flush.console() }
    write_csv <- function(x,n) write.csv(x,file.path(out,n),row.names=FALSE,na="")
    b <- readRDS(file.path(base_dir,"WGCNA_analysis_objects.rds"))
    X <- as.matrix(b$datExpr)
    genes <- colnames(X)
    modules <- setNames(as.character(b$moduleLabels),genes)
    colors <- setNames(as.character(b$mergedColors),genes)
    module_names <- paste0("M",1:23)
    hubs <- read.csv(file.path(base_dir,"WGCNA_hub_genes.csv"))
    stopifnot(identical(dim(X),c(32L,13176L)), !anyNA(X), !anyDuplicated(genes))
    stopifnot(sum(modules=="Unassigned")==14L, nrow(hubs)==230L)
    md5 <- tools::md5sum(c(robustness_input_file,file.path(base_dir,"WGCNA_analysis_objects.rds"),robustness_script_path))
    write_csv(data.frame(path=names(md5),md5=unname(md5)),"source_checksums.csv")
    
    ari <- function(a,b) {
      z <- table(a,b); choose2 <- function(x) x*(x-1)/2
      n <- sum(z); s <- sum(choose2(z)); r <- sum(choose2(rowSums(z))); c <- sum(choose2(colSums(z)))
      expected <- r*c/choose2(n); denom <- (r+c)/2-expected
      if(abs(denom)<1e-12) return(if(identical(as.integer(factor(a)),as.integer(factor(b)))) 1 else NA_real_)
      (s-expected)/denom
    }
    
    # Recompute module PCs and align their arbitrary signs to the reference PCs.
    fixed_kme <- function(keep) {
      xx <- X[keep,,drop=FALSE]
      me <- moduleEigengenes(as.data.frame(xx),colors=colors,excludeGrey=TRUE)$eigengenes
      reference <- b$MEs[rownames(xx),,drop=FALSE]
      result <- vector("list",length(module_names)); names(result)<-module_names
      for(m in module_names) {
        ids <- genes[modules==m]; color <- unique(colors[ids]); eig <- me[[paste0("ME",color)]]
        alignment <- stats::cor(eig,reference[[m]])
        if(alignment<0) eig <- -eig
        k <- as.numeric(suppressWarnings(stats::cor(xx[,ids,drop=FALSE],eig))); names(k)<-ids
        result[[m]] <- list(k=k,alignment=abs(alignment))
      }
      result
    }
    
    if(stage=="hubs") {
      stamp("Reference hub reconstruction")
      ref <- fixed_kme(seq_len(nrow(X)))
      ref_all <- list(); validation <- list()
      for(m in module_names) {
        k<-ref[[m]]$k; ranked<-names(sort(k,decreasing=TRUE)); original<-hubs$Gene[hubs$Module==m]
        validation[[m]]<-data.frame(Module=m,Genes=length(k),Original_hubs_recovered=sum(head(ranked,10)%in%original),Max_abs_kME_difference=max(abs(k-as.numeric(b$kME[names(k),paste0("kME",unique(colors[names(k)]))]))))
        ref_all[[m]]<-data.frame(Module=m,Gene=names(k),Reference_kME=k,Reference_rank=rank(-k,ties.method="min"))
      }
      validation<-do.call(rbind,validation); write_csv(validation,"baseline_hub_reproduction.csv")
      stopifnot(all(validation$Original_hubs_recovered==10),max(validation$Max_abs_kME_difference)<1e-8)
      write_csv(do.call(rbind,ref_all),"reference_gene_kME.csv")
      per<-list(); each_hub<-list(); index<-0
      for(tissue in rownames(X)) {
        stamp("Leave-one-profile-out hubs:",tissue)
        kk<-fixed_kme(which(rownames(X)!=tissue))
        for(m in module_names) {
          index<-index+1; r<-ref[[m]]$k; k<-kk[[m]]$k; original<-hubs$Gene[hubs$Module==m]
          ordered<-names(sort(k,decreasing=TRUE)); top<-head(ordered,10)
          ranks<-rank(-k,ties.method="min",na.last="keep")
          ranks[!is.finite(k)]<-length(k)+1L
          per[[index]]<-data.frame(Omitted_tissue=tissue,Module=m,N_genes=length(k),Undefined_kME=sum(!is.finite(k)),Spearman_kME=stats::cor(r,k,method="spearman",use="complete.obs"),Top10_retained=sum(top%in%original),Top10_Jaccard=length(intersect(top,original))/length(union(top,original)),Eigengene_alignment=kk[[m]]$alignment,Median_original_hub_rank=median(ranks[original]),Max_original_hub_rank=max(ranks[original]))
          each_hub[[index]]<-data.frame(Omitted_tissue=tissue,Module=m,Gene=original,kME=k[original],Rank=ranks[original],In_top10=original%in%top)
        }
      }
      per<-do.call(rbind,per); each_hub<-do.call(rbind,each_hub)
      write_csv(per,"hub_leave_one_tissue_out.csv"); write_csv(each_hub,"original_hub_leave_one_tissue_out.csv")
      per_mod<-lapply(module_names,function(m) {
        z<-per[per$Module==m,]; w<-which.min(z$Spearman_kME)
        data.frame(Module=m,N_genes=unique(z$N_genes),Median_rank_correlation=median(z$Spearman_kME),Min_rank_correlation=min(z$Spearman_kME),Worst_tissue_for_ranks=z$Omitted_tissue[w],Median_top10_retained=median(z$Top10_retained),Min_top10_retained=min(z$Top10_retained),Mean_top10_retained=mean(z$Top10_retained),Top10_retained_no20DAP=z$Top10_retained[z$Omitted_tissue=="20DAP"],Rank_correlation_no20DAP=z$Spearman_kME[z$Omitted_tissue=="20DAP"])
      })
      write_csv(do.call(rbind,per_mod),"hub_stability_by_module.csv")
      byhub<-lapply(seq_len(nrow(hubs)),function(i) {
        z<-each_hub[each_hub$Module==hubs$Module[i]&each_hub$Gene==hubs$Gene[i],]
        data.frame(Module=hubs$Module[i],Gene=hubs$Gene[i],Original_kME=hubs$kME[i],Top10_retention_frequency=mean(z$In_top10),Median_rank=median(z$Rank),Worst_rank=max(z$Rank),Undefined_kME=sum(!is.finite(z$kME)),Median_kME=median(z$kME,na.rm=TRUE),Minimum_kME=min(z$kME,na.rm=TRUE))
      })
      write_csv(do.call(rbind,byhub),"hub_stability_by_gene.csv")
      saveRDS(ref,file.path(out,"reference_hub_values.rds"))
      stamp("Completed hub jackknife")
    }
    
    compare_partition <- function(newcolors,tag, elapsed, n_profiles, power, deep, minimum, merge_height) {
      scenario_genes<-names(newcolors)
      stopifnot(length(scenario_genes)==length(newcolors),!anyDuplicated(scenario_genes))
      newcolors<-setNames(as.character(newcolors),scenario_genes)
      shared<-intersect(genes,scenario_genes)
      assigned_old<-modules[shared]!="Unassigned"; assigned_new<-newcolors!="grey"
      details<-list()
      for(m in module_names) {
        oldfull<-genes[modules==m];oldset<-intersect(oldfull,shared); target_colors<-setdiff(unique(newcolors),"grey")
        jacc<-sapply(target_colors,function(cc) {newset<-shared[newcolors[shared]==cc];length(intersect(oldset,newset))/length(union(oldset,newset))})
        best<-target_colors[which.max(jacc)]; newset<-shared[newcolors[shared]==best]
        details[[m]]<-data.frame(Scenario=tag,Reference_module=m,Reference_total_size=length(oldfull),Reference_size=length(oldset),Reference_eligibility_fraction=length(oldset)/length(oldfull),Matched_color=best,Matched_size=length(newset),Intersection=length(intersect(oldset,newset)),Reference_recall=if(length(oldset))length(intersect(oldset,newset))/length(oldset)else NA_real_,Matched_precision=length(intersect(oldset,newset))/length(newset),Jaccard=max(jacc),Original_hubs_eligible=sum(hubs$Gene[hubs$Module==m]%in%scenario_genes),Original_hubs_in_matched_module=sum(hubs$Gene[hubs$Module==m]%in%newset))
      }
      details<-do.call(rbind,details)
      summary<-data.frame(Scenario=tag,N_profiles=n_profiles,N_genes=length(newcolors),Shared_reference_genes=length(shared),Reference_genes_excluded=length(setdiff(genes,scenario_genes)),New_genes=length(setdiff(scenario_genes,genes)),Power=power,DeepSplit=deep,Minimum_size=minimum,Merge_height=merge_height,Modules=length(setdiff(unique(newcolors),"grey")),Unassigned=sum(!assigned_new),ARI_all=ari(modules[shared],newcolors[shared]),ARI_reference_assigned=ari(modules[shared][assigned_old],newcolors[shared][assigned_old]),Median_module_Jaccard=median(details$Jaccard),Minimum_module_Jaccard=min(details$Jaccard),Reference_size_weighted_Jaccard=weighted.mean(details$Jaccard,details$Reference_size),Median_reference_recall=median(details$Reference_recall,na.rm=TRUE),Elapsed_seconds=elapsed)
      write_csv(summary,paste0("summary_",tag,".csv"));write_csv(details,paste0("overlap_",tag,".csv"))
      write_csv(data.frame(Gene=scenario_genes,Reference_module=modules[scenario_genes],Scenario_color=newcolors),paste0("assignments_",tag,".csv"))
      print(summary);flush.console()
      list(summary=summary,details=details,colors=newcolors)
    }
    
    network <- function(xx,power=9L,tag="network",cuts=list(c(2,30,.2)),save_diss=FALSE) {
      started<-proc.time()[3];stamp("Network start",tag,"profiles",nrow(xx),"power",power)
      good<-goodSamplesGenes(xx,verbose=0)
      if(any(!good$goodGenes)) {
        write_csv(data.frame(Gene=colnames(xx)[!good$goodGenes]),paste0("unestimable_genes_",tag,".csv"))
        xx<-xx[,good$goodGenes,drop=FALSE]
      }
      stopifnot(all(good$goodSamples))
      a<-adjacency(xx,power=power,type="unsigned")
      stamp("Adjacency complete",tag)
      tom<-TOMsimilarity(a,TOMType="unsigned",verbose=1);rm(a);gc()
      stamp("TOM complete",tag)
      diss<-1-tom;rm(tom);gc()
      tree<-fastcluster::hclust(as.dist(diss),method="average")
      if(save_diss) saveRDS(list(tree=tree,diss=diss),file.path(out,paste0("cache_",tag,".rds")),compress=FALSE)
      for(cut in cuts) {
        deep<-as.integer(cut[1]);minimum<-as.integer(cut[2]);height<-cut[3]
        dd<-cutreeDynamic(dendro=tree,distM=diss,deepSplit=deep,pamRespectsDendro=FALSE,minClusterSize=minimum,verbose=0)
        cc<-labels2colors(dd)
        merged<-mergeCloseModules(xx,cc,cutHeight=height,verbose=0)
        this_tag<-if(length(cuts)==1)tag else paste0(tag,"_d",deep,"_n",minimum,"_m",gsub("\\.","p",as.character(height)))
        if(tag=="baseline"&&deep==2L&&minimum==30L&&height==.2) {
          stopifnot(ari(modules,merged$colors)>0.999999)
          stamp("Baseline partition reproduced exactly")
        }
        compare_partition(setNames(merged$colors,colnames(xx)),this_tag,unname(proc.time()[3]-started),nrow(xx),power,deep,minimum,height)
      }
      rm(diss,tree);gc();stamp("Network finished",tag)
    }
    
    if(stage=="baseline") network(X,tag="baseline",cuts=list(c(2,30,.2),c(1,30,.2),c(3,30,.2),c(2,20,.2),c(2,50,.2),c(2,30,.15),c(2,30,.25)))
    if(stage=="no20DAP") network(X[rownames(X)!="20DAP",,drop=FALSE],tag="without20DAP_fixed_genes")
    if(stage=="power8") network(X,power=8L,tag="power8")
    if(stage=="power10") network(X,power=10L,tag="power10")
    if(stage=="no20DAP_refiltered") network(readRDS(file.path(out,"no20DAP_refiltered_expression.rds")),tag="without20DAP_refiltered")
    if(stage=="all_loto") {
      for(tissue in rownames(X)) {
        tag<-paste0("loto_",tissue)
        if(file.exists(file.path(out,paste0("summary_",tag,".csv")))) {stamp("Resume: completed",tag);next}
        network(X[rownames(X)!=tissue,,drop=FALSE],tag=tag)
      }
    }
    if(stage=="powers") {
      network(X,power=8L,tag="power8")
      network(X,power=10L,tag="power10")
    }
    
    if(stage=="raw_audit") {
      suppressPackageStartupMessages(library(readxl))
      raw_path<-robustness_input_file
      headers<-names(read_excel(raw_path,sheet="Sheet1",skip=1,n_max=0))
      select<-headers%in%c("Transcript id",rownames(X))
      types<-rep("skip",length(headers));types[select]<-"numeric";types[headers=="Transcript id"]<-"text"
      raw<-as.data.frame(read_excel(raw_path,sheet="Sheet1",skip=1,col_types=types))
      names(raw)[1]<-"Gene";raw<-raw[!duplicated(raw$Gene),]
      stopifnot(nrow(raw)==38818L)
      logx<-log2(as.matrix(raw[,rownames(X)])+1);rownames(logx)<-raw$Gene
      stopifnot(max(abs(t(logx[genes,])-X))<1e-10)
      no20<-logx[,colnames(logx)!="20DAP",drop=FALSE]
      keep32<-apply(2^logx-1,1,max)>=.1 & apply(logx,1,var)>1.5
      keep31<-apply(2^no20-1,1,max)>=.1 & apply(no20,1,var)>1.5
      stopifnot(setequal(rownames(logx)[keep32],genes))
      write_csv(data.frame(Gene=rownames(logx),Retained_32=keep32,Retained_31=keep31,Variance_32=apply(logx,1,var),Variance_31=apply(no20,1,var)),"20DAP_filter_sensitivity.csv")
      write_csv(data.frame(Tissue=colnames(logx),Mean_FPKM=colMeans(2^logx-1),Mean_log2_FPKM1=colMeans(logx),Median_FPKM=apply(2^logx-1,2,median)),"tissue_expression_scale.csv")
      write_csv(data.frame(Original_32_retained=sum(keep32),Without20DAP_retained=sum(keep31),Shared=sum(keep32&keep31),Lost=sum(keep32&!keep31),Gained=sum(!keep32&keep31)),"20DAP_filter_summary.csv")
      saveRDS(t(no20[keep31,,drop=FALSE]),file.path(out,"no20DAP_refiltered_expression.rds"))
      stamp("Raw expression verified; retained without20DAP",sum(keep31),"lost",sum(keep32&!keep31),"gained",sum(!keep32&keep31))
    }
    
    if(stage=="collect") {
      ff<-list.files(out,pattern="^summary_.*\\.csv$",full.names=TRUE)
      write_csv(do.call(rbind,lapply(ff,read.csv)),"network_sensitivity_summary.csv")
      ff<-list.files(out,pattern="^overlap_.*\\.csv$",full.names=TRUE)
      write_csv(do.call(rbind,lapply(ff,read.csv)),"network_sensitivity_by_module.csv")
    }
  }
  summarize_outputs <- function() {
    out <- robustness_dir
    suppressPackageStartupMessages(library(WGCNA))
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(gridExtra))
    options(stringsAsFactors=FALSE)
    read_csv<-function(n)read.csv(file.path(out,n),check.names=FALSE)
    write_csv<-function(x,n)write.csv(x,file.path(out,n),row.names=FALSE,na="")
    ff<-list.files(out,pattern="^summary_.*\\.csv$",full.names=TRUE)
    summary<-do.call(rbind,lapply(ff,read.csv));stopifnot(nrow(summary)==42L)
    write_csv(summary,"network_sensitivity_summary.csv")
    ff<-list.files(out,pattern="^overlap_.*\\.csv$",full.names=TRUE)
    overlap<-do.call(rbind,lapply(ff,read.csv));stopifnot(nrow(overlap)==42L*23L)
    write_csv(overlap,"network_sensitivity_by_module.csv")
    b<-readRDS(file.path(reference_dir,"WGCNA_analysis_objects.rds"))
    X<-as.matrix(b$datExpr);genes<-colnames(X);tissues<-rownames(X);modules<-paste0("M",1:23)
    hubs<-read.csv(file.path(reference_dir,"WGCNA_hub_genes.csv"))
    jack<-read_csv("hub_leave_one_tissue_out.csv")
    hubgene<-read_csv("hub_stability_by_gene.csv")
    hubmod<-read_csv("hub_stability_by_module.csv")
    stopifnot(nrow(jack)==736L,nrow(hubgene)==230L,nrow(read_csv("original_hub_leave_one_tissue_out.csv"))==7360L)
    stopifnot(all(read_csv("baseline_hub_reproduction.csv")$Original_hubs_recovered==10))
    baseline<-summary[summary$Scenario=="baseline_d2_n30_m0p2",];stopifnot(baseline$ARI_all==1)
    loto<-summary[grepl("^loto_",summary$Scenario),];stopifnot(nrow(loto)==32L)
    loto$Omitted_tissue<-sub("^loto_","",loto$Scenario)
    loto<-loto[match(tissues,loto$Omitted_tissue),]
    write_csv(loto,"module_leave_one_tissue_out_summary.csv")
    
    # Re-rank hubs within each independently reconstructed, overlap-matched module.
    # PC signs are aligned to the original module PC on the retained profiles.
    rebuilt<-list();index<-0
    assignment_retention<-matrix(0L,nrow=length(genes),ncol=length(tissues),dimnames=list(genes,tissues))
    assignment_retention[as.character(b$moduleLabels)=="Unassigned",]<-NA_integer_
    for(tissue in tissues) {
      tag<-paste0("loto_",tissue)
      assignment<-read_csv(paste0("assignments_",tag,".csv"))
      xx<-X[rownames(X)!=tissue,assignment$Gene,drop=FALSE]
      me<-moduleEigengenes(as.data.frame(xx),colors=assignment$Scenario_color,excludeGrey=TRUE)$eigengenes
      details<-overlap[overlap$Scenario==tag,]
      for(m in modules) {
        index<-index+1;row<-details[details$Reference_module==m,]
        ids<-assignment$Gene[assignment$Scenario_color==row$Matched_color]
        original_members<-genes[as.character(b$moduleLabels)==m]
        assignment_retention[intersect(original_members,ids),tissue]<-1L
        pc<-me[[paste0("ME",row$Matched_color)]]
        alignment<-stats::cor(pc,b$MEs[rownames(xx),m])
        if(alignment<0)pc<- -pc
        k<-as.numeric(stats::cor(xx[,ids,drop=FALSE],pc));names(k)<-ids
        top<-names(head(sort(k,decreasing=TRUE),10))
        original<-hubs$Gene[hubs$Module==m]
        rebuilt[[index]]<-data.frame(Omitted_tissue=tissue,Module=m,Matched_color=row$Matched_color,Matched_module_size=length(ids),Eigengene_alignment=abs(alignment),Original_hubs_in_matched_module=sum(original%in%ids),Original_hubs_in_rebuilt_top10=sum(original%in%top),Top10_Jaccard=length(intersect(original,top))/length(union(original,top)))
      }
    }
    rebuilt<-do.call(rbind,rebuilt)
    write_csv(rebuilt,"hub_stability_after_network_reconstruction.csv")
    write_csv(data.frame(Gene=genes,Reference_module=as.character(b$moduleLabels),Matched_module_retention_frequency=rowMeans(assignment_retention),Retained_in_matched_module_without20DAP=assignment_retention[,"20DAP"],Original_hub=genes%in%hubs$Gene),"module_assignment_stability_by_gene.csv")
    mod_summary<-do.call(rbind,lapply(modules,function(m) {
      z<-overlap[overlap$Reference_module==m&grepl("^loto_",overlap$Scenario),]
      q<-rebuilt[rebuilt$Module==m,]
      data.frame(Module=m,Original_size=unique(z$Reference_total_size),Median_module_Jaccard=median(z$Jaccard),Minimum_module_Jaccard=min(z$Jaccard),Worst_omission=sub("^loto_","",z$Scenario[which.min(z$Jaccard)]),Median_module_recall=median(z$Reference_recall),Minimum_module_recall=min(z$Reference_recall),Median_rebuilt_top10_retained=median(q$Original_hubs_in_rebuilt_top10),Minimum_rebuilt_top10_retained=min(q$Original_hubs_in_rebuilt_top10),Rebuilt_top10_retained_without20DAP=q$Original_hubs_in_rebuilt_top10[q$Omitted_tissue=="20DAP"])
    }))
    write_csv(mod_summary,"module_stability_by_module.csv")
    combined<-merge(hubmod,mod_summary,by="Module",sort=FALSE)
    combined<-combined[match(modules,combined$Module),]
    write_csv(combined,"Supplementary_WGCNA_module_robustness.csv")
    
    # Preserve the complete refiltered sensitivity result separately from Data S11.
    # S labels are deliberately distinct from the original M1-M23 reference labels.
    alternative<-read_csv("assignments_without20DAP_refiltered.csv")
    xx_alt<-readRDS(file.path(out,"no20DAP_refiltered_expression.rds"))
    xx_alt<-xx_alt[,alternative$Gene,drop=FALSE]
    alt_sizes<-sort(table(alternative$Scenario_color[alternative$Scenario_color!="grey"]),decreasing=TRUE)
    alt_key<-data.frame(Sensitivity_module=paste0("S",seq_along(alt_sizes)),Color=names(alt_sizes),Transcripts=as.integer(alt_sizes))
    alt_map<-setNames(alt_key$Sensitivity_module,alt_key$Color)
    alternative$Sensitivity_module<-ifelse(alternative$Scenario_color=="grey","Unassigned",alt_map[alternative$Scenario_color])
    alt_me<-moduleEigengenes(as.data.frame(xx_alt),colors=alternative$Scenario_color,excludeGrey=TRUE)$eigengenes
    alt_hubs<-do.call(rbind,lapply(seq_len(nrow(alt_key)),function(i) {
      ids<-alternative$Gene[alternative$Scenario_color==alt_key$Color[i]]
      k<-as.numeric(stats::cor(xx_alt[,ids,drop=FALSE],alt_me[[paste0("ME",alt_key$Color[i])]]));names(k)<-ids
      top<-head(sort(k,decreasing=TRUE),10)
      data.frame(Sensitivity_module=alt_key$Sensitivity_module[i],Color=alt_key$Color[i],Gene=names(top),Signed_kME=as.numeric(top))
    }))
    write_csv(alternative,"without20DAP_refiltered_module_assignments.csv")
    write_csv(alt_key,"without20DAP_refiltered_module_key.csv")
    write_csv(alt_hubs,"without20DAP_refiltered_hub_candidates.csv")
    
    filter<-read_csv("20DAP_filter_summary.csv")
    scale<-read_csv("tissue_expression_scale.csv")
    without20<-loto[loto$Omitted_tissue=="20DAP",]
    without20refiltered<-summary[summary$Scenario=="without20DAP_refiltered",]
    parameters<-summary[!grepl("^loto_",summary$Scenario)&summary$Scenario!="without20DAP_refiltered"&summary$Scenario!="baseline_d2_n30_m0p2",]
    global<-list(
      baseline=list(profiles=nrow(X),transcripts=ncol(X),assigned=sum(b$moduleLabels!="Unassigned"),unassigned=sum(b$moduleLabels=="Unassigned"),modules=23,hubs=nrow(hubs),reproduced_ARI=baseline$ARI_all),
      fixed_module_hubs=list(comparisons=nrow(jack),median_spearman=median(jack$Spearman_kME),minimum_spearman=min(jack$Spearman_kME),median_top10_retained=median(jack$Top10_retained),original_hubs_retained_in_all_32=sum(hubgene$Top10_retention_frequency==1),original_hubs_retained_in_at_least_90_percent=sum(hubgene$Top10_retention_frequency>=.9),without20DAP_retained=sum(jack$Top10_retained[jack$Omitted_tissue=="20DAP"]),undefined_kME_comparisons=sum(jack$Undefined_kME)),
      reconstructed_networks=list(jackknife_runs=nrow(loto),median_ARI=median(loto$ARI_all),min_ARI=min(loto$ARI_all),max_ARI=max(loto$ARI_all),worst_omission=loto$Omitted_tissue[which.min(loto$ARI_all)],min_modules=min(loto$Modules),max_modules=max(loto$Modules),median_module_Jaccard=median(overlap$Jaccard[grepl("^loto_",overlap$Scenario)]),median_top10_retained=median(rebuilt$Original_hubs_in_rebuilt_top10),median_total_hubs_retained=median(tapply(rebuilt$Original_hubs_in_rebuilt_top10,rebuilt$Omitted_tissue,sum)),without20DAP_ARI=without20$ARI_all,without20DAP_modules=without20$Modules,without20DAP_hubs_retained=sum(rebuilt$Original_hubs_in_rebuilt_top10[rebuilt$Omitted_tissue=="20DAP"])),
      no20DAP_refiltered=list(transcripts=filter$Without20DAP_retained,shared=filter$Shared,lost=filter$Lost,gained=filter$Gained,modules=without20refiltered$Modules,unassigned=without20refiltered$Unassigned,ARI_common_genes=without20refiltered$ARI_all,median_module_Jaccard_common=without20refiltered$Median_module_Jaccard),
      parameter_sensitivity=list(runs=nrow(parameters),min_ARI=min(parameters$ARI_all),max_ARI=max(parameters$ARI_all),min_modules=min(parameters$Modules),max_modules=max(parameters$Modules)),
      raw20DAP_mean_ratio_to_next_highest=scale$Mean_FPKM[scale$Tissue=="20DAP"]/max(scale$Mean_FPKM[scale$Tissue!="20DAP"])
    )
    jsonlite::write_json(global,file.path(out,"results_summary.json"),pretty=TRUE,auto_unbox=TRUE,digits=12)
    print(global)
    
    # Publication figures. All axes retain their full natural range.
    theme_pub<-theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face="bold"),axis.text=element_text(colour="black"))
    loto$Tissue<-factor(loto$Omitted_tissue,levels=tissues)
    loto$Highlight<-ifelse(loto$Omitted_tissue=="20DAP","20DAP omitted","Other profile omitted")
    p1<-ggplot(loto,aes(Tissue,ARI_all,fill=Highlight))+geom_col(width=.75)+scale_fill_manual(values=c("20DAP omitted"="#d95f02","Other profile omitted"="#2166ac"))+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.2))+labs(title="A  Module assignments after network reconstruction",x=NULL,y="Adjusted Rand index",fill=NULL)+theme_pub+theme(axis.text.x=element_text(angle=60,hjust=1),legend.position="top")
    total_hubs<-aggregate(Top10_retained~Omitted_tissue,jack,sum)
    total_hubs$Tissue<-factor(total_hubs$Omitted_tissue,levels=tissues)
    total_hubs$Fraction<-total_hubs$Top10_retained/230
    total_hubs$Highlight<-ifelse(total_hubs$Omitted_tissue=="20DAP","20DAP omitted","Other profile omitted")
    p2<-ggplot(total_hubs,aes(Tissue,Fraction,fill=Highlight))+geom_col(width=.75)+scale_fill_manual(values=c("20DAP omitted"="#d95f02","Other profile omitted"="#2166ac"))+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.2))+labs(title="B  Original top-ten hubs retained within fixed modules",x="Omitted tissue or developmental-stage profile",y="Fraction of 230 original hubs",fill=NULL)+theme_pub+theme(axis.text.x=element_text(angle=60,hjust=1),legend.position="none")
    g1<-arrangeGrob(p1,p2,ncol=1)
    ggsave(file.path(out,"Figure_WGCNA_profile_sensitivity.pdf"),g1,width=13,height=8.5)
    ggsave(file.path(out,"Figure_WGCNA_profile_sensitivity.png"),g1,width=13,height=8.5,dpi=300)
    
    heat_module<-overlap[grepl("^loto_",overlap$Scenario),]
    heat_module$Tissue<-factor(sub("^loto_","",heat_module$Scenario),levels=tissues)
    heat_module$Module<-factor(heat_module$Reference_module,levels=rev(modules))
    heat_hubs<-jack;heat_hubs$Tissue<-factor(heat_hubs$Omitted_tissue,levels=tissues);heat_hubs$Module<-factor(heat_hubs$Module,levels=rev(modules))
    p3<-ggplot(heat_module,aes(Tissue,Module,fill=Jaccard))+geom_tile(colour="white",linewidth=.15)+scale_fill_gradient(limits=c(0,1),low="white",high="#2166ac",name="Module\nJaccard")+labs(title="A  Module overlap after independent reconstruction",x=NULL,y=NULL)+theme_pub+theme(axis.text.x=element_text(angle=60,hjust=1),panel.grid=element_blank())
    p4<-ggplot(heat_hubs,aes(Tissue,Module,fill=Top10_retained))+geom_tile(colour="white",linewidth=.15)+scale_fill_gradient(limits=c(0,10),low="white",high="#2166ac",name="Original hubs\nretained / 10")+labs(title="B  Hub retention when module membership is held fixed",x="Omitted tissue or developmental-stage profile",y=NULL)+theme_pub+theme(axis.text.x=element_text(angle=60,hjust=1),panel.grid=element_blank())
    g2<-arrangeGrob(p3,p4,ncol=1)
    ggsave(file.path(out,"Figure_WGCNA_module_hub_stability.pdf"),g2,width=14,height=13)
    ggsave(file.path(out,"Figure_WGCNA_module_hub_stability.png"),g2,width=14,height=13,dpi=300)
    
    # Cross-check ARI with an independent contingency-table implementation.
    py_checks<-do.call(rbind,lapply(c("baseline_d2_n30_m0p2","loto_20DAP","without20DAP_refiltered","power8"),function(tag) {
      aa<-read_csv(paste0("assignments_",tag,".csv"));aa<-aa[!is.na(aa$Reference_module)&nzchar(aa$Reference_module),]
      z<-table(aa$Reference_module,aa$Scenario_color);C2<-function(x)x*(x-1)/2
      expected<-sum(C2(rowSums(z)))*sum(C2(colSums(z)))/C2(sum(z))
      independent<-(sum(C2(z))-expected)/(.5*(sum(C2(rowSums(z)))+sum(C2(colSums(z))))-expected)
      reported<-summary$ARI_all[summary$Scenario==tag]
      data.frame(Scenario=tag,Reported_ARI=reported,Recomputed_ARI=independent,Absolute_difference=abs(reported-independent))
    }))
    stopifnot(all(py_checks$Absolute_difference<1e-12))
    write_csv(py_checks,"verification_ARI_recalculation.csv")
    write_csv(data.frame(Check=c("Baseline partition reproduced","Baseline 230 hubs reproduced","32 independent tissue omissions completed","736 fixed-module comparisons completed","736 reconstructed-module hub comparisons completed","42 scenario summaries complete","Independent ARI arithmetic matches"),Passed=TRUE),"verification_checks.csv")
  }
  for (stage_name in stages) {
    cat("\nRobustness stage:", stage_name, "\n")
    run_stage(stage_name)
    gc()
  }
  if (isTRUE(summarize)) summarize_outputs()
  cat("\nRobustness analysis completed. Outputs:", robustness_dir, "\n")
  invisible(robustness_dir)
}

# Release the reference network's dense matrices before sequential rebuilds.
# The saved reference objects and other reference outputs remain available.
if (!("--skip-robustness" %in% args)) {
  dense_objects <- intersect(c("TOM", "dissTOM"), ls(envir = .GlobalEnv))
  if (length(dense_objects)) rm(list = dense_objects, envir = .GlobalEnv)
  gc()
  run_coe_robustness(file_path, output_dir)
}
