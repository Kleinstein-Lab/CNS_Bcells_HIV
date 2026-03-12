library(Signac)
library(Seurat)
library(ArchR)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(IRanges)
library(S4Vectors)
library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(igraph)
library(ggraph)
library(tidygraph)
library(ComplexHeatmap)
library(circlize)

addArchRThreads(threads = 128) 
addArchRGenome("hg38")


set.seed(42)
#-----------------------------------------------------------------------------
# Differential genes from RNA
#-----------------------------------------------------------------------------
#diff_genes <- ATAC_candidate_genes_Jan29$gene
diff_genes <- c(ATAC_candidate_genes_Jan29$gene[which(ATAC_candidate_genes_Jan29$p_val_adj_multiome < 0.05)],
                ATAC_candidate_genes_Jan29$gene[which(ATAC_candidate_genes_Jan29$p_val_adj_5prime < 0.05)])
diff_genes <- unique(diff_genes)
length(diff_genes)

#-----------------------------------------------------------------------------
# Integrate RNA data to ATAC ArchR Project
#-----------------------------------------------------------------------------
rna_seurat <- readRDS("./refined_rna_seurat.RDS")
rna_seurat_Bcell <- rna_seurat[,which(rna_seurat$tissue == "CSF")]
rna_seurat_Bcell <- rna_seurat_Bcell[,which(rna_seurat_Bcell$annotated_clusters_all == "B cells")]

proj_atac_csf<- loadArchRProject("./", force = TRUE)
idxPass <- which(proj_atac_csf$rna_celltype == "B cells")
atac_pass_filtered <- proj_atac_csf$cellNames[idxPass]
proj_atac_celltype <- proj_atac_csf[atac_pass_filtered, ]

RM_Bcell <- rna_seurat$new_cell_name[which(rna_seurat$Bcell_subtype %in% c("Resting Memory"))]
proj_atac_celltype <- proj_atac_celltype[which(proj_atac_celltype$cellNames %in% RM_Bcell)]

seRNA <- as.SingleCellExperiment(rna_seurat_Bcell)
seRNA <- seRNA[,which(seRNA$new_cell_name %in% proj_atac_celltype$cellNames)]
genomeAnnotation <- getGenomeAnnotation(proj_atac_celltype)
geneAnnotation <- getGeneAnnotation(proj_atac_celltype)
genes <- geneAnnotation$genes
sce_genes <- rownames(seRNA)
match_idx <- match(sce_genes, genes$symbol)
n_matched <- sum(!is.na(match_idx))
message("Matched ", n_matched, " out of ", length(sce_genes), " genes (", 
        round(n_matched/length(sce_genes)*100, 2), "%)")
matched_genes <- sce_genes[!is.na(match_idx)]
matched_idx <- match_idx[!is.na(match_idx)]
matched_gene_info <- genes[matched_idx, ]
gr <- GRanges(
  seqnames = seqnames(matched_gene_info),
  ranges = IRanges(start = start(matched_gene_info), end = end(matched_gene_info)),
  strand = strand(matched_gene_info)
)
names(gr) <- matched_genes

unmatched_genes <- sce_genes[is.na(match_idx)]
if (length(unmatched_genes) > 0) {
  message("Creating placeholder entries for ", length(unmatched_genes), " unmatched genes")
  
  placeholder_gr <- GRanges(
    seqnames = rep("chr1", length(unmatched_genes)),
    ranges = IRanges(start = rep(1, length(unmatched_genes)), 
                     end = rep(1, length(unmatched_genes))),
    strand = rep("*", length(unmatched_genes))
  )
  names(placeholder_gr) <- unmatched_genes
  
  # Combine matched and placeholder GRanges
  all_gr <- c(gr, placeholder_gr)
  
  # Sort to match original order
  all_gr <- all_gr[sce_genes]
} else {
  all_gr <- gr
}

# Add the GRanges to the SCE object
rowRanges(seRNA) <- all_gr
seRNA <- seRNA[matched_genes,]
colnames(seRNA) <- seRNA$new_cell_name
seRNA <- seRNA[, proj_atac_celltype$cellNames]
table(seRNA$new_cell_name == proj_atac_celltype$cellNames)
proj_atac_celltype <- addGeneExpressionMatrix(input = proj_atac_celltype, seRNA = seRNA, strictMatch = TRUE, force = TRUE)


proj_atac_celltype <- addIterativeLSI(
  ArchRProj = proj_atac_celltype,
  useMatrix = "PeakMatrix", 
  name = "IterativeLSI_peak", 
  iterations = 2, 
  clusterParams = list( 
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  force = TRUE
)


#-----------------------------------------------------------------------------
# Add peak to gene links
#-----------------------------------------------------------------------------
proj_atac_celltype <- addPeak2GeneLinks(
  ArchRProj = proj_atac_celltype,
  reducedDims = "IterativeLSI_peak",
  useMatrix = "GeneExpressionMatrix",
  k = 20,
  maxDist = 50000
)


#=============================================================================
# PART 1: chromVAR - Differential TF Motif Analysis
#=============================================================================

#-----------------------------------------------------------------------------
# 1.1 Add motif annotations
#-----------------------------------------------------------------------------

proj_atac_celltype <- addMotifAnnotations(
  ArchRProj = proj_atac_celltype,
  motifSet = "cisbp",
  name = "Motif",
  force = FALSE
)

#-----------------------------------------------------------------------------
# 1.2 Add background peaks for chromVAR
#-----------------------------------------------------------------------------

proj_atac_celltype <- addBgdPeaks(proj_atac_celltype, force = TRUE)

#-----------------------------------------------------------------------------
# 1.3 Add chromVAR deviations matrix
#-----------------------------------------------------------------------------

proj_atac_celltype <- addDeviationsMatrix(
  ArchRProj = proj_atac_celltype,
  peakAnnotation = "Motif",
  force = TRUE
)
#-----------------------------------------------------------------------------
# 1.4 Differential TF motif accessibility: PWH vs Control
#-----------------------------------------------------------------------------

markerMotifs <- getMarkerFeatures(
  ArchRProj = proj_atac_celltype,
  useMatrix = "MotifMatrix",
  groupBy = "HIV_status",
  useGroups = "PWH",
  bgdGroups = "control",
  bias = c("TSSEnrichment", "log10(nFrags)"),  # Adjust for technical biases
  testMethod = "wilcoxon"
)

motif_results <- getMarkers(markerMotifs, cutOff = "FDR <= 1")  # Get all results
motif_df <- as.data.frame(motif_results$PWH)

# View top differential motifs
head(motif_df[order(motif_df$FDR), ], 20)

motif_df$TF <- gsub("_.*", "", motif_df$name)

p <- ggplot(motif_df, aes(x = MeanDiff, y = -log10(FDR))) +
  geom_point(aes(color = FDR < 0.1), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("grey50", "red"), name = "FDR < 0.1") +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "blue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_text_repel(
    data = motif_df %>% filter(FDR < 0.1),
    aes(label = TF),
    size = 3, max.overlaps = 30
  ) +
  labs(
    x = "MeanDiff (PWH - Control)",
    y = "-log10(FDR)",
    title = "Differential TF Motif Accessibility in CSF Resting Memory B cells",
    subtitle = ""
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom"
  )

write_csv(motif_df, "./CSF_B_Motif_disease.csv")

png(filename = "./Bcell-differential_motif_CSF_disease.png",  bg = "transparent",width = 7, height = 6, units = "in", res = 600)
p
dev.off()
#-----------------------------------------------------------------------------
# 2.1 Extract unique TFs from  marker motifs
#-----------------------------------------------------------------------------
diff_motif_df <- motif_df %>% filter(FDR < 0.1)

rna_seurat_Bcell_motif <- rna_seurat_Bcell[diff_motif_df$TF,]
rowSums(rna_seurat_Bcell_motif@assays$RNA$counts)


#=============================================================================
# STEP 2: PEAK-TO-GENE LINKS
# Find which peaks are correlated with which genes
#=============================================================================
# Extract peak-to-gene links

p2g <- getPeak2GeneLinks(proj_atac_celltype,corCutOff = 0,FDRCutOff = 0.05,returnLoops = TRUE)
p2g$Peak2GeneLinks <- p2g$Peak2GeneLinks[which(abs(p2g$Peak2GeneLinks$value) > 0.25),]

p2g <- getPeak2GeneLinks(proj_atac_celltype, returnLoops = FALSE)

p2geneDF <- metadata(proj_atac_celltype@peakSet)$Peak2GeneLinks
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]


p2geneDF <-  as.data.frame(p2geneDF)
p2geneDF <- p2geneDF %>% filter(FDR < 0.1)
p2geneDF <- p2geneDF %>% filter(abs(Correlation) > 0.25)
p2geneDF <- p2geneDF %>% filter(VarQATAC > 0.25)
p2geneDF <- p2geneDF %>% filter(VarQRNA > 0.25)


# Summary statistics
cat("\n=== Peak-to-Gene Link Summary ===\n")
cat("Total links:", nrow(p2geneDF), "\n")
cat("Unique peaks with links:", length(unique(p2geneDF$peakName)), "\n")
cat("Unique genes with links:", length(unique(p2geneDF$gene)), "\n")

# Check correlation distribution
cat("\nCorrelation distribution:\n")
print(summary(p2geneDF$Correlation))

# View some examples
cat("\nExample peak-to-gene links:\n")
print(head(p2geneDF %>% dplyr::select(peakName, geneName, Correlation), 10))
diff_genes_with_links <- intersect(diff_genes, unique(p2geneDF$geneName))
cat("\n===  Differential genes with peak links ===\n")
cat(length(diff_genes_with_links), "out of", length(diff_genes), "have peak-to-gene links\n")
cat("Genes with links:", paste(diff_genes_with_links, collapse = ", "), "\n")


#=============================================================================
# STEP 1: MOTIF MATCHING
# Find peaks that contain each TF's motif
#=============================================================================

#=============================================================================
# STEP 1: MOTIF MATCHING
# Use  differential motifs 
#=============================================================================

#  differential motifs from the analysis
diff_motif_names <- diff_motif_df$name

# Get the motif-peak match matrix
motif_positions <- getMatches(proj_atac_celltype, name = "Motif")
motif_match_mat <- assay(motif_positions)

rownames(motif_match_mat) <- paste0(
  seqnames(rowRanges(motif_positions)), ":",
  start(rowRanges(motif_positions)), "-",
  end(rowRanges(motif_positions))
)

cat("Motif match matrix:", nrow(motif_match_mat), "peaks x", ncol(motif_match_mat), "motifs\n")

# Check which of differential motifs are in the matrix
motifs_found <- intersect(diff_motif_names, colnames(motif_match_mat))
motifs_missing <- setdiff(diff_motif_names, colnames(motif_match_mat))

cat("\nDifferential motifs found in matrix:", length(motifs_found), "/", length(diff_motif_names), "\n")

if (length(motifs_missing) > 0) {
  cat("Missing motifs:", paste(motifs_missing, collapse = ", "), "\n")
}

# Count peaks per differential motif
cat("\n=== Peaks containing each differential motif ===\n")

motif_peak_counts <- data.frame(
  motif = motifs_found,
  TF = gsub("_.*", "", motifs_found),
  n_peaks = sapply(motifs_found, function(m) sum(motif_match_mat[, m]))
)

motif_peak_counts <- motif_peak_counts %>% arrange(desc(n_peaks))
print(motif_peak_counts)

# Summary
cat("\nTotal peaks across all differential motifs:\n")
cat("  Range:", min(motif_peak_counts$n_peaks), "-", max(motif_peak_counts$n_peaks), "\n")
cat("  Median:", median(motif_peak_counts$n_peaks), "\n")


#=============================================================================
# STEP 3: TF TARGET INFERENCE
# Combine: TF motif in peak + peak linked to gene = TF target
#=============================================================================

cat("Columns in p2geneDF:\n")
print(colnames(p2geneDF))

# Add peak_name column if not present (adjust column names as needed)
if (!"peak_name" %in% colnames(p2geneDF)) {
  if ("idxATAC" %in% colnames(p2geneDF)) {
    peak_set <- getPeakSet(proj_atac_celltype)
    p2geneDF$peak_name <- paste0(
      seqnames(peak_set)[p2geneDF$idxATAC], ":",
      start(peak_set)[p2geneDF$idxATAC], "-",
      end(peak_set)[p2geneDF$idxATAC]
    )
  }
}

cat("\nFirst few rows of p2geneDF:\n")
print(head(p2geneDF))

#-----------------------------------------------------------------------------
# Function to get target genes for each differential motif
#-----------------------------------------------------------------------------

get_motif_targets <- function(motif_name, motif_match_mat, p2geneDF) {
  
  # Find peaks containing this motif
  peaks_with_motif <- rownames(motif_match_mat)[motif_match_mat[, motif_name]]
  
  # Find genes linked to those peaks
  gene_col <- ifelse("geneName" %in% colnames(p2geneDF), "geneName", "gene")
  peak_col <- ifelse("peak_name" %in% colnames(p2geneDF), "peak_name", "peakName")
  
  target_genes <- p2geneDF %>%
    filter(.data[[peak_col]] %in% peaks_with_motif) %>%
    pull(.data[[gene_col]]) %>%
    unique()
  
  list(
    motif = motif_name,
    TF = gsub("_.*", "", motif_name),
    n_peaks_with_motif = length(peaks_with_motif),
    n_target_genes = length(target_genes),
    targets = target_genes
  )
}

#-----------------------------------------------------------------------------
# Build target lists for all differential motifs
#-----------------------------------------------------------------------------

cat("=== Building motif → target gene lists ===\n\n")

motif_targets_list <- list()

for (motif in diff_motif_names) {
  result <- get_motif_targets(motif, motif_match_mat, p2geneDF)
  motif_targets_list[[motif]] <- result
  
  cat(result$TF, "(", motif, "):", 
      result$n_peaks_with_motif, "peaks →", 
      result$n_target_genes, "target genes\n")
}

#-----------------------------------------------------------------------------
# Summary table
#-----------------------------------------------------------------------------

motif_target_summary <- data.frame(
  motif = names(motif_targets_list),
  TF = sapply(motif_targets_list, function(x) x$TF),
  n_peaks = sapply(motif_targets_list, function(x) x$n_peaks_with_motif),
  n_targets = sapply(motif_targets_list, function(x) x$n_target_genes)
) %>% arrange(desc(n_targets))

cat("\n=== Motif Target Summary ===\n")
print(motif_target_summary)

#-----------------------------------------------------------------------------
# Define background
#-----------------------------------------------------------------------------

all_linked_genes <- unique(p2geneDF$geneName)
n_background <- length(all_linked_genes)

diff_genes_in_background <- intersect(diff_genes, all_linked_genes)
n_diff_in_bg <- length(diff_genes_in_background)

cat("=== Enrichment Test Setup ===\n")
cat("Background (all genes with peak links):", n_background, "\n")
cat("Diff genes in background:", n_diff_in_bg, "out of", length(diff_genes), "\n")
cat("Diff genes with links:", paste(diff_genes_in_background, collapse = ", "), "\n\n")

#-----------------------------------------------------------------------------
# Run Fisher's Exact Test for each motif
#-----------------------------------------------------------------------------

cat("=== Running Fisher's Exact Test ===\n\n")

enrichment_results <- lapply(names(motif_targets_list), function(motif) {
  
  targets <- motif_targets_list[[motif]]$targets
  tf <- motif_targets_list[[motif]]$TF
  
  # Overlap with differential genes
  overlap_genes <- intersect(targets, diff_genes)
  
  a <- length(overlap_genes)                                    # TF targets AND diff genes
  b <- length(setdiff(targets, diff_genes))                     # TF targets but NOT diff genes
  c <- length(setdiff(diff_genes_in_background, targets))       # Diff genes but NOT TF targets
  d <- n_background - a - b - c                                 # Neither
  
  # Expected overlap by chance
  expected <- (a + b) * (a + c) / n_background
  
  # Fisher's exact test
  if (a + b > 0 & a + c > 0 & d > 0) {
    fisher_res <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
    
    data.frame(
      motif = motif,
      TF = tf,
      n_targets = length(targets),
      n_diff_targeted = a,
      expected = round(expected, 2),
      fold_enrichment = round(a / max(expected, 0.01), 2),
      pval = fisher_res$p.value,
      odds_ratio = round(fisher_res$estimate, 2),
      overlap_genes = paste(overlap_genes, collapse = "; ")
    )
  } else {
    NULL
  }
}) %>% bind_rows()

# Adjust p-values for multiple testing
enrichment_results$padj <- p.adjust(enrichment_results$pval, method = "BH")

# Sort by p-value
enrichment_results <- enrichment_results %>% arrange(padj)

#-----------------------------------------------------------------------------
# Display results
#-----------------------------------------------------------------------------

cat("=== ENRICHMENT RESULTS ===\n\n")
print(enrichment_results %>% 
        dplyr::select(motif, TF, n_targets, n_diff_targeted, expected, 
                      fold_enrichment, pval, padj) %>%
        head(25))

#-----------------------------------------------------------------------------
# Highlight significant results
#-----------------------------------------------------------------------------

sig_enriched <- enrichment_results %>% filter(padj < 0.1)

cat("\n=== SIGNIFICANT ENRICHMENTS (p < 0.05) ===\n")

if (nrow(sig_enriched) > 0) {
  for (i in 1:nrow(sig_enriched)) {
    cat("\n----------------------------------------\n")
    cat(sig_enriched$TF[i], "(", sig_enriched$motif[i], ")\n")
    cat("  Total targets:", sig_enriched$n_targets[i], "\n")
    cat("  Diff genes targeted:", sig_enriched$n_diff_targeted[i], 
        "(expected:", sig_enriched$expected[i], ")\n")
    cat("  Fold enrichment:", sig_enriched$fold_enrichment[i], "x\n")
    cat("  P-value:", format(sig_enriched$pval[i], digits = 3), "\n")
    cat("  Adjusted p-value:", format(sig_enriched$padj[i], digits = 3), "\n")
    cat("  Target diff genes:\n    ", sig_enriched$overlap_genes[i], "\n")
  }
} else {
  cat("\nNo motifs with p < 0.05\n")
  cat("\nTop 5 motifs by p-value:\n")
  print(enrichment_results %>% 
          select(TF, n_diff_targeted, expected, fold_enrichment, pval) %>%
          head(5))
}

#-----------------------------------------------------------------------------
# Summary statistics
#-----------------------------------------------------------------------------

cat("\n=== SUMMARY ===\n")
cat("Total motifs tested:", nrow(enrichment_results), "\n")
cat("Motifs with p < 0.05:", sum(enrichment_results$pval < 0.05), "\n")
cat("Motifs with p < 0.1:", sum(enrichment_results$pval < 0.1), "\n")
cat("Motifs with padj < 0.1:", sum(enrichment_results$padj < 0.1), "\n")

write.csv(enrichment_results, "./target_gene_enrichment_result_differential_motifs_0.05.csv", row.names = FALSE)
#=============================================================================
# Analyze commonly targeted genes
#=============================================================================

# Get all significant TFs
sig_tfs <- enrichment_results %>% filter(padj < 0.1) %>% pull(motif)

# Count how many TFs target each diff gene
gene_tf_counts <- data.frame(gene = diff_genes_in_background) %>%
  rowwise() %>%
  mutate(
    n_tfs_targeting = sum(sapply(sig_tfs, function(m) {
      gene %in% motif_targets_list[[m]]$targets
    })),
    targeting_tfs = paste(
      sapply(sig_tfs, function(m) {
        if (gene %in% motif_targets_list[[m]]$targets) {
          motif_targets_list[[m]]$TF
        } else NA
      }) %>% na.omit(),
      collapse = ", "
    )
  ) %>%
  ungroup() %>%
  arrange(desc(n_tfs_targeting))

cat("=== Diff genes targeted by multiple significant TFs ===\n")
print(gene_tf_counts %>% filter(n_tfs_targeting > 0))

#=============================================================================
# STEP 5: VISUALIZATION
#=============================================================================

library(ggplot2)
library(ggrepel)

#-----------------------------------------------------------------------------
# 5a. Enrichment bar plot
#-----------------------------------------------------------------------------

plot_data <- enrichment_results %>%
  mutate(
    significant = padj < 0.1,
    TF_label = paste0(TF, " (", n_diff_targeted, ")")
  )

p <- ggplot(plot_data, aes(x = reorder(TF, -log10(padj)), y = -log10(padj))) +
  geom_col(aes(fill = fold_enrichment), width = 0.7) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_fill_gradient(low = "grey70", high = "darkred", name = "Fold\nEnrichment") +
  coord_flip() +
  labs(
    x = "Transcription Factor",
    y = "-log10(p_adj)",
    title = "TF Target Enrichment Among Differential Genes",
    subtitle = "Red line = p_adj < 0.1"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold")
  )

png(filename = "./Bcell_TF_target_enrichment_barplot.png", width = 6, height = 6, units = "in", res = 600)
p
dev.off()

#-----------------------------------------------------------------------------
# 5b. Dot plot (enrichment vs motif accessibility change)
#-----------------------------------------------------------------------------

# Merge enrichment results with motif differential results
plot_data2 <- enrichment_results %>%
  left_join(motif_df %>% select(name, MeanDiff, FDR), 
            by = c("motif" = "name")) %>%
  mutate(significant = padj < 0.1)

p <- ggplot(plot_data2, aes(x = MeanDiff, y = -log10(padj))) +
  geom_point(aes(size = n_diff_targeted, color = fold_enrichment), alpha = 0.8) +
  geom_text_repel(
    data = plot_data2 %>% filter(padj < 0.1),
    aes(label = TF),
    size = 3.5, max.overlaps = 20
  ) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_gradient(low = "blue", high = "red", name = "Fold\nEnrichment") +
  scale_size_continuous(name = "# Diff genes\ntargeted", range = c(3, 12)) +
  labs(
    x = "Motif Accessibility Change (PWH - Control)",
    y = "-log10(Enrichment p_adj)",
    title = "TF Regulatory Activity in CSF B cells",
    subtitle = "X-axis: chromatin accessibility | Y-axis: target gene enrichment"
  ) +
  theme_bw()

png(filename = "./Bcell_TF_target_enrichment_dotplot.png", width = 6, height = 6, units = "in", res = 600)
p
dev.off()


#-----------------------------------------------------------------------------
# 5c. Network visualization
#-----------------------------------------------------------------------------

library(igraph)

# Build edge list from significant TFs
sig_enriched <- enrichment_results %>% filter(padj < 0.1)

edges <- lapply(1:nrow(sig_enriched), function(i) {
  tf <- sig_enriched$TF[i]
  motif <- sig_enriched$motif[i]
  overlap <- intersect(motif_targets_list[[motif]]$targets, diff_genes)
  
  if (length(overlap) > 0) {
    data.frame(from = tf, to = overlap)
  } else NULL
}) %>% bind_rows()

# Create network
library(igraph)

g <- graph_from_data_frame(edges, directed = TRUE)

# Node attributes
V(g)$type <- ifelse(V(g)$name %in% sig_enriched$TF, "TF", "Target")
V(g)$is_hla <- grepl("^HLA", V(g)$name)

V(g)$color <- case_when(
  #V(g)$is_hla ~ "#66C2A5",
  V(g)$type == "TF" ~ "#FC8D62",
  TRUE ~ "#8DA0CB"
)

V(g)$frame.color <- case_when(
  #V(g)$is_hla ~ "#1B9E77",
  V(g)$type == "TF" ~ "#D95F02",
  TRUE ~ "#7570B3"
)

# Size and labels
V(g)$size <- ifelse(V(g)$type == "TF", 12, 6)
V(g)$label.cex <- ifelse(V(g)$type == "TF", 0.8, 0.65)
V(g)$label.color <- "grey20"
V(g)$label.dist <- ifelse(V(g)$type == "TF", 0, 0)

# Edge attributes
E(g)$color <- adjustcolor("grey50", alpha.f = 0.5)
E(g)$width <- 0.8

set.seed(42)
layout <- layout_with_fr(g, niter = 1000)


png(filename = "./Bcell_differential_motif_network.png", width = 6, height = 6, units = "in", res = 600)
# Plot with margins
library(ggraph)
library(tidygraph)
library(ggrepel)

tg <- as_tbl_graph(edges, directed = TRUE) %>%
  mutate(
    type = ifelse(name %in% sig_enriched$TF, "TF", "Target"),
    #type = ifelse(grepl("^HLA", name), "HLA", type)
  )

set.seed(42)
ggraph(tg, layout = "fr") +
  geom_edge_link(
    arrow = arrow(length = unit(1.2, "mm"), type = "closed"),
    end_cap = circle(3, "mm"),
    color = "grey60",
    alpha = 0.5
  ) +
  
  geom_node_point(aes(color = type, size = type)) +
  geom_node_text(
    aes(label = name),
    repel = TRUE,           # <-- automatically avoids overlap
    size = 3,
    max.overlaps = 50,      # allow more labels to show
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "grey50",
    segment.size = 0.2
  ) +
  scale_color_manual(
    values = c("TF" = "#E07A5F", "Target" = "#81A4CD", "HLA" = "#52B788")
  ) +
  scale_size_manual(values = c("TF" = 6, "Target" = 3, "HLA" = 4), guide = "none") +
  theme_void() +
  theme(legend.position = "bottom")

dev.off()
#-----------------------------------------------------------------------------
# 5d. Heatmap of TF-gene connections
#-----------------------------------------------------------------------------
library(ComplexHeatmap)
library(circlize)

sig_tfs <- sig_enriched$motif
target_diff_genes <- unique(unlist(strsplit(sig_enriched$overlap_genes, "; ")))

tf_gene_mat <- matrix(0, 
                      nrow = length(sig_tfs), 
                      ncol = length(target_diff_genes),
                      dimnames = list(sig_enriched$TF, target_diff_genes))

for (i in 1:nrow(sig_enriched)) {
  tf <- sig_enriched$TF[i]
  motif <- sig_enriched$motif[i]
  genes <- unlist(strsplit(sig_enriched$overlap_genes[i], "; "))
  tf_gene_mat[tf, genes] <- 1
}

# Color scheme
col_fun <- colorRamp2(c(0, 1), c("#F5F5F5", "#B2182B"))

# Highlight HLA genes in column names
is_hla <- grepl("^HLA", target_diff_genes)
col_labels <- target_diff_genes
#col_colors <- ifelse(is_hla, "#2D6A4F", "grey20")
col_colors <- "grey20"

# Column annotation for HLA genes
col_anno <- HeatmapAnnotation(
  Type = ifelse(is_hla, "HLA", "Other"),
  col = list(Type = c("HLA" = "#52B788", "Other" = "#81A4CD")),
  annotation_height = unit(3, "mm"),
  show_legend = TRUE,
  annotation_name_side = "left",
  simple_anno_size_adjust = TRUE
)

# Row annotation - number of targets per TF
n_targets <- rowSums(tf_gene_mat)
row_anno <- rowAnnotation(
  `# Targets` = anno_barplot(
    n_targets,
    width = unit(1.5, "cm"),
    gp = gpar(fill = "#E07A5F", col = NA),
    border = FALSE
  ),
  annotation_name_rot = 0
)

# Main heatmap
p <- Heatmap(tf_gene_mat,
             name = "Regulated",
             col = col_fun,
             
             # Clustering
             cluster_rows = TRUE,
             cluster_columns = TRUE,
             clustering_distance_rows = "binary",
             clustering_distance_columns = "binary",
             
             # Appearance
             rect_gp = gpar(col = "grey85", lwd = 0.5),  # cell borders
             border = TRUE,
             
             # Row (TF) settings
             row_names_gp = gpar(fontsize = 10, fontface = "italic"),
             row_title = "Transcription Factors",
             row_title_gp = gpar(fontsize = 12, fontface = "bold"),
             row_dend_width = unit(15, "mm"),
             
             # Column (Gene) settings
             column_names_gp = gpar(fontsize = 8, col = col_colors),
             column_names_rot = 45,
             column_title = "Differential Genes",
             column_title_gp = gpar(fontsize = 12, fontface = "bold"),
             column_dend_height = unit(15, "mm"),
             
             # Annotations
             #top_annotation = col_anno,
             right_annotation = row_anno,
             
             # Legend
             heatmap_legend_param = list(
               title = "Regulation",
               at = c(0, 1),
               labels = c("No", "Yes"),
               legend_height = unit(2, "cm"),
               title_gp = gpar(fontsize = 10, fontface = "bold"),
               labels_gp = gpar(fontsize = 9)
             ),
             
             # Size
             width = ncol(tf_gene_mat) * unit(4, "mm"),
             height = nrow(tf_gene_mat) * unit(8, "mm")
)


png(filename = "./Bcell_differential_motif_heatmap.png", width = 15, height = 7, units = "in", res = 600)
draw(p, 
     padding = unit(c(3, 3, 3, 3), "mm"),
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()