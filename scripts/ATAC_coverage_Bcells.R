# =============================================================================
# B Cell ATAC-seq Coverage Plot Analysis
# =============================================================================

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Hsapiens.v75)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(patchwork)
library(ArchR)
library(stringr)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(harmony)
library(cowplot)
library(dplyr)
library(readxl)
library(colorspace)
library(grDevices)
library(paletteer)

set.seed(1234)
addArchRThreads(threads = 128)
addArchRGenome("hg38")

# -----------------------------------------------------------------------------
# Paths 
# -----------------------------------------------------------------------------
DATA_DIR   <- "."
ATAC_DIR   <- file.path(DATA_DIR, "ATAC/ArchR_combined")
RNA_DIR    <- file.path(DATA_DIR, "RNA")
OUTPUT_DIR <- ATAC_DIR

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
setwd(ATAC_DIR)

proj_atac_csf <- loadArchRProject("./", force = TRUE)
rna_seurat    <- readRDS(file.path(RNA_DIR, "refined_rna_seurat.RDS"))

rna_seurat_Bcell_annotated <- readRDS(
  file.path(RNA_DIR, "6-gex_multi_Bsubclust_annotated.rds")
)

# Transfer B cell subtype annotations to the full Seurat object
rna_seurat$Bcell_subtype <- "Other"
cells_to_update <- Cells(rna_seurat_Bcell_annotated)
rna_seurat$Bcell_subtype[cells_to_update] <- as.character(
  rna_seurat_Bcell_annotated$annotated_clusters_B
)

# -----------------------------------------------------------------------------
# Subset to B cells in ATAC, then further to Resting Memory B cells
# -----------------------------------------------------------------------------
idxPass <- which(proj_atac_csf$rna_celltype == "B cells")
atac_pass_filtered <- proj_atac_csf$cellNames[idxPass]
proj_atac_celltype <- proj_atac_csf[atac_pass_filtered, ]

RM_Bcell <- rna_seurat$new_cell_name[
  which(rna_seurat$Bcell_subtype %in% c("Resting Memory"))
]
proj_atac_celltype <- proj_atac_celltype[
  which(proj_atac_celltype$cellNames %in% RM_Bcell)
]

# -----------------------------------------------------------------------------
# Prepare GeneExpressionMatrix for ArchR
# Match RNA genes to ArchR gene annotations; create placeholder GRanges for
# unmatched genes so that rowRanges can be assigned to the SCE object.
# -----------------------------------------------------------------------------
seRNA <- as.SingleCellExperiment(rna_seurat)
seRNA <- seRNA[, which(seRNA$new_cell_name %in% proj_atac_celltype$cellNames)]

geneAnnotation <- getGeneAnnotation(proj_atac_celltype)
genes     <- geneAnnotation$genes
sce_genes <- rownames(seRNA)
match_idx <- match(sce_genes, genes$symbol)

n_matched <- sum(!is.na(match_idx))
message(
  "Matched ", n_matched, " out of ", length(sce_genes), " genes (",
  round(n_matched / length(sce_genes) * 100, 2), "%)"
)

matched_genes    <- sce_genes[!is.na(match_idx)]
matched_gene_info <- genes[match_idx[!is.na(match_idx)], ]

gr <- GRanges(
  seqnames = seqnames(matched_gene_info),
  ranges   = IRanges(start = start(matched_gene_info),
                     end   = end(matched_gene_info)),
  strand   = strand(matched_gene_info)
)
names(gr) <- matched_genes

unmatched_genes <- sce_genes[is.na(match_idx)]
if (length(unmatched_genes) > 0) {
  message("Creating placeholder entries for ", length(unmatched_genes),
          " unmatched genes")
  placeholder_gr <- GRanges(
    seqnames = rep("chr1", length(unmatched_genes)),
    ranges   = IRanges(start = rep(1, length(unmatched_genes)),
                       end   = rep(1, length(unmatched_genes))),
    strand   = rep("*", length(unmatched_genes))
  )
  names(placeholder_gr) <- unmatched_genes
  all_gr <- c(gr, placeholder_gr)[sce_genes]
} else {
  all_gr <- gr
}

rowRanges(seRNA) <- all_gr
seRNA <- seRNA[matched_genes, ]
colnames(seRNA) <- seRNA$new_cell_name

proj_atac_celltype <- addGeneExpressionMatrix(
  input = proj_atac_celltype, seRNA = seRNA,
  strictMatch = TRUE, force = TRUE
)

# -----------------------------------------------------------------------------
# Iterative LSI on peak matrix
# -----------------------------------------------------------------------------
proj_atac_celltype <- addIterativeLSI(
  ArchRProj   = proj_atac_celltype,
  useMatrix   = "PeakMatrix",
  name        = "IterativeLSI_peak",
  iterations  = 2,
  clusterParams = list(
    resolution  = c(0.2),
    sampleCells = 10000,
    n.start     = 10
  ),
  varFeatures = 25000,
  dimsToUse   = 1:30,
  force       = TRUE
)

# -----------------------------------------------------------------------------
# Peak-to-gene linkage
# -----------------------------------------------------------------------------
proj_atac_celltype <- addPeak2GeneLinks(
  ArchRProj   = proj_atac_celltype,
  reducedDims = "IterativeLSI_peak",
  useMatrix   = "GeneExpressionMatrix",
  k           = 20,        # number of nearest neighbors for KNN smoothing
  maxDist     = 50000      # max peak-gene distance (bp)
)

# Extract loops for visualization (|correlation| > 0.25, FDR < 0.05)
p2g_loops <- getPeak2GeneLinks(
  proj_atac_celltype,
  corCutOff = 0, FDRCutOff = 0.05,
  returnLoops = TRUE
)
p2g_loops$Peak2GeneLinks <- p2g_loops$Peak2GeneLinks[
  which(abs(p2g_loops$Peak2GeneLinks$value) > 0.25),
]

# Extract and filter peak-to-gene data frame
p2geneDF <- metadata(proj_atac_celltype@peakSet)$Peak2GeneLinks
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>%
  {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]

p2geneDF <- as.data.frame(p2geneDF) %>%
  filter(FDR < 0.1, abs(Correlation) > 0.25,
         VarQATAC > 0.25, VarQRNA > 0.25)

# -----------------------------------------------------------------------------
# Relabel HIV status and sample IDs for plotting
# -----------------------------------------------------------------------------
proj_atac_celltype$HIV_status[
  which(proj_atac_celltype$HIV_status == "control")
] <- "PWoH"

sample_map <- read.csv("./sample_map.csv")
proj_atac_celltype$new_sample <- sample_map[proj_atac_celltype$Sample]

# -----------------------------------------------------------------------------
# Coverage plot: PWH vs PWoH
# -----------------------------------------------------------------------------
STATUS_COLORS <- c("PWH" = "#9e9ac8", "PWoH" = "#41ab5d")
GENES_TO_PLOT <- c("HLA-DQA1", "TGFBR2", "NFKBIA", "UBE2K")

p <- plotBrowserTrack(
  ArchRProj    = proj_atac_celltype,
  groupBy      = "HIV_status",
  geneSymbol   = GENES_TO_PLOT,
  upstream     = 20000,
  downstream   = 20000,
  title        = "",
  minCells     = 3,
  pal          = unname(STATUS_COLORS),
  loops        = p2g_loops,
  facetbaseSize = 10
)

output_file <- file.path(OUTPUT_DIR, "Bcell_ATAC_coverage_plot_HIV_status.pdf")
pdf(file = output_file, width = 10, height = 8)
for (i in seq_along(p)) {
  grid::grid.newpage()
  grid::grid.draw(p[[i]])
}
dev.off()