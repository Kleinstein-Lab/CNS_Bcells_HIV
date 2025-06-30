library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(patchwork)
library(readr)
library(hdf5r)
library(ArchR)
library(stringr)
library(rtracklayer)
library(foreach)
library(doParallel)
library(GenomicRanges)
library(ggbio)
library(GenomicFeatures)
library(harmony)
library(VennDiagram)
library(dplyr)
library(UpSetR)
set.seed(1234)

addArchRThreads(threads = 128) 
addArchRGenome("hg38")


# ***********************
# Create ArchR Project #
# ***********************

setwd("/banach2/chang/CSF_HIV_Shelli/DC_multiome_paper/ATAC")

sample_dir_all <- list.dirs("/banach2/chang/CSF_HIV_Shelli/data/paper_data", full.names = TRUE, recursive = FALSE)
sample_name <- list.dirs("/banach2/chang/CSF_HIV_Shelli/data/paper_data", full.names = FALSE, recursive = FALSE)
h5_file_path <- paste0(sample_dir_all,"/","filtered_feature_bc_matrix.h5")
input_file_path <- paste0(sample_dir_all,"/","atac_fragments.tsv.gz")
print(input_file_path)

ArrowFiles <- createArrowFiles(
  inputFiles = input_file_path,
  sampleNames = sample_name,
  #validBarcodes = rna_no_doublet$cell_id_unique,
  minTSS = 0, # will do the quality control in the further steps
  minFrags = 0,
  maxFrags = 10e7,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = TRUE,
  excludeChr = c("chrM", "chrY"),
)

# If used 3077 before, run the following command to remove it
# ArrowFiles <- ArrowFiles[c(1:6,8:14)]

proj_atac <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "/banach2/chang/CSF_HIV_Shelli/DC_multiome_paper/ATAC/ArchR_combined",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

saveArchRProject(ArchRProj = proj_atac, outputDirectory = "/banach2/chang/CSF_HIV_Shelli/DC_multiome_paper/ATAC/ArchR_combined/", load = FALSE)
save(ArrowFiles, file = "/banach2/chang/CSF_HIV_Shelli/DC_multiome_paper/ATAC/ArchR_combined/ArrowFiles.RData")


# ******************
# Quality Control #
# ******************

# *******************************************************
##### Load cellranger cell names from RNA and filter ATAC
sample_dir_all <- list.dirs("/banach2/chang/CSF_HIV_Shelli/data/paper_data", full.names = TRUE, recursive = FALSE)
sample_name <- list.dirs("/banach2/chang/CSF_HIV_Shelli/data/paper_data", full.names = FALSE, recursive = FALSE)
h5_file_path <- paste0(sample_dir_all,"/","filtered_feature_bc_matrix.h5")
input_file_path <- paste0(sample_dir_all,"/","atac_fragments.tsv.gz")

sample_Seurat_list <- list()
for(i in 1:length(sample_name)){
  data <- Read10X_h5(filename = h5_file_path[i])
  sample_Seurat_list[sample_name[i]] <- CreateSeuratObject(counts = data$`Gene Expression`, project = sample_name[i],assay = "RNA")
}

# Add sample names to the cell names
for(i in names(sample_Seurat_list)){
  sample_Seurat_list[[i]] <- RenameCells(sample_Seurat_list[[i]], add.cell.id = i)
  #sample_Seurat_list[[i]] <- RenameCells(sample_Seurat_list[[i]], new.names =gsub("_","#",Cells(sample_Seurat_list[[i]])))
}

for(i in names(sample_Seurat_list)){
  sample_Seurat_list[[i]] <- RenameCells(sample_Seurat_list[[i]], new.names =gsub("_","#",Cells(sample_Seurat_list[[i]])))
}

# Create the list of cell names of all samples
cellnames_before_filter <- c()
for(i in names(sample_Seurat_list)){
  cellnames_before_filter <- c(cellnames_before_filter,Cells(sample_Seurat_list[[i]]))
}

idxPass <- which(proj_atac$cellNames %in% cellnames_before_filter)
cellsPass <- proj_atac$cellNames[idxPass]
proj_atac <- proj_atac[cellsPass, ]


idxPass <- which(log10(proj_atac$nFrags) >= 3 & proj_atac$TSSEnrichment >= 5)
atac_pass_filtered <- proj_atac$cellNames[idxPass]

proj_atac <- proj_atac[atac_pass_filtered, ]

# *******************************************************
#### Load RNA(before removing doublets) and ATAC object (after all QCs)
proj_atac <- loadArchRProject("./")
rna_seurat <- readRDS("/banach2/chang/CSF_HIV_Shelli/DC_multiome_paper/RNA/refined_rna_seurat.RDS")

# ********************************************************
#### Add RNA metadata (cell type, HIV status, sex, race, ethnicity)
rna_celltypes <- as.character(rna_seurat$annotated_clusters_all)
rna_cell_names <- rna_seurat$new_cell_name
cell_type_map <- setNames(rna_celltypes, rna_cell_names)
atac_cell_names <- proj_atac$cellNames
atac_celltypes <- ifelse(atac_cell_names %in% names(cell_type_map),
                         cell_type_map[atac_cell_names],
                         "unknown")
proj_atac$rna_celltype <- atac_celltypes


sample_to_condition <- unique(rna_seurat@meta.data[, c("sample_final", "sex","age","ethnicity","race","HIV_status","tissue","subject_id")])
rownames(sample_to_condition) <- sample_to_condition$sample_final
proj_atac$sex <- sample_to_condition[proj_atac$Sample,]$sex
proj_atac$age <- sample_to_condition[proj_atac$Sample,]$age
proj_atac$ethnicity <- sample_to_condition[proj_atac$Sample,]$ethnicity
proj_atac$race <- sample_to_condition[proj_atac$Sample,]$race
proj_atac$HIV_status <- sample_to_condition[proj_atac$Sample,]$HIV_status
proj_atac$tissue <- sample_to_condition[proj_atac$Sample,]$tissue
proj_atac$subject_id <- sample_to_condition[proj_atac$Sample,]$subject_id
proj_atac$harc_id <- substr(proj_atac$Sample,3,6)

# Add myeloid information
rna_celltypes <- as.character(rna_seurat$celltype_res1)
rna_cell_names <- rna_seurat$new_cell_name
cell_type_map <- setNames(rna_celltypes, rna_cell_names)
atac_cell_names <- proj_atac$cellNames
atac_celltypes <- ifelse(atac_cell_names %in% names(cell_type_map),
                         cell_type_map[atac_cell_names],
                         "unknown")
proj_atac$finer_myeloid_celltype <- atac_celltypes

rna_celltypes <- as.character(rna_seurat$whether_myeloid)
rna_cell_names <- rna_seurat$new_cell_name
cell_type_map <- setNames(rna_celltypes, rna_cell_names)
atac_cell_names <- proj_atac$cellNames
atac_celltypes <- ifelse(atac_cell_names %in% names(cell_type_map),
                         cell_type_map[atac_cell_names],
                         "unknown")
proj_atac$whether_myeloid <- atac_celltypes

proj_atac$celltype_final_res1 <- proj_atac$finer_myeloid_celltype
proj_atac$celltype_final_res1[which(proj_atac$celltype_final_res1 %in% c("MG CCL2+","MG SPP1+","MG transitional"))] <- "Microglia-like"
proj_atac$celltype_final_res1[which(proj_atac$celltype_final_res1 %in% c("DC1","DC2","DC3/Mono-derived","DC5/asDC"))] <- "cDC"

proj_atac$celltype_final_res2 <- proj_atac$finer_myeloid_celltype
proj_atac$celltype_final_res2[which(proj_atac$celltype_final_res2 %in% c("MG CCL2+"))] <- "Microglia-like CCL2+"
proj_atac$celltype_final_res2[which(proj_atac$celltype_final_res2 %in% c("MG SPP1+","MG transitional"))] <- "Microglia-like CCL2-"

saveArchRProject(proj_atac)

# Filter the atac doublets
proj_atac <- filterDoublets(proj_atac,filterRatio=1)

# filter the cells without rna cell type label
idxPass <- which(!proj_atac$rna_celltype == "unknown")
atac_pass_filtered <- proj_atac$cellNames[idxPass]
proj_atac <- proj_atac[atac_pass_filtered, ]

saveArchRProject(proj_atac)

# **************
# Peak Calling #
# **************

proj_atac <- addGroupCoverages(
  ArchRProj = proj_atac,
  groupBy = "celltype_final_res1",
  useLabels = TRUE,
  minCells = 40,
  maxCells = 500,
  maxFragments = 25 * 10^6,
  minReplicates = 2,
  maxReplicates = 5,
  sampleRatio = 0.8,
  kmerLength = 6,
  threads = getArchRThreads(),
  returnGroups = FALSE,
  parallelParam = NULL,
  force = TRUE,
  verbose = TRUE,
  logFile = createLogFile("addGroupCoverages")
)

pathToMacs2 <- "/home/chang/miniconda3/bin/macs2"
proj_atac <- addReproduciblePeakSet(
  ArchRProj = proj_atac,
  groupBy = "celltype_final_res1",
  minCells = 25,
  pathToMacs2 = pathToMacs2,
  excludeChr = c("chrM", "chrY")
)

proj_atac <- addPeakMatrix(proj_atac)
getAvailableMatrices(proj_atac)

saveArchRProject(proj_atac)

proj_atac <- addIterativeLSI(
  ArchRProj = proj_atac,
  useMatrix = "PeakMatrix", 
  name = "IterativeLSI_peak", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  force = TRUE
)

print("finish LSI")
saveArchRProject(proj_atac)

proj_atac <- addHarmony(
  ArchRProj = proj_atac,
  reducedDims = "IterativeLSI_peak",
  name = "Harmony_peak",
  groupBy = "Sample",
  force = TRUE
)


print("finish Harmony")
saveArchRProject(proj_atac)

proj_atac <- addUMAP(
  ArchRProj = proj_atac, 
  reducedDims = "Harmony_peak", 
  name = "UMAP_peak", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force=TRUE
)

saveArchRProject(proj_atac)

# ******************
# B cell analysis #
# ******************
proj_atac_csf <- proj_atac[which(proj_atac$tissue == "CSF")]

idxPass <- which(proj_atac_csf$celltype_final_res1 == "B cells")
atac_pass_filtered <- proj_atac_csf$cellNames[idxPass]
proj_atac_celltype <- proj_atac_csf[atac_pass_filtered, ]

# Add peak to gene links
seRNA <- as.SingleCellExperiment(rna_seurat)
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
  # Create placeholder GRanges for unmatched genes
  # Using chromosome 1 position 1 as placeholder (these won't be used by ArchR)
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
proj_atac_celltype <- addGeneExpressionMatrix(input = proj_atac_celltype, seRNA = seRNA, strictMatch = TRUE, force = TRUE)

proj_atac_celltype <- addIterativeLSI(
  ArchRProj = proj_atac_celltype,
  useMatrix = "PeakMatrix", 
  name = "IterativeLSI_peak", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  force = TRUE
)

proj_atac_celltype <- addPeak2GeneLinks(
  ArchRProj = proj_atac_celltype,
  reducedDims = "IterativeLSI_peak",
  useMatrix = "GeneExpressionMatrix",
  k = 20,
  maxDist = 50000
)

p2geneDF <- metadata(proj_atac_celltype@peakSet)$Peak2GeneLinks
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]


p2geneDF <-  as.data.frame(p2geneDF)
p2geneDF <- p2geneDF %>% filter(FDR < 0.1)
p2geneDF <- p2geneDF %>% filter(abs(Correlation) > 0.25)
p2geneDF <- p2geneDF %>% filter(VarQATAC > 0.25)
p2geneDF <- p2geneDF %>% filter(VarQRNA > 0.25)

gene_name <- Genes_Bcells_ATAC_explore$gene[Genes_Bcells_ATAC_explore$gene %in% p2geneDF$geneName]
