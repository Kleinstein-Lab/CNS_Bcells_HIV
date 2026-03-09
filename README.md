# Multimodal single cell analysis reveals persistent memory B cell dysfunction in the CNS of ART-treated people with HIV

This repository contains the full analysis pipeline and code used in the publication: **Multimodal single cell analysis reveals persistent memory B cell dysfunction in the CNS of ART-treated people with HIV**. This work investigates B cell phenotypes and clonality in blood, cerebrospinal fluid (CSF), and choroid plexus tissues of people with HIV (PWH) and people without HIV (PWoH), using single-cell RNA-seq (scRNA-seq), single-nucleus RNA/ATAC-seq (snRNA/snATAC-seq), and B cell receptor sequencing (scBCR-seq).

## Overview

We performed integrative multi-omic profiling of:
- CSF and blood from antiretroviral therapy (ART)-suppressed, neurologically asymptomatic PWH and matched PWoH.
- Postmortem choroid plexus from ART-treated PWH and PWoH.

### Data Types
- **scRNA-seq & scBCR-seq:** 10x Genomics 5′ V(D)J platform.
- **snRNA-seq & snATAC-seq:** 10x Genomics Multiome (ATAC + Gene Expression).
- **Proteomics:** Human Cytokine 96-Plex Discovery Assay (HD96, Eve Technologies)

### Data Processing
- **RNA-seq:**  Cell Ranger v7.0.0, Seurat v5.2.0.
- **BCR-seq:** Cell Ranger VDJ, nf-core/airrflow v4.0 pipeline.
- **ATAC-seq:** Cell Ranger ARC v2.0.2, ArchR v1.0.3.

See the manuscript's Methods for full experimental and computational details.

## Data Availability

Data are available in public repositories:

- **CSF Group 1 (part of BCR-seq data), CSF Group 2 & Choroid Plexus:**
  - SRA Accession: PRJNA1289878
- **CSF Group 1 (Previously Published):**
  - RNA-seq: SRA PRJNA717310 (study PRJNA717310), GEO GSE243905
    - C1_BLD_RNA: SRR14076861
    - C1_CSF_RNA: SRR14076871
    - C2_BLD_RNA: SRR14076860
    - C2_CSF_RNA: SRR14076869
    - C3_BLD_RNA: SRR14076858
    - C3_CSF_RNA: SRR14076868
  - BCR-seq: SRA PRJNA717310
    - PBMC_10 BCR (C1): SRS8582419
    - CSF BCR (C1): SRS8582409
    - PBMC_11 BCR (C2): SRS8582418
    - CSF BCR (C2): SRS8582412
    - PBMC_12 BCR (C3): SRS8582420
    - CSF BCR (C3): SRS8582413
