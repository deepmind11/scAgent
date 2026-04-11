# Single-Cell RNA-Seq on 10x Chromium: Experimental Paradigms, Analysis Workflows, and Parameter Selection

**Research Brief — April 2026**

---

## Table of Contents

1. [Experimental Paradigms](#1-experimental-paradigms)
2. [Experimental Designs Within Each Paradigm](#2-experimental-designs-within-each-paradigm)
3. [Canonical Analysis Workflow](#3-canonical-analysis-workflow)
4. [Tools Landscape](#4-tools-landscape)
5. [Parameter Selection Guidance](#5-parameter-selection-guidance)
6. [Automated Cell Type Annotation](#6-automated-cell-type-annotation)
7. [Agentic Workflow Implications](#7-agentic-workflow-implications)
8. [Sources](#8-sources)

---

## 1. Experimental Paradigms

These are the **major research paradigms** in which 10x Chromium scRNA-seq is routinely deployed. Each paradigm has distinct analysis needs.

### Paradigm 1: Cell Atlas / Tissue Characterization
**Goal:** Catalog all cell types in a tissue or organ.
**Key output:** Cell type taxonomy, marker gene lists, relative proportions.
**Examples:** Human Cell Atlas, Tabula Muris, Tabula Sapiens.

### Paradigm 2: Disease vs. Healthy Comparison
**Goal:** Identify disease-associated cell states, composition changes, and differentially expressed genes.
**Key output:** Differential abundance, cross-condition DE, disease-specific cell states.
**Examples:** Tumor microenvironment profiling, autoimmune tissue biopsies, COVID-19 PBMC studies.

### Paradigm 3: Developmental / Differentiation Trajectory
**Goal:** Order cells along a pseudotemporal axis to map differentiation or developmental processes.
**Key output:** Trajectory graphs, pseudotime ordering, branch-specific gene dynamics.
**Examples:** Embryonic organogenesis, hematopoietic differentiation, neural crest development.

### Paradigm 4: Perturbation Screens (Perturb-seq / CRISPR Screening)
**Goal:** Link genetic perturbations (CRISPR KO/KD/activation) to transcriptomic phenotypes at single-cell resolution.
**Key output:** Perturbation–gene regulatory networks, gene program identification, epistatic interactions.
**Examples:** Genome-scale Perturb-seq in T cells, in vivo brain perturbation screens, drug mechanism studies.
**10x product:** Chromium Single Cell CRISPR Screening.

### Paradigm 5: Immune Repertoire / V(D)J Profiling
**Goal:** Pair full-length TCR/BCR clonotype sequences with gene expression in the same cell.
**Key output:** Clonotype frequencies, clonal expansion analysis, phenotype-to-clonotype mapping.
**Examples:** Vaccine response tracking, tumor-infiltrating lymphocyte characterization, autoimmune repertoire studies.
**10x product:** Chromium Single Cell 5' + V(D)J.

### Paradigm 6: Multi-Modal / Surface Protein + Transcriptome (CITE-seq style)
**Goal:** Simultaneously measure surface protein markers (via oligo-tagged antibodies) and whole transcriptome.
**Key output:** Protein-level cell classification, weighted nearest neighbor (WNN) integration, protein-RNA correlations.
**Examples:** PBMC immunophenotyping with antibody panels, tumor antigen characterization.
**10x product:** Chromium Single Cell Gene Expression with Feature Barcode technology.

### Paradigm 7: Multiome (ATAC + Gene Expression)
**Goal:** Simultaneously capture chromatin accessibility and gene expression from the same cell.
**Key output:** Linked regulatory elements, gene regulatory networks, epigenetic cell state classification.
**Examples:** Developmental enhancer dynamics, cancer epigenetic heterogeneity, cell fate commitment.
**10x product:** Chromium Single Cell Multiome ATAC + Gene Expression.

### Paradigm 8: Spatial Integration (Chromium + Visium/Xenium)
**Goal:** Combine dissociated single-cell data with spatial transcriptomics to recover spatial context.
**Key output:** Cell type deconvolution of spatial spots, spatial cell communication maps, spatially-resolved gene programs.
**Examples:** Tumor spatial architecture, tissue microenvironment mapping, organ development.

### Paradigm 9: Temporal / Longitudinal Studies
**Goal:** Track transcriptomic changes across time points (treatment, disease progression, aging).
**Key output:** Dynamic cell state transitions, treatment response trajectories, cell composition shifts over time.
**Examples:** Drug response time courses, aging atlases, infection dynamics.

---

## 2. Experimental Designs Within Each Paradigm

### 2.1 Cell Atlas / Tissue Characterization

| Design | Description |
|--------|-------------|
| **Single-tissue deep profiling** | One tissue, high cell count (50K–500K), aim for rare cell discovery |
| **Multi-donor atlas** | Same tissue from multiple individuals to capture population-level variation |
| **Cross-tissue comparison** | Multiple tissues from the same individual for systemic characterization |
| **Species comparison** | Same tissue from human vs. mouse for evolutionary conservation |
| **Developmental stage series** | Same tissue at multiple developmental stages (e.g., fetal vs. adult) |
| **Dissociation protocol comparison** | Test enzymatic digestion methods to check dissociation artifacts |

### 2.2 Disease vs. Healthy Comparison

| Design | Description |
|--------|-------------|
| **Case-control** | Matched diseased vs. healthy donors, same tissue |
| **Treatment response** | Pre/post treatment from same patient (paired design) |
| **Disease subtype stratification** | Multiple disease subtypes profiled to find subtype-specific signatures |
| **Graded severity** | Samples along disease severity spectrum (mild/moderate/severe) |
| **Tumor vs. adjacent normal** | Intra-patient comparison of tumor and adjacent normal tissue |
| **Multi-region tumor** | Multiple biopsies from different regions of same tumor |

### 2.3 Developmental / Differentiation Trajectory

| Design | Description |
|--------|-------------|
| **Time-course sampling** | Samples at defined time points during differentiation |
| **In vitro organoid differentiation** | Organoids sampled at progressive stages |
| **Lineage tracing + scRNA-seq** | Cells carrying heritable barcodes profiled by scRNA-seq |
| **Embryonic staging** | Whole embryos or tissues at multiple embryonic ages |
| **Reprogramming dynamics** | iPSC reprogramming sampled over time |

### 2.4 Perturbation Screens

| Design | Description |
|--------|-------------|
| **Pooled CRISPR KO screen** | Thousands of guides pooled, one guide per cell, read out by scRNA-seq |
| **CRISPRi/CRISPRa screen** | Knockdown or activation screens with transcriptomic readout |
| **Combinatorial perturbation** | 2-guide combinations for epistatic interactions |
| **Drug perturbation** | Cells treated with different drugs/doses, profiled individually |
| **Perturb-seq + protein** | CRISPR + CITE-seq for multi-modal perturbation readout |
| **In vivo Perturb-seq** | Perturbation delivered to animals, tissues harvested and profiled |

### 2.5 Immune Repertoire / V(D)J

| Design | Description |
|--------|-------------|
| **Pre/post vaccination** | Track clonal expansion after vaccination |
| **Tumor-infiltrating vs. peripheral** | Compare TIL and PBMC repertoires |
| **Longitudinal infection monitoring** | Serial blood draws during infection |
| **Autoimmune clonotype tracking** | Identify autoreactive clones in affected tissue |
| **Antigen-specific enrichment** | Tetramer-sorted cells + V(D)J + GEX |

### 2.6 Multi-Modal (CITE-seq style)

| Design | Description |
|--------|-------------|
| **Deep immunophenotyping** | Large antibody panels (100+ proteins) + transcriptome |
| **Surface marker discovery** | Screen novel surface proteins for cell identity |
| **Signaling state profiling** | Phospho-protein proxies + transcriptome |
| **Protein-guided cell classification** | Use protein data to define cell types, validate with RNA |
| **Cross-sample multiplexed CITE-seq** | Hashtag antibodies for multiplexing + protein + RNA |

### 2.7 Multiome (ATAC + GEX)

| Design | Description |
|--------|-------------|
| **Enhancer-gene linking** | Identify cis-regulatory elements controlling specific genes |
| **Epigenetic cell type definition** | Define cell types by chromatin + transcriptome |
| **TF activity inference** | Infer transcription factor activity from motif accessibility |
| **Disease epigenetic rewiring** | Compare chromatin landscape in disease vs. healthy |
| **Developmental chromatin dynamics** | Track opening/closing of regulatory elements over differentiation |

### 2.8 Spatial Integration

| Design | Description |
|--------|-------------|
| **Cell type deconvolution** | Use Chromium scRNA-seq as reference for Visium spot decomposition |
| **Spatial niche analysis** | Map cell neighborhoods and their communication |
| **Ligand–receptor spatial mapping** | Combine scRNA-seq communication inference with spatial data |
| **Validation design** | Use Chromium for discovery, Xenium for in situ validation |
| **Multi-modal spatial** | Chromium Multiome + Visium/Xenium for integrated spatial-epigenomic maps |

### 2.9 Temporal / Longitudinal

| Design | Description |
|--------|-------------|
| **Dense time-course** | Many closely spaced time points for continuous dynamics |
| **Sparse longitudinal** | Widely spaced time points (months) for long-term changes |
| **Multi-condition × time** | Treatment vs. control sampled at each time point |
| **Biological replicate time series** | Multiple individuals sampled at each time point |
| **Pulse-chase + scRNA-seq** | Metabolic labeling (e.g., scSLAM-seq on Chromium) + time |

---

## 3. Canonical Analysis Workflow

The following is the standard analysis pipeline for 10x Chromium scRNA-seq (gene expression), ordered from raw data to biological interpretation.

### Step 0: Raw Data Processing (Cell Ranger)
- **Tool:** Cell Ranger `multi` or `count` (10x proprietary)
- **Input:** FASTQ files
- **Output:** Feature-barcode count matrix (filtered and raw), web_summary.html, .cloupe file

### Step 1: Quality Control (QC) and Cell Filtering
- **Metrics:** UMI counts per cell, genes per cell, % mitochondrial reads, % ribosomal reads
- **Tools:** Scanpy (`sc.pp.filter_cells`, `sc.pp.filter_genes`), Seurat (`PercentageFeatureSet`)
- **Optional advanced:** SoupX / CellBender (ambient RNA removal), scDblFinder / DoubletFinder (doublet detection)
- **Typical thresholds (PBMC):** >200 genes, <5000 genes, <10% mito (dataset-dependent)

### Step 2: Normalization
- **Methods:**
  - **SCTransform** (Seurat): Recommended by 10x for Chromium data
  - **scran pooling** (Bioconductor): Pool-based size factors
  - **Log-normalization + scaling** (Scanpy default)

### Step 3: Feature Selection (Highly Variable Genes)
- **Typical count:** 2,000–5,000 HVGs

### Step 4: Dimensionality Reduction — PCA
- **Typical parameters:** Top 30–50 PCs retained

### Step 5: Batch Effect Correction / Data Integration
- **Tools:** Harmony (recommended by 10x), scVI, Scanorama, BBKNN

### Step 6: Neighbor Graph Construction
- **Parameters:** `n_neighbors` (10–30), `n_pcs` (20–50), `metric` (cosine/euclidean)

### Step 7: Clustering
- **Algorithm:** Leiden (preferred) or Louvain
- **Key parameter:** `resolution` (see §5 for guidance)

### Step 8: Visualization — UMAP / t-SNE

### Step 9: Differential Expression Analysis
- **Cluster markers:** Wilcoxon (Scanpy) / FindAllMarkers (Seurat)
- **Cross-condition (pseudobulk):** DESeq2, edgeR

### Step 10: Cell Type Annotation
- **Automated:** CellTypist, Azimuth, SingleR, Cell Ranger (v9+)
- **Manual:** Marker gene comparison against literature/databases

### Step 11+: Advanced / Downstream Analyses
- Trajectory inference (Monocle 3, Slingshot, PAGA)
- RNA velocity (scVelo)
- Cell–cell communication (CellChat, CellPhoneDB, NicheNet)
- Gene set enrichment (GSEA, clusterProfiler, fgsea)
- Transcription factor analysis (pySCENIC, SCENIC+)
- Compositional analysis (scCODA, miloR)

---

## 4. Tools Landscape

### Primary Frameworks

| Framework | Language | Strengths | Weaknesses |
|-----------|----------|-----------|------------|
| **Scanpy** | Python | Scalable, GPU support, AnnData ecosystem | Fewer statistical methods for DE |
| **Seurat** | R | Mature, extensive vignettes, large community | Memory-heavy for >500K cells |
| **Cell Ranger** | CLI | Official 10x pipeline, QC metrics | Black-box, limited customization |
| **Loupe Browser** | Desktop | Interactive, no-code, directly reads .cloupe | Limited to 10x analysis |

---

## 5. Parameter Selection Guidance

### 5.1 Leiden Clustering Resolution

**Key findings from Sciaraffa et al. (2025), Frontiers in Bioinformatics:**
- Higher resolution (up to 2.0) consistently improved clustering accuracy across all datasets
- Best configurations: resolution = 2, n_neighbors = 10, UMAP method, cosine metric
- Worst configurations: resolution = 0.5, n_neighbors = 30, Gauss method
- **WC_dispersion and Banfield-Raftery index** can serve as unsupervised proxies for clustering quality

**Automated resolution selection tools:**
- **chooseR** (Patterson-Cross et al., 2021): Subsampling-based robustness metric
- **callback** (knockoff calibration): Generates knockoff features to avoid over-clustering
- **MultiK** (Liu et al., 2021): Consensus across K-means runs
- **scSHC** (Grabski et al.): Statistical hypothesis testing for cluster significance

### 5.2 Number of Nearest Neighbors (k)
- k = 10–15 for fine-grained, k = 20–30 for coarser analysis

### 5.3 Number of Principal Components
- Typical range: 20–50, dataset-dependent

### 5.4 Distance Metric
- Cosine preferred for default; Euclidean works well with UMAP method + high resolution

### 5.5 QC Thresholds
- Must be adjusted per tissue, sample quality, and species

---

## 6. Automated Cell Type Annotation

| Tool | Approach | Language | Strengths |
|------|----------|----------|-----------|
| **Cell Ranger (v9+)** | Reference-based | Built-in | No setup; human/mouse |
| **CellTypist** | Supervised logistic regression | Python | Fast; 100+ pre-trained models |
| **Azimuth** | Reference mapping (Seurat) | R/Web | Web interface; PBMC reference very good |
| **SingleR** | Correlation-based | R | Simple; many references |
| **scAgent / CellTypeAgent** | LLM-based | Python | Novel cell type discovery (emerging) |

---

## 7. Agentic Workflow Implications

**scBench** (Workman et al., 2026, arXiv:2602.09063) evaluated 8 frontier LLMs on 394 verifiable scRNA-seq problems:
- Best model accuracy: 52.8% (Claude Opus 4.6)
- Normalization tasks: ~82–84% accuracy
- QC tasks: ~61–64% accuracy
- Cell typing: ~48% accuracy
- Differential expression: ~41% accuracy
- Platform matters: 40+ pp accuracy drop on unfamiliar platforms

---

## 8. Sources

### Primary Sources Consulted

- [Best Practices for Analysis of 10x Genomics Single Cell RNA-seq Data](https://www.10xgenomics.com/analysis-guides/best-practices-analysis-10x-single-cell-rnaseq-data)
- [Chromium Single Cell Platform](https://www.10xgenomics.com/single-cell)
- Sciaraffa et al. (2025). *Frontiers in Bioinformatics* 5:1562410. https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2025.1562410/full
- Patterson-Cross et al. (2021). *BMC Bioinformatics* 22:39.
- Xu et al. (2025). "scCluBench." arXiv:2512.02471.
- Workman et al. (2026). "scBench." arXiv:2602.09063.
- CellAgent (2024). arXiv:2407.09811.
- Luecken & Theis (2019). *Molecular Systems Biology* 15:e8746.
- Heumos et al. (2023). *Nature Reviews Genetics* 24:550–572.
- Füllgrabe et al. (2020). *Nature Biotechnology* 38:1384–1386.
