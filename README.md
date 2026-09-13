# Comparative WGCNA: Bud vs. Flower Tissue Sex-Biased Modules

An R-based systems biology workflow for weighted gene co-expression network analysis (WGCNA) comparing two tissue developmental stages (**Bud** and **Flower/Flw**). This pipeline constructs co-expression networks, evaluates scale-free topology, correlates module eigengenes with sex traits, and intersects module members with differentially expressed gene (DEG) datasets.

## Computational Workflow

1. **Data Preprocessing & QC:**
   * Transposes raw FPKM expression matrices to standard WGCNA format (samples $\times$ genes).
   * Runs `goodSamplesGenes` checks to identify sparse/missing expression values.
   * Performs average-linkage hierarchical clustering (`hclust`) to evaluate sample outliers.

2. **Network Topology & Module Detection:**
   * Fits scale-free topology ($R^2 > 0.8$) across a range of soft-thresholding powers ($\beta$).
   * Executes blockwise network construction (`blockwiseModules`) with dynamic tree cutting (`minModuleSize = 30`, `mergeCutHeight = 0.25`).

3. **Trait Correlation & Visualization:**
   * Formats categorical sample metadata (Sex: Female = 0, Male = 1).
   * Calculates Module Eigengenes (MEs) and Student asymptotic p-values for module-trait relationships.
   * Exports publication-ready PDF visual plots (dendrogram color bars and module-trait heatmaps).

4. **DEG Integration & Export:**
   * Extracts target gene lists from significant modules (`midnightblue`, `purple`, `pink`, `green`, `brown`).
   * Intersects module members with statistically significant DEGs ($q \le 0.05$).
   * Saves workspace snapshots (`.RData`) prepped for downstream Cytoscape network visualization.

## Tech Stack & Packages
* **Language:** R
* **Core Libraries:** `WGCNA`, `tidyverse`

---
*Developed as part of research in the Guerrero Lab at North Carolina State University.*
