# AML-profiling
This code was written to make the various analyses made in the study “Blood Proteome Profiling Using Proximity Extension Assay in Patients with Acute Myeloid Leukemia”.

# The study

Acute Myeloid Leukemia (AML) is the most common form of acute leukemias in adults with low survival rate. Given the complexity of the disease, plasma proteomic profiling presents an opportunity to discover diagnostic biomarkers that could facilitate a shorter and simpler workflow for arriving at the diagnosis. This study was conducted by analyzing 1 158 unique proteins, measured in blood plasma, in 52 AML patients at diagnosis using the Olink Explore 1536 platform. Both differential expression analysis and feature selection were applied to find the most significant proteins to distinguish AML from 867 healthy individuals and 1 734 patients of varying cancer types, including different hematological malignancies. The analysis identified proteins with significant altered expression in AML patients, highlighting potential biomarkers.
The results from the study are published in insert journal*: Johansson, E, et al. Blood Proteome Profiling Using Proximity Extension Assay in Patients with Acute Myeloid Leukemia. Insert journal here* (2026). insert DOI here*

# Content

This repository includes the code to generate the results describe above:
1. `/scripts`: contains all necessary scripts to reproduce the analysis.
2. `AML-profiling.Rproj`: R project file.

## Directory structure

```bash
📁 AML-profiling/
├── 📂 scripts/                # Core analysis scripts (Quarto)
│   ├── 01_exploratory_analysis.qmd            # Exploration of the cohorts used in the project
│   ├── 02_differential_analysis.qmd           # Differential analysis by limma
│   ├── 03_enrichment_analysis.qmd             # Enrichment analysis by GSEA
│   ├── 04_feature_selection.qmd               # Feature selection by LASSO
│   ├── 05_analysis_summary.qmd                # Summarizing results
│   ├── 06_ukb_summary                         # Overview of the UK-Biobank data
│   └── 📂 functions/         # Core function scripts (Quarto)
│       └── custom_plots.R             # Theme and palettes
├── 📄 README.md               # Overview of the project
└── 📜 LICENSE                 # License file
```
# Citation

Johansson, E, et al. Blood Proteome Profiling using Proximity Extension Assay in Patients with Acute Myeloid Leukemia. Insert journal here* (2026). insert DOI here*
