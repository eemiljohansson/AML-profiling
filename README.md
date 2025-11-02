# AML-profiling
This code was generated in the context of the study “Blood proteome profiling using proximity extension assay in patients with acute myeloid leukemia”, were a comprehensive analysis of the plasma proteome was performed on AML patients using a pan-cancer cohort representing the major cancer types and a healthy cohort as controls.

# The study

Acute Myeloid Leukemia (AML) is the most common form of acute leukemias in adults. Given the complexity of the disease, proteomic profiling presents an opportunity to discover early diagnostic biomarkers that could facilitate a shorter and simpler workflow for arriving at the diagnosis. This study was conducted by analyzing 1 463 plasma proteins in 52 AML patients at diagnosis using the Olink Explore 1536 platform. Both differential expression analysis and feature selection by machine learning were applied to find the most significant proteins to distinguish AML from 867 healthy individuals and 1 734 patients of varying cancer types, including different hematological malignancies. The analysis identified proteins with significant altered expression in AML patients, highlighting potential biomarkers.
The results from the study are published in insert journal*: Johansson, E, et al. Blood proteome profiling using proximity extension assay in patients with acute myeloid leukemia. Insert journal here* (2024). insert DOI here*

# Synthetic data

The data repository contains synthetic Olink data (`olink_data.csv` & `ukb_olink_data.csv`) and metadata (`metadata.xlsx` & `ukb_metadata.xlsx`) scrambled for testing the code that comes with the publication. **DO NOT** use the data to replicate this study or for any further work!

# Content

This repository includes the code to generate the results describe above, as well as a synthetic dataset to test the code:
1. `/data`: contains example data to test the code, as well as additional data files to reproduce the analysis. Note that this is not real data and should not be used in any research.
2. `/scripts`: contains all necessary scripts to reproduce the analysis.
3. `/results`: all the plots resulting from the analysis will be stored in this directory. Note that the results are based on synthetic data and should not be interpreted as valid biological resutlts.
4. `AML-profiling.Rproj`: R project file.

# Citation

To use this code in your own research, please cite our code and/or study:
Johansson, E. eemiljohansson/AML-profiling: AML-profiling (insert version here*). Zenodo. insert DOI here* (2024).
Johansson, E, et al. Blood proteome profiling using proximity extension assay in patients with acute myeloid leukemia. Insert journal here* (2024). insert DOI here*
