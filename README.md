# Dissertation Analysis Scripts

This repository contains the code and workflows used for my PhD dissertation, organised by chapter.  
It includes all scripts for data preprocessing, analysis, and visualisation for the two main chapters: **glioblastoma (Chapter 3)** and **pond ecosystem dynamics (Chapter 4)**.

---

## Chapter 3 – Glioblastoma Analysis (`Chapter3_tumour`)

This folder is organised into two main components: **single-cell analysis** and **bulk RNA-seq analysis**.

### Single-cell analysis (`phd_chapter_sc`)
- `prep/` – Data cleaning and integration steps  
- `scenic_*` – Regulon enrichment analysis:  
  - One script runs the SCENIC pipeline  
  - One script performs downstream analysis  
- `degs/` – Differentially expressed gene (DEG) analysis  
- `fig_*` – Scripts for generating main figures (e.g., Figure 1, Figure 2)  
- `liana/` – Ligand–receptor interaction analysis, including deeper exploration modules  

### Bulk RNA-seq analysis (`phd_chapter_survival`)
- `ols_R_final.Rmd` – Integration of cohort data and Kaplan–Meier survival analysis (including supplementary figures)  
- `ol_specificity/` – Testing specificity of derived signatures  
- `momf/` – Inferring oligodendrocyte composition in bulk data using single-cell-derived signatures  

---

## Chapter 4 – Pond Ecosystem Dynamics (`Chapter4_ponds`)

This folder contains two main subfolders:

- `models/` – Scripts for modelling; each script is named according to its function  
- `plots/` – Visualisation scripts  
- `input/` – Data input and initial preprocessing (including cleaning steps)  

---

## Notes
- Each chapter is self-contained and corresponds directly to dissertation chapters.  
- Figure-generating scripts are labelled with `fig_*` for traceability.  
- Supplementary analyses are included within relevant subfolders.
