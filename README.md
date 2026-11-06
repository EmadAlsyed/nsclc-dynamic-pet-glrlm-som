# Dynamic 18F-FDG PET GLRLM features in NSCLC: Self-Organizing Map analysis

Analysis code for the manuscript *"Clustering Dynamic 18F-FDG PET GLRLM Radiomic Features and Their Association with Overall Survival in NSCLC: An Exploratory Self-Organizing Map Study"* (submitted to *Diagnostics*, MDPI).

The script trains a Kohonen self-organizing map (SOM) on 16 dynamic GLRLM features of 35 NSCLC lesions (34 patients), groups the nine map nodes into two super-clusters, assesses map quality and stability, compares the result with PCA, k-means and hierarchical clustering, and tests overall survival (OS) between the super-clusters.

## Data

This is a secondary analysis of public data from:

> Noortman WA, Vriens D, Slump CH, et al. Adding the temporal domain to PET radiomic features. PLoS ONE 2020;15(9):e0239438. https://doi.org/10.1371/journal.pone.0239438

The radiomic features were **not** recomputed. They are the 16 dynamic GLRLM features distributed with that study (S1 Data), based on sixteen 150-s frames between 10 and 50 min after injection.

| File | Content |
|---|---|
| `Onlinedata.csv` | OS (months) and the 16 z-scoreable GLRLM features, 35 rows (one per lesion), extracted from S1 Data |
| `pone.0239438.s002.xlsx` | S1 Data of the parent study (sheet "Clinical data": vital status and lesion ID, e.g. 26A and 26B) |

S1 Data is published under the Creative Commons Attribution (CC BY) license of PLOS ONE. If you prefer not to redistribute it, delete the file from the repository and download it from the article page instead.

## Requirements

- R 4.6.1 (earlier 4.x versions should also work, but the reported numbers come from 4.6.1)
- R packages (installed automatically by the script if missing): `kohonen` 3.0.13, `survival` 3.8-6, `survminer` 0.5.2, `ggplot2`, `RColorBrewer`, `readxl`, `cluster`
- Exact versions of the run that produced the manuscript are in `NSCLC_SOM_Outputs/sessionInfo.txt`.

## How to run

1. Put `NSCLC_SOM_revised_analysis.R`, `Onlinedata.csv` and `pone.0239438.s002.xlsx` in one folder.
2. Set that folder as the working directory and run:

```r
source("NSCLC_SOM_revised_analysis.R")
# or from a terminal:  Rscript NSCLC_SOM_revised_analysis.R
```

The run takes several minutes because of the 200 random initializations, 200 bootstrap SOMs, grid-size runs, 10,000 permutations and 2000 RMST bootstrap resamples.

All random steps use fixed seeds (main map: `set.seed(2024)`), so the results are reproducible with the same package versions. If the script prints the warning "The map differs from the one in the submitted manuscript", your kohonen or R version gives a different map, and the numbers in the manuscript will differ.

## Outputs (`NSCLC_SOM_Outputs/`)

| Output | Manuscript item |
|---|---|
| `Fig2_SOM_training.png` | Figure 2 |
| `Fig3_codebook_heatmap.png`, `Table1_2_node_zscores.csv` | Figure 3, Tables 2 and 3 |
| `Fig4_counts.png` | Figure 4 |
| `Fig5_component_planes.png` | Figure 5 |
| `Fig6_Umatrix.png` | Figure 6 |
| `Fig8_superclusters.png` | Figure 7 |
| `Fig9_KM_primary.png` | Figure 8 |
| `TableS_stability_seeds.csv`, `TableS_stability_bootstrap.csv`, `TableS_grid_size_*.csv`, `TableS_choice_of_k.csv` | Table 4, Figures S2–S5 |
| `Table5_benchmark_methods.csv` | Table 5 |
| `Table4_survival_analyses.csv` | Table 6 |
| `Table3_node_summary_censoring_aware.csv` | Table 7 |
| `Table_KM_steps_by_supercluster.csv`, `Table_KM_summary_by_supercluster.csv` | Supplementary Table S1 |
| `FigS_*.png` | Figures S1–S8 |
| `manuscript_numbers.txt` | Every number quoted in the text |
| `sessionInfo.txt` | R and package versions |

Output file names keep the numbering of the first submission. The mapping to the revised manuscript is given in the table above.

## Script structure

0. Settings; 1. Import, vital status, patient ID; 2. Main SOM; 3. Quality (quantization and topographic error); 4. Stability (seeds, bootstrap, grid sizes); 5. Choice of k and U-matrix; 6. Benchmarks (PCA, k-means, hierarchical); 7. Survival (Kaplan–Meier, log-rank, Cox, Schoenfeld, RMST, permutation); 8. Node summaries; 9. Figures; 10. Numbers file and session information.

## Citation

If you use this code, please cite the manuscript (reference to be added after publication) and the parent study above.

## License

Code: MIT License (see `LICENSE`). Data: as published by Noortman et al. (CC BY 4.0).
