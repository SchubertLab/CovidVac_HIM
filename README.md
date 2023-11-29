# CovidVac - MVP
In this study, we analyze the effect of multiple vaccinations (n>200) on the human immune response repertoire. This repository contains the single-cell analysis of T cell receptor (TCR), Gene expression, Surface Proteins (Antibody Captured), and dextramer staining across 6 donors (1 frequently vaccinated patients, 5 controls).

## Data
- The raw sequencing data can be downloaded from GEO TODO link (Accession Number).
- The cellranger output can be downloaded from the same repository and should be stored as `./data/20231017/GEX/gex_mixed_run_{1-3}_feature_bc_matrics.h5` and `./data/20231017/VDJ/vdj_mixed_run_{1-3}_filtered_contigs.csv`
- The processed and annotated data can be downloaded from Zenodo (doi: ) and stored as `./data/mvp/02_mvp_annotated_cd4.h5ad` an `./data/mvp/02_mvp_annotated_cd8.h5ad` (entry point notebooks 03)

Additionally, we use the following external resources (not provided by us, store in `./data/scores/`):
- Genes for various T cell scores by Szabo et al. (Nature Communications, 2019) (https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-12464-3/MediaObjects/41467_2019_12464_MOESM7_ESM.xlsx)
- Genes for Cell Scores by Meckiff et al. (Cell, 2020) (https://ars.els-cdn.com/content/image/1-s2.0-S0092867420313076-mmc3.xlsx)

## Reproducibility
To recreate the results of the paper:
```
git clone https://github.com/SchubertLab/CovidVac_MVP.git
cd CovidVac_MVP
conda env create -f requirements.yml
```
Following, the notebooks must be run (ideally in this order). Note, that there are issues with reproducing UMAPs across different machines, even when the same seeds and package versions are used. Results might therefore look slightly different. To fully reproduce the paper results, use the annotated data.

To separate multiple sequencing runs from the cellranger output:
- `01_1_mixed_runs_preprocessing`
- `01_2_mixed_runs_preprocessing`
- `01_2_mixed_runs_preprocessing`

To annotate the data:
- `02_mvp_annotation`

To generate the paper results (if you use the annotated data directly you can start at this point):
- `03A_mvp_visualization_cd4`
- `03B_mvp_visualization_cd8`
- `05_mvp_tables`

## Citation
If you refer to this work, please consider citing the following paper:

```
@article{
}
```
