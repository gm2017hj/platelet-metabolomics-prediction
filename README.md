# Platelet Metabolomics Prediction

R code for **"Donor-wise metabolomic prediction of platelet storage quality: a proof-of-concept study"**.

## Data

Download the three files below from Metabolomics Workbench (study ST003937) and place them in the `data/` folder:

- `MSdata_ST003937_1.txt`
- `MSdata_ST003937_2.txt`
- `ST003937_AN006465.txt`

Access the data at: https://www.metabolomicsworkbench.org (search for ST003937).

## Requirements

- R ≥ 4.5.3
- Packages: `tidyverse`, `caret`, `randomForest`, `xgboost`, `pROC`, `glmnet`, `fastshap`, `shapviz`, `cluster`, `factoextra`, `torch`, `patchwork`

Install all packages with:
```r
install.packages(c("tidyverse", "caret", "randomForest", "xgboost", "pROC", "glmnet", "fastshap", "shapviz", "cluster", "factoextra", "torch", "patchwork"))
If using torch, run torch::install_torch() after installation.

How to run
Place the three data files in the data/ folder.

Open main_analysis.R in RStudio.

Run the entire script. All figures and tables will be generated in the results/ folder.

Or simply run:

r
source("main_analysis.R")
Output
Figures: Fig1–5 (main), FigS1–S3 (supplementary)

Tables: Table1, Supplementary Tables S1–S4

CSV files: AUC values, SHAP importance, donor variability

License
MIT License. See LICENSE for details.

Citation
If you use this code, please cite our paper (forthcoming) and the original data source:

Zhang B, Gao M. Donor-wise metabolomic prediction of platelet storage quality: a proof-of-concept study. Transfusion. (in press)

van Wonderen SF, et al. Blood Adv. 2025;9(20):5164-76.