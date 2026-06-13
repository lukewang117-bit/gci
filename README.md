# gCI: Group-Level Fairness Assessment for Clinical Risk Prediction

This repository contains the analysis code and organized notebooks for the paper on **group-level discrimination metrics** for fairness assessment in predictive models. The paper proposes **gCI** and **gAUC**, which summarize subgroup-specific discrimination relative to the broader population, and contrasts them with traditional **within-group CI/AUC** and pairwise **xCI/xAUC** analyses.

## Paper summary

In many fairness analyses, subgroup performance is reported using within-group discrimination metrics such as the concordance index (CI) or the AUC. Those metrics only evaluate ranking **within** a subgroup. In practice, however, clinical decisions are often made across a shared population, so it is also important to understand how well people in a given group are ranked **relative to everyone else**.

This project introduces:

- **xCI / xAUC**: pairwise cross-group ranking metrics
- **gCI / gAUC**: group-level summary metrics that aggregate all relevant comparisons involving a subgroup

The motivating use case is the **PREVENT** equation for 10-year ASCVD risk prediction.

## Suggested workflow

### 1) Data preparation
Use the R notebooks in `notebooks/01_data_preparation/` to assemble the study cohort, curate variables, and generate analytic artifacts in the Truveta environment.

### 2) Exploratory analysis
Use `notebooks/02_exploratory_analysis/01_prevent_eda.ipynb` for exploratory checks, distributions, and intermediate validation.

### 3) Fairness metric analysis
Use the notebooks in `notebooks/03_metric_analysis/` to compute:

- within-group CI / AUC
- pairwise cross-group xCI / xAUC
- proposed group-level gCI / gAUC

## Data availability

Due to data use restrictions, source EHR data cannot be shared publicly. This repository contains code used for cohort construction, risk score computation, and fairness metric evaluation.



## Citation

```bibtex
@misc{wang_gci_prevent,
  title  = {A simplified metric to streamline between-group fairness assessment for predictive models},
  author = {Wang, Haoyuan and Hong, Chuan and Pencina, Michael and Engelhard, Matthew},
  note   = {Manuscript in revision}
}
```

## Contact

For questions about the code or reproducibility details, please contact Haoyuan Wang (lukewang@sas.upenn.edu).
