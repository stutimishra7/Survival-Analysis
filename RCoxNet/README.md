# RCoxNet

Code for the paper:

**RCoxNet: A deep learning framework integrating random walk with restart, mutation, and clinical data for cancer survival prediction**  
Stuti Kumari, Sakshi Gujral, Abhishek Halder, Smruti Panda, Bernadette Mathew, Prashant Gupta, Ralf Herwig, Gaurav Ahuja, Debarka Sengupta  
*Journal of Computational Biology*, 2025  
DOI: https://doi.org/10.1101/2024.09.17.613428

---

## Overview

RCoxNet propagates somatic mutation profiles across a ConsensusPathDB protein–protein interaction (PPI) network using Random Walk with Restart (RWR). The network-smoothed scores are combined with clinical covariates (age, TMB, MSI) and passed through a deep Cox proportional hazards model to predict overall survival.

We evaluated on four TCGA cohorts — BRCA, LUNG, GBM, and OV — using 20 independent random splits per cohort.

---

## Repository layout

```
RCoxNet/
├── pipeline_rcoxnet.ipynb    # step-by-step walkthrough (Steps 1-7)
├── km_curves.ipynb           # Kaplan-Meier survival curves (Figure 2)
├── best_params.json          # best hyperparameters from our grid search
├── run.py                    # command-line entry point
├── requirements.txt
├── rcoxnet/                  # core library
│   ├── config.py             # constants, cancer configs, hyperparameter grid
│   ├── model.py              # CoxPhRWRNet, loss function, C-index
│   ├── rwr.py                # Random Walk with Restart
│   ├── network.py            # PPI transition matrix
│   ├── features.py           # log-rank gene selection
│   ├── preprocessing.py      # normalisation and stratified splitting
│   ├── train.py              # training loop and 20-repeat evaluation
│   ├── mutations.py          # MAF parsing
│   ├── clinical.py           # clinical feature extraction
│   ├── data_loader.py        # file I/O
│   └── pipeline.py           # end-to-end grid-search pipeline
├── scripts/
│   ├── run_ablation.py       # component ablation study
│   ├── plot_km.py            # KM curves from saved risk scores
│   ├── run_benchmarks.py     # baseline comparisons
│   └── plot_benchmark.py     # Figure 4 boxplot
└── data/
    ├── BRCA/                 # included via Git LFS
    └── README.md             # download instructions for other cancers
```

---

## Requirements

Python 3.9 or higher. A GPU is recommended but not required.

```bash
pip install -r requirements.txt
```

Core dependencies: `torch>=2.0`, `numpy>=1.23`, `pandas>=1.5`, `scipy>=1.9`, `lifelines>=0.27`.  
Baseline comparisons additionally require `scikit-survival` and `pycox`.

---

## Data

BRCA data (mutations, clinical files, PPI network) is included via Git LFS and will be pulled automatically on `git clone`.

For LUNG, GBM, and OV, run:

```bash
python download_data.py --cancer LUNG
python download_data.py --cancer GBM
python download_data.py --cancer OV
```

This downloads mutation and clinical files from cBioPortal. Tissue-specific PPI files are available on request (contact: stutik@iiitd.ac.in). See `data/README.md` for the full file list and formats.

---

## Getting started

### Interactive notebook

Open `pipeline_rcoxnet.ipynb` in Jupyter and run the cells top to bottom. The notebook covers all steps from loading raw data to a trained model and C-index results. Set `CANCER = 'BRCA'` at the top (or switch to another cohort once you have the data).

### Command line

Reproduce the paper C-index numbers using the hyperparameters in `best_params.json`:

```bash
python run.py --reproduce --cancer BRCA
python run.py --reproduce --cancer ALL   # runs all four cancers
```

Run the full pipeline with joint hyperparameter search (α × K × architecture × LR × L2):

```bash
python run.py --cancer BRCA
```

### Kaplan-Meier curves (Figure 2)

First run Steps 1–8 in `pipeline_rcoxnet.ipynb` to generate per-patient risk scores, then:

```bash
python scripts/plot_km.py --cancer ALL
```

Or open `km_curves.ipynb` for an interactive version.

---

## Results

C-index (mean ± std, 20 splits):

| Cancer | C-index       |
|--------|---------------|
| BRCA   | 0.807 ± 0.045 |
| GBM    | 0.704 ± 0.042 |
| LUNG   | 0.750 ± 0.040 |
| OV     | 0.668 ± 0.037 |

---

## Citation

```bibtex
@article{kumari2025rcoxnet,
  title   = {RCoxNet: A deep learning framework integrating random walk with restart,
             mutation, and clinical data for cancer survival prediction},
  author  = {Kumari, Stuti and Gujral, Sakshi and Halder, Abhishek and
             Panda, Smruti and Mathew, Bernadette and Gupta, Prashant and
             Herwig, Ralf and Ahuja, Gaurav and Sengupta, Debarka},
  journal = {Journal of Computational Biology},
  year    = {2025},
  doi     = {10.1101/2024.09.17.613428}
}
```
