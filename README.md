# RCoxNet: Network-diffusion powered Deep Cox Proportional Hazard model for improved survival prediction based on cancer mutations

To address the data sparsity issues of mutation data, we introduce RCoxNet (RCoxNet stands for Random walk with restart Cox-proportional Hazard Network model), a novel deep learning network-based model integrating the RWR algorithm utilizing TCGA mutation data, clinical and patient data. In this model, we have showcased the efficacy of the RWR algorithm in prioritizing candidate genes associated with diseases by leveraging the comprehensive network topology within a Protein-Protein Interaction (PPI) network. The hidden layers consist of three tiers, with decreasing nodes per layer from the first to the third. They aim to conduct feature extraction and discern salient information from the data. Additionally, the hidden layers address the non-linear effects present in high-dimensional data and allow for the capture of more intricate data. The 'tanh' activation function is employed within these hidden layers. The clinical layer integrates clinical information such as age, TMB score, and MSI score, serving as covariates directly incorporated into the Cox layer. The Cox layer incorporates the Cox Proportional Hazards (PH) model, a widely employed technique for analyzing survival data, enabling the simultaneous inclusion and assessment of the impact of multiple covariates, encompassing variables of specific research interest, such as treatment groups, and potential confounders including demographic and clinical factors.

![image](https://github.com/stutimishra7/Survival-Analysis/assets/70698090/afe7e14d-99c1-487f-9c7e-f4bebc6f71c8)


# Usage

The **RCoxNet** package is for improving survival outcomes using TCGA mutation data of cancer patients. the user can perform:

1. Preprocessing of TCGA cancer patient, clinical and mutation data and its usage for RWR algorithm
2. Perform model training and hyperparameter tuning
3. Calculate Prognostic index for patients of each cancer types.
4. Visualize survival curves for prognostic index, RWR score patterns and clinical features.

# Installation

## System Requirements
The deep learning models were trained on a system with the following specifications:

RAM: 8+ GB

CPU: 2 cores

GPU: 15360MiB

## Platforms used for training models
Linux (Ubuntu 22.04.3 LT)

Microsoft Windows 10

## Prerequisites
- pandas=1.5.3
- numpy=1.25.0
- lifelines==0.28.0
- matplotlib==3.8.2
- rpy2==3.5.15
- pytorch=2.0.1
- scikit-learn=1.2.2
- scipy=1.10.1

## Download Github repository
git clone https://github.com/stutimishra7/Survival-Analysis

# Contributing
For further information contact debarka@iiitd.ac.in or stutik@iiitd.ac.in.
