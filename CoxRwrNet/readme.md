**RWR_data_bulder.py** 

- This Python script preprocesses TCGA breast cancer mutation and clinical data, identifies mutated genes, performs Random walk with restart analysis using R, integrates additional clinical features, and splits the data for the  model

**dataloader.py**
- This Python script defines functions for sorting genomic and clinical data based on survival time,Converts the sorted data to PyTorch tensors with the specified data type

**model.py**
- This PyTorch script defines a Cox proportional hazards regression model (`CoxPhRWRNet`) with a random walk restart (RWR) layer, utilizing genomic inputs, MSI scores, TMB scores, and patient age. The model comprises two hidden layers with tanh activation and a Cox layer incorporating additional features for predicting hazard ratios.

**train.py**
- The script defines a training function (`train_cox_rwrnnet`) for a Cox proportional hazards regression model with a random walk restart network, using negative partial log-likelihood as the loss function and Adam optimizer. The function returns training and evaluation losses, along with concordance indices during training.
  
**Run_train.py**
- This script loads data, conducts hyperparameter grid search, trains a Cox proportional hazards model using optimal parameters, and prints the best hyperparameters along with the concordance index on the test set.

**CostFunc_CIndex.py**
- The script provides PyTorch functions for creating an indicator matrix of risk sets, computing the negative partial log-likelihood, and evaluating the concordance index.
  
**Run_for_entiredata.py**
- Conducts interpretation of model on the entire dataset, saving key model details, weights, and node values. Linear predictions and survival data are combined and stored for analysis.

**Run_Survival.py**
- This script generates Kaplan-Meier survival plots for high and low-risk groups based on a median cutoff of a prognostic index (PI) calculated from a Cox proportional hazards model. It includes statistical analysis using the log-rank test and saves the plots in a PDF file.
