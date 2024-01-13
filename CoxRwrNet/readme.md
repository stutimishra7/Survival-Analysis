**RWR_data_bulder.py** 

- This Python script preprocesses TCGA breast cancer mutation and clinical data, identifies mutated genes, performs Random walk with restart analysis using R, integrates additional clinical features, and splits the data for the  model

**dataloader.py**
- This Python script defines functions for sorting genomic and clinical data based on survival time,Converts the sorted data to PyTorch tensors with the specified data type

**model.py**
- This PyTorch script defines a Cox proportional hazards regression model (`CoxPhRWRNet`) with a random walk restart (RWR) layer, utilizing genomic inputs, MSI scores, TMB scores, and patient age. The model comprises two hidden layers with tanh activation and a Cox layer incorporating additional features for predicting hazard ratios.

**train.py**
- The script defines a training function (`train_cox_rwrnnet`) for a Cox proportional hazards regression model with a random walk restart network, using negative partial log-likelihood as the loss function and Adam optimizer. The function returns training and evaluation losses, along with concordance indices during training.
