**RWR_data_bulder.py** 

- This Python script preprocesses TCGA breast cancer mutation and clinical data, identifies mutated genes, performs Random walk with restart analysis using R, integrates additional clinical features, and splits the data for the  model

**dataloader.py**
- This Python script defines functions for sorting genomic and clinical data based on survival time,Converts the sorted data to PyTorch tensors with the specified data type
