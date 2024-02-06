
import numpy as np
import pandas as pd
import torch

def sort_genomic_clinical_data(file_path):
    '''Arrange genomic and clinical data based on survival time (OS_MONTHS) in descending order.
Parameters:
file_path: File path to the input dataset in CSV format.
Returns:
genomic_inputs: Genomic inputs sorted accordingly.
survival_time: Survival time (OS_MONTHS) sorted in relation to 'genomic_inputs'.
censoring_status: Censoring status (OS_EVENT) sorted corresponding to 'genomic_inputs' (1 for deceased, 0 for censored).
patient_age: Age data sorted in relation to 'genomic_inputs'.
msi_score: MSI (Microsatellite Instability) data sorted in relation to 'genomic_inputs'.
tmb_score: TMB (Tumor Mutational Burden) data sorted in relation to 'genomic_inputs'.
'''

    data = pd.read_csv(file_path)

    data.sort_values("OS_MONTHS", ascending=False, inplace=True)

    genomic_inputs = data.drop(["SAMPLE_ID", "OS_MONTHS", "OS_EVENT", "AGE", "tmb_score", "msi_score"], axis=1).values
    # genomic_inputs = data.drop(["OS_MONTHS", "OS_EVENT", "AGE", "tmb_score", "msi_score"], axis=1).values

    survival_time = data.loc[:, ["OS_MONTHS"]].values
    censoring_status = data.loc[:, ["OS_EVENT"]].values
    patient_age = data.loc[:, ["AGE"]].values
    msi_score = data.loc[:, ["msi_score"]].values
    tmb_score = data.loc[:, ["tmb_score"]].values

    return genomic_inputs, survival_time, censoring_status, patient_age, msi_score, tmb_score


def load_sorted_data(file_path, tensor_dtype):
    '''Load the sorted data and convert it to a Pytorch tensor.
    Input:
        file_path: File path to the input dataset (expected to be a CSV file).
        tensor_dtype: Define the data type of the tensor (e.g., tensor_dtype=torch.FloatTensor).
    Output:
        X: Pytorch tensor of 'genomic_inputs' from sort_genomic_clinical_data().
        YTIME: Pytorch tensor of 'survival_time' from sort_genomic_clinical_data().
        YEVENT: Pytorch tensor of 'censoring_status' from sort_genomic_clinical_data().
        AGE: Pytorch tensor of 'patient_age' from sort_genomic_clinical_data().
        MSI: Pytorch tensor of 'msi_score' from sort_genomic_clinical_data().
        TMB: Pytorch tensor of 'tmb_score' from sort_genomic_clinical_data().
    '''
    genomic_inputs, survival_time, censoring_status, patient_age, msi_score, tmb_score = sort_genomic_clinical_data(file_path)
    ###

    X = torch.from_numpy(genomic_inputs).type(tensor_dtype)
    YTIME = torch.from_numpy(survival_time).type(tensor_dtype)
    YEVENT = torch.from_numpy(censoring_status).type(tensor_dtype)
    AGE = torch.from_numpy(patient_age).type(tensor_dtype)
    MSI = torch.from_numpy(msi_score).type(tensor_dtype)
    TMB = torch.from_numpy(tmb_score).type(tensor_dtype)
    ###if GPU is being used
    if torch.cuda.is_available():
        X = X.cuda()
        YTIME = YTIME.cuda()
        YEVENT = YEVENT.cuda()
        AGE = AGE.cuda()
        MSI = MSI.cuda()
        TMB = TMB.cuda()
    ###
    return X, YTIME, YEVENT, AGE, MSI, TMB
