import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import pandas as pd
from sklearn.model_selection import train_test_split
import itertools
import numpy as np

from sklearn.model_selection import train_test_split
from IPython.display import display

import pandas as pd
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt
from lifelines.statistics import logrank_test
from matplotlib.backends.backend_pdf import PdfPages






# from google.colab import drive
# drive.mount('/content/drive')



def load_data():
    df_clinical_patient = pd.read_csv("data_clinical_patient.txt", sep="\t", skiprows=4)
    df_clinical_sample = pd.read_csv("data_clinical_sample.txt", sep="\t", comment='#', skiprows=1)
    df_data_mut = pd.read_csv("data_mutations.txt", sep="\t")

    return df_clinical_patient, df_clinical_sample, df_data_mut

def clean_data(df_data_mut):
    s = []
    l1 = df_data_mut['Tumor_Sample_Barcode'].tolist()
    for i in l1:
        s.append(str(i).rsplit('-01', 1)[0])
    d = pd.DataFrame(s, columns=["A"])
    d.to_csv("clean_tcga.csv")

    df_data_mut['Tumor_Sample_Barcode'] = d['A']
    # df_data_mut_final = df_data_mut[~df_data_mut['Hugo_Symbol'].isin(get_deleted_genes(df_data_mut))]

    return df_data_mut

def merge_data(df_clinical_patient, df_clinical_sample, df_data_mut_final):
    merged_inner = pd.merge(left=df_clinical_patient, right=df_clinical_sample, left_on='PATIENT_ID', right_on='PATIENT_ID')
    merged_inner_final = pd.merge(left=merged_inner, right=df_data_mut_final, left_on='PATIENT_ID', right_on='Tumor_Sample_Barcode')
    df_data_mut_final = merged_inner_final[merged_inner_final['CANCER_TYPE'] == 'Breast Cancer']

    return df_data_mut_final


def get_all_mutated_genes(df_data_mut):
    df_data_mut_patient_all_mutated = df_data_mut.groupby('Tumor_Sample_Barcode').apply(lambda x: x['Hugo_Symbol'].unique())
    return df_data_mut_patient_all_mutated



def get_deleted_genes(df_data_mut):
    df_data_genes = pd.read_csv("output_9606.protein.links.full.v11.5.txt")
    list_of_genes = df_data_genes['Gene_symbol'].to_list()
    item = '-'
    list_of_genes = remove_items(list_of_genes, item)
    # print(len(list_of_genes))

    unique_mutated_genes_all = set(list(itertools.chain(*get_all_mutated_genes(df_data_mut))))
    deleted_genes_all = np.setdiff1d(list(unique_mutated_genes_all), list_of_genes)
    del_genes_all = deleted_genes_all.tolist()

    return del_genes_all,list_of_genes



def remove_items(test_list, item):
    res = [i for i in test_list if i != item]
    return res

def calculate_seed(df_data_mut_patient_all_mutated,df_new, del_genes_all):

  # list_of_mutated_genes = df_data_mut_patient_all_mutated.to_list()
  a= del_genes_all
  df_new=df_new[~df_new['Hugo_Symbol'].isin(a)]


  df_data_mut_patient_all_data = df_new.groupby('Tumor_Sample_Barcode').apply(lambda x: x['Hugo_Symbol'].unique())
  df_data_mut_patient_all_data .to_csv("patient_id_updated.csv")
  list_of_mutated_genes=df_data_mut_patient_all_data.to_list()
  tcga_seed = pd.DataFrame(list_of_mutated_genes)
  tcga_seed.to_csv("seed_tcga_brca_all.csv",index=False,header=False)
  df_data_mut_patient_all_mutated .to_csv("patient_id_updated.csv")
  df_data_mut_patient_all_data_new=df_new.groupby('Tumor_Sample_Barcode')

  display(df_new)



  return list_of_mutated_genes, df_data_mut_patient_all_data_new


def calculate_sg_score(list_of_genes, list_of_mutated_genes):
    res_train = []
    gene_mutated_symbol_train = []

    for x in range(len(list_of_genes)):
        count = 0
        aa1 = 0

        for y in range(len(list_of_mutated_genes)):
            if list_of_genes[x] in list_of_mutated_genes[y]:
                count += 1
                aa1 += len(list_of_mutated_genes[y])

        if count > 0:
            avg_aa1 = aa1 / count
            res_train.append(avg_aa1)
            gene_mutated_symbol_train.append(list_of_genes[x])

    score_genes_all = pd.DataFrame(list(zip(gene_mutated_symbol_train, res_train)), columns=['NodeName', 'Score'])
    score_genes_all.to_csv("score_sg_all_tcga_brca_all.csv", index=False)
    return gene_mutated_symbol_train, res_train




def run_r_script(input_csv_sg, input_csv_seed, output_csv):
    r_script = f"""
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("RandomWalkRestartMH")

library(RandomWalkRestartMH)
library(tidyverse)
library(igraph)

sg_data <- read.csv('{input_csv_sg}', header=TRUE, na.strings="")
string_data <- read.csv('Gene_interaction_list1.csv', header=TRUE, na.strings="")
# print(string_data)
colnames(string_data)[3] <- 'weight'
delta = 0.3
string_graph = graph_from_data_frame(string_data, directed=F, vertices=NULL)
string_MultiplexObject = create.multiplex(list(STRING=string_graph))
AdjMatrix_string <- compute.adjacency.matrix(string_MultiplexObject, delta=delta)
AdjMatrixNorm_string <- normalize.multiplex.adjacency(AdjMatrix_string)

seed_data <- read.csv('{input_csv_seed}', header=FALSE, na.strings="")

Col_gene_list <- read.csv("Cosmic.csv", header = TRUE, na.strings = "")
final_gene_list <- unique(Col_gene_list$Genes)

score_data <- data.frame(matrix(nrow = nrow(seed_data), ncol = length(final_gene_list)))
colnames(score_data) <- final_gene_list

# Preallocate memory
l = vector("list", length = nrow(seed_data))
probs <- vector(mode = 'numeric', length = nrow(seed_data))
p <- vector(mode = 'list', length = nrow(seed_data))
l2 <- vector(mode = 'list', length = nrow(seed_data))
l3 <- vector(mode = 'list', length = nrow(seed_data))
l4 <- vector(mode = 'list', length = nrow(seed_data))

# Random Walk Restart Multiplex
for (i in 1:nrow(seed_data)) {{
  g <- seed_data[i, ]
  l[[i]] <- g[which(g != "")]
  rwr_result <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_string, string_MultiplexObject, l[[i]])
  probs[i] <- rwr_result$RWRM_Results["Score"]
  p[i] <- rwr_result$RWRM_Results["NodeNames"]

  df <- as.data.frame(rbind(probs[[i]], p[[i]]))
  l2[[i]] <- df
  df1 <- as.data.frame(t(df))
  l3[[i]] <- df1
  colnames(l3[[i]]) <- c("Score", "NodeName")
}}

# Sort data frames
l4 <- lapply(l3, function(df) df[order(df$Score, decreasing = TRUE), ])

# Merge data and calculate length
merge_data <- lapply(l4, function(df) {{
  merged <- merge(df, sg_data, by.x = "NodeName", by.y = "NodeName")
  len <- nrow(merged)
  list(data = merged, length = len)
}})

# Preallocate memory for len_of_list
len_of_list <- numeric(length = length(merge_data))

# Fill in score_data
for (i in seq_along(merge_data)) {{
  a <- merge_data[[i]]$data
  len_of_list[i] <- merge_data[[i]]$length
  for (j in seq_along(a$NodeName)) {{
    n <- a[j, 1]
    s <- a[j, 2]
    if (n %in% colnames(score_data)) {{
      pos <- which(colnames(score_data) == n)
      score_data[i, pos] <- s
      print("OK")
    }}
  }}
}}


result_path <- '{output_csv}'
write.csv(score_data, result_path)
result_path
"""

    return str(robjects.r(r_script))

def Extract(lst):
    return [item[0] for item in lst]

def process_and_clean_data(df_data_mut_patient_all_mutated,df_data_mut_patient_all_data_new, csv_output_path,data):
    AGE = df_data_mut_patient_all_data_new['AGE'].apply(list)
    os_months = df_data_mut_patient_all_data_new['OS_MONTHS'].apply(list)
    os_status = df_data_mut_patient_all_data_new['OS_STATUS'].apply(list)


    new_age = Extract(AGE)
    new_os = Extract(os_months)         # Extract first element from each list
    new_os_status = Extract(os_status)


    tmb = df_data_mut_patient_all_data_new['TMB_NONSYNONYMOUS'].apply(list)
    tmb_score = Extract(tmb)                                                     # Extract and clean additional columns

    msi = df_data_mut_patient_all_data_new['MSI_SENSOR_SCORE'].apply(list)
    msi_score = Extract(msi)

    all_index = df_data_mut_patient_all_mutated.index.to_list()

    # Load the existing dataset
    # data = pd.read_csv('/content/tcga_brca_cosmic_all.csv')


    null_columns = data.columns[data.isnull().all()]
    data = data.drop(axis=1, labels=null_columns)

    # Iterate over each column and fill null values
    for column in data.columns:
        if data[column].isnull().any():
            data[column].fillna(data[column].mean(), inplace=True)

    # Add last 4 columns
    data['SAMPLE_ID'] = all_index
    data['OS_MONTHS'] = new_os
    data['OS_STATUS'] = new_os_status
    data['AGE'] = new_age
    data['tmb_score'] = tmb_score
    data['msi_score'] = msi_score

    # Write cleaned data to a new CSV file
    data['OS_STATUS'] = data['OS_STATUS'].replace('0:LIVING', '0')
    data['OS_STATUS'] = data['OS_STATUS'].replace('1:DECEASED', '1')
    data.to_csv(csv_output_path, index=False)
    return data




def split_data(data, train_ratio=0.7, test_ratio=0.2, validation_ratio=0.1, random_state=42):

    train_data, remaining_data = train_test_split(data, train_size=train_ratio, random_state=random_state)


    test_data, validation_data = train_test_split(remaining_data, test_size=validation_ratio / (test_ratio + validation_ratio), random_state=random_state)

    print(f"Train set size: {len(train_data)}")
    print(f"Test set size: {len(test_data)}")
    print(f"Validation set size: {len(validation_data)}")

    return train_data, test_data, validation_data

def save_data_to_csv(data, filename):
    data.to_csv(filename, index=False)

# original_data = pd.read_csv("entire_data_tcga_Brca_all_cosmic_tmb_an_all.csv")
# train_data, test_data, validation_data = split_data(original_data)



def cll_funs():

    df_clinical_patient, df_clinical_sample, df_data_mut = load_data()

    df_data_mut = clean_data(df_data_mut)
    # display(df_data_mut_final )
    df_data_mut_final = merge_data(df_clinical_patient, df_clinical_sample, df_data_mut)
    # display(df_data_mut_final)
    df_data_mut_patient_all_mutated=get_all_mutated_genes(df_data_mut_final)
    # display(df_data_mut_patient_all_mutated)
    del_genes_all,list_of_genes=get_deleted_genes(df_data_mut)
    # display(len(del_genes_all))

    list_of_mutated_genes, df_data_mut_patient_all_data_new = calculate_seed(df_data_mut_patient_all_mutated,df_data_mut_final, del_genes_all)
    print(df_data_mut_patient_all_data_new)
    gene_mutated_symbol_train, res_train = calculate_sg_score(list_of_genes, list_of_mutated_genes)
    csv_output_path = "entire_data_tcga_Brca_all_cosmic_tmb.csv"

    input_csv_sg = "score_sg_all_tcga_brca_all.csv"
    input_csv_seed = "seed_tcga_brca_all.csv"
    output_csv = "output_tcga_brca_cosmic_all.csv"
    rwr_data = run_r_script(input_csv_sg, input_csv_seed, output_csv)
    data = process_and_clean_data(df_data_mut_patient_all_mutated,df_data_mut_patient_all_data_new, csv_output_path,rwr_data)
    train_data, test_data, validation_data= split_data(data, train_ratio=0.7, test_ratio=0.2, validation_ratio=0.1, random_state=42)
    # Save to CSV files
    save_data_to_csv(train_data, 'train_tcga_Brca_all_cosmic_tmb.csv')
    save_data_to_csv(test_data, 'test_final_tcga_Brca_all_cosmic_tmb.csv')
    save_data_to_csv(validation_data, 'validation_data_tcga_Brca_all_cosmic_tmb.csv')
    # run_survival_analysis()



