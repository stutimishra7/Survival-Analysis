Explanation of RCoxNet parameters:

<b>RWR_data_bulder.py</b>

load_data():

Reads clinical patient, clinical sample, and mutation data from three different files (data_clinical_patient.txt, data_clinical_sample.txt, and data_mutations.txt).
Returns the dataframes.
clean_data(df_data_mut):

Cleans mutation data by extracting the prefix from the 'Tumor_Sample_Barcode' column and saving it to a new CSV file named "clean_tcga.csv".
Modifies the 'Tumor_Sample_Barcode' column in the original dataframe.
merge_data(df_clinical_patient, df_clinical_sample, df_data_mut_final):

Merges clinical patient, clinical sample, and cleaned mutation data based on patients ID.
Returns the merged and filtered dataframe.
get_all_mutated_genes(df_data_mut):

Groups mutation data by 'Tumor_Sample_Barcode' and returns a dataframe with unique mutated genes for each sample.
get_deleted_genes(df_data_mut):

Reads data from "output_9606.protein.links.full.v11.5.txt" and identifies deleted genes.
Returns a list of deleted genes and the original list of genes.
remove_items(test_list, item):

Removes specified items from a list.
calculate_seed(df_data_mut_patient_all_mutated, df_new, del_genes_all):

Filters out deleted genes from the mutation data and saves the result to "patient_id_updated.csv".
Groups the updated mutation data by 'Tumor_Sample_Barcode' and returns a dataframe.
calculate_sg_score(list_of_genes, list_of_mutated_genes):

Calculates a score for each gene based on its presence in the mutation data.
Saves the results to "score_sg_all_tcga_brca_all.csv".
Returns lists of mutated genes and their scores.
scores has not been use in code
run_r_script(input_csv_sg, input_csv_seed, output_csv):

Executes an R script using rpy2, performing calculations based on input CSV files and saving the results to an output CSV file.
Returns generated CSV file and dataframe.
Extract(lst):

Extracts the first element from each list.
process_and_clean_data(df_data_mut_patient_all_mutated, df_data_mut_patient_all_data_new, csv_output_path, data):

Processes and cleans additional columns from the mutation data, adds new columns, and saves the cleaned data to a CSV file.
Returns the cleaned dataframe.
split_data(data, train_ratio, test_ratio, validation_ratio, random_state):

Splits the data into training, testing, and validation sets based on specified ratios.
Prints the size of each set and returns the sets.
save_data_to_csv(data, filename):

Saves the provided dataframe to a CSV file.
call_funs():

Calls all the defined functions in sequence to load, clean, merge, calculate, and process data.
Saves the training, testing, and validation sets to CSV files.
<p><b>dataloader.py</b></p>

load_sorted_data(file_path, tensor_dtype):

Calls the sort_genomic_clinical_data function to obtain sorted data.

Converts the sorted data to PyTorch tensors with the specified tensor_dtype.

Returns the following PyTorch tensors: X: Genomic inputs. YTIME: Survival time. YEVENT: Censoring status. AGE: Age data. MSI: MSI data. TMB: TMB data.

<p><b>model.py</b></p>

RCoxNet Class:
Constructor (__init__):

Initializes the neural network architecture.
Parameters:
input_nodes: Number of input nodes (features) for genomic data.
hidden_nodes1: Number of nodes in the first hidden layer.
hidden_nodes2: Number of nodes in the second hidden layer.
output_nodes: Number of nodes in the output layer.
Attributes:

tanh: The hyperbolic tangent activation function (nn.Tanh()).

rwr_layer: The linear layer for the Random Walk Restart (RWR) method, transforming genomic input.

hidden_layer1: The first hidden linear layer.

hidden_layer2: The second hidden linear layer.

cox_layer: The linear layer for Cox Proportional Hazard model, combining the RWR output and additional features (age, MSI, TMB).

Methods (forward):

Defines the forward pass of the neural network.
Takes genomic data (x_genomic), MSI data (x_msi), TMB data (x_tmb), and age data (x_age) as inputs.
Applies the tanh activation function to the RWR layer and the two hidden layers.
Concatenates the hidden layer 2 output with additional features.
Passes the combined features through the Cox Proportional Hazard layer.
Returns the output of the Cox Proportional Hazard layer.
This architecture is designed for survival analysis, combining genomic data with additional features and using the Random Walk Restart method in the initial layers. The final layer produces the output for the Cox Proportional Hazard model.

<p><b>train.py</b></p>

Trains a Cox proportional hazards model with a neural network structure.

Parameters:

train_x, eval_x: Training and evaluation input features.

train_age, eval_age: Age information for training and evaluation.

train_ytime, eval_ytime: Time-to-event for training and evaluation.

train_yevent, eval_yevent: Event indicator for training and evaluation.

train_msi, eval_msi: Microsatellite instability for training and evaluation.

train_tmb, eval_tmb: Tumor mutational burden for training and evaluation.

In_Nodes, hidden_nodes1, hidden_nodes2, Out_Nodes: Neural network architecture parameters.

Learning_Rate: Learning rate for optimization.

L2: L2 regularization parameter.

Num_Epochs: Number of training epochs.

Returns:

Tuple containing training loss, evaluation loss, training concordance index, and evaluation concordance index.

<p><b>Run_train.py</b></p>
Train RCoxNet with optimal hyperparameters using train data, and then evaluate the trained model with test data
Note that test data are only used to evaluate the trained RCoxNet.

<b><p>Run_Survival.py</p></b>
plot survival plot from output of the RCoxNet model


