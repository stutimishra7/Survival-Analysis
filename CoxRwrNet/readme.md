**RWR_data_bulder.py**


2. **load_data():**
   - Reads clinical patient, clinical sample, and mutation data from three different files (`data_clinical_patient.txt`, `data_clinical_sample.txt`, and `data_mutations.txt`).
   - Returns the dataframes.

3. **clean_data(df_data_mut):**
   - Cleans mutation data by extracting the prefix from the 'Tumor_Sample_Barcode' column and saving it to a new CSV file named "clean_tcga.csv".
   - Modifies the 'Tumor_Sample_Barcode' column in the original dataframe.
 
4. **merge_data(df_clinical_patient, df_clinical_sample, df_data_mut_final):**
   - Merges clinical patient, clinical sample, and cleaned mutation data based on patients ID.
   - Returns the merged and filtered dataframe.

5. **get_all_mutated_genes(df_data_mut):**
   - Groups mutation data by 'Tumor_Sample_Barcode' and returns a dataframe with unique mutated genes for each sample.

6. **get_deleted_genes(df_data_mut):**
   - Reads data from "output_9606.protein.links.full.v11.5.txt" and identifies deleted genes.
   - Returns a list of deleted genes and the original list of genes.

7. **remove_items(test_list, item):**
   - Removes specified items from a list.

8. **calculate_seed(df_data_mut_patient_all_mutated, df_new, del_genes_all):**
   - Filters out deleted genes from the mutation data and saves the result to "patient_id_updated.csv".
   - Groups the updated mutation data by 'Tumor_Sample_Barcode' and returns a dataframe.
   - 
  9. **calculate_sg_score(list_of_genes, list_of_mutated_genes):**
   - Calculates a score for each gene based on its presence in the mutation data.
   - Saves the results to "score_sg_all_tcga_brca_all.csv".
   - Returns lists of mutated genes and their scores.
   - scores has not been use in code 

10. **run_r_script(input_csv_sg, input_csv_seed, output_csv):**
    - Executes an R script using rpy2, performing calculations based on input CSV files and saving the results to an output CSV file.
    - Returns generated CSV file and dataframe.

11. **Extract(lst):**
    - Extracts the first element from each list.

12. **process_and_clean_data(df_data_mut_patient_all_mutated, df_data_mut_patient_all_data_new, csv_output_path, data):**
    - Processes and cleans additional columns from the mutation data, adds new columns, and saves the cleaned data to a CSV file.
    - Returns the cleaned dataframe.

13. **split_data(data, train_ratio, test_ratio, validation_ratio, random_state):**
    - Splits the data into training, testing, and validation sets based on specified ratios.
    - Prints the size of each set and returns the sets.

14. **save_data_to_csv(data, filename):**
    - Saves the provided dataframe to a CSV file.

15. **call_funs():**
    - Calls all the defined functions in sequence to load, clean, merge, calculate, and process data.
    - Saves the training, testing, and validation sets to CSV files.
