#!/usr/bin/env python

__author__ = "Raquel Parrondo-Pizarro"
__date__ = "20250220"
__copyright__ = "Copyright 2024, Chemotargets"
__license__ = ""
__credits__ = ["Data Science & Translational Research Group"]
__maintainer__ = "Raquel Parrondo-Pizarro"
__version__ = "20250220"
__deprecated__ = False

### Imports

import pandas as pd
import os

from DR_Utils import classify_skewness

### 
def __getSummaryDataFrame(self, n_mols_total, n_mols_sources_df, statistics_entire_dataset, statistics_sources_df, 
                          prop_ref_mols_total, prop_ref_mols_sources_df, skewness_total, skewness_sources_df, 
                          kurtosis_total, kurtosis_sources_df, outliers_total, outliers_sources_df, oor_total, oor_sources_df, 
                          endpoint_distribution_results_df, feature_similarity_results):
    
    """
    Generates the summary DataFrame by concatenating and merging the results from all
    conducted analysis.
    """

    # Define a dictionary containing the "entire-dataset" and "each-source" information to concatenate
    raw_dfs = {'n_mols': [n_mols_total, n_mols_sources_df], 'statistics': [statistics_entire_dataset, statistics_sources_df],
               'prop_ref_mols': [prop_ref_mols_total, prop_ref_mols_sources_df],
               'skewness': [skewness_total, skewness_sources_df], 'kurtosis': [kurtosis_total, kurtosis_sources_df],
               'outliers': [outliers_total, outliers_sources_df], 'oor': [oor_total, oor_sources_df]}
    concat_dfs = []

    for item, dataframes in raw_dfs.items():
        # Concatenate the 'prop_ref_mols' dataframes (if a reference set is defined)
        if item == 'prop_ref_mols' and self.reference_set is None:
            pass

        # Concatenate the 'skewness', 'kurtosis', 'outliers' and 'oor' dataframes (only for regression endpoints) 
        if item in ['skewness', 'kurtosis', 'outliers', 'oor'] and self.task == CONF_TASK_CLASSIFICATION:
            pass

        else:
            total_df = pd.DataFrame(data=[['entire_dataset']+dataframes[0]], columns=dataframes[1].columns)
            concat_df = pd.concat([total_df, dataframes[1]], ignore_index=True)
            concat_dfs.append(concat_df)
                                                                
    # Create a DataFrame summarizing feature similarity metrics for each dataset
    chemical_similarity_df = pd.DataFrame(data=[(item[0], sum(item[1]) / len(item[1]), sum(item[2]) / len(item[2])) for item in feature_similarity_results],
                                          columns=['dataset','average_within_similarity','average_between_similarity']) # TODO: Should we change the 'similarity' concept to 'distance'?
    concat_dfs.extend([endpoint_distribution_results_df, chemical_similarity_df])

    # Build the final dataframe
    for i in range(1, len(concat_dfs)):
        if i == 1:
            summary_df = pd.merge(concat_dfs[i-1], concat_dfs[i], on='dataset')
        else:
            summary_df = pd.merge(summary_df, concat_dfs[i], on='dataset', how='left')

    # Remove specific columns
    if self.task == CONF_TASK_REGRESSION and self.upper_bound is None:
        summary_df = summary_df.drop(columns=['n_upper_oor'])
    if self.task == CONF_TASK_REGRESSION and self.lower_bound is None:
        summary_df = summary_df.drop(columns=['n_lower_oor'])
    if self.task == CONF_TASK_CLASSIFICATION:
        summary_df = summary_df.drop(columns=['meandiff'])

    # Apply last changes
    summary_df[CONFIGS.NAMES.ENDPOINT_ID] = self.endpoint_name # add endpoint column
    summary_df = summary_df[[CONFIGS.NAMES.ENDPOINT_ID] + [col for col in summary_df if col != CONFIGS.NAMES.ENDPOINT_ID]] # reorder columns
    summary_df['dataset'] = pd.Categorical(summary_df['dataset'], categories=['entire_dataset'] + self.sources_list, ordered=True) # reorder rows
    summary_df = summary_df.sort_values('dataset').reset_index()

    return summary_df

###
def __writeToTSV(self, dataframe):

    """
    Formats numeric data and saves the summary DataFrame to TSV.
    """

    # Format values
    def format_values(val):
        try:
            return f"{val:.6e}" if 0 < abs(val) < 1e-6 else f"{round(val, 6)}" # TODO: Verify the correct value formatting 
        except (ValueError, TypeError):
            return val
    dataframe = dataframe.applymap(format_values)

    # Write the Pandas DataFrame as TSV 
    dataframe.to_csv(os.path.join(self.directory, self.endpoint_name, 'dataset_reporting.tsv'), sep='\t', index=False)

###
def __generateOutputSummary(self, results_df):

    """
    Generates a summary report for the resulting dataframe if all data sources were merged.
    
    It includes the total number of molecules, as well as several statistics specific to the 
    modeling task.  
    """
    
    # Select entire dataset data
    data = results_df.loc[results_df['dataset'] == 'entire_dataset']

    # Retrieve general information
    n_mols = f'\tNumber of molecules: {data["num_mols"].iloc[0]}\n'

    # Retrieve endpoint statistics and Create the final summary report
    description = f'\n### SUMMARY REPORT ###\nResulting {self.endpoint_name} dataset if all individual data sources were merged:\n'
    if self.task == CONF_TASK_REGRESSION:
        # Classify the skewness of the entire dataset
        data['skewness_category'] = data['skewness_value'].apply(classify_skewness)

        mean_std = f'\t{self.endpoint_name} Distribution: {data["mean"].iloc[0]:.3f} ± {data["standard_deviation"].iloc[0]:.3f}\n'
        range_iqr = f'\t{self.endpoint_name} Range: [{data["minimum"].iloc[0]:.3f} ; {data["maximum"].iloc[0]:.3f}] (IQR: [{data["1st_quartile"].iloc[0]:.3f} ; {data["3rd_quartile"].iloc[0]:.3f}])\n'
        n_outliers = f'\tNumber of outliers: {data["n_outliers_Zscore"].iloc[0]} ({data["prop_outliers_Zscore"].iloc[0]:.2f}%)\n'
        skewness = f'\tSkewness: {data["skewness_value"].iloc[0]:.2f} ({data["skewness_category"].iloc[0]})\n'
        
        n_upper_oor = n_lower_oor = ""
        if self.upper_bound is not None:
            n_upper_oor = f'\tNumber of upper OOR data points: {data["n_upper_oor"].iloc[0]}\n'
        if self.lower_bound is not None:
            n_lower_oor = f'\tNumber of lower OOR data points: {data["n_lower_oor"].iloc[0]}\n'

        output_summary = description + n_mols + mean_std + range_iqr + n_outliers + skewness + n_upper_oor + n_lower_oor

    elif self.task == CONF_TASK_CLASSIFICATION:
        class_counts = f'\tClass counts: {data["class_counts"].iloc[0]}\n'
        class_ratio = f'\tClass ratio: {data["class_ratio"].iloc[0]}\n'

        output_summary = description + n_mols + class_counts + class_ratio

    writeToLog('info', oriMessage=f'{output_summary}', toFile=True)