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
import numpy as np

from AI_Utils import classify_skewness, formatWarningTitle

### Configs
import os
import json

config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'configs.json')

with open(config_path, "r") as file:
    config = json.load(file)

class InsightReport():
    """
    Class to generate the insight report file containig several warnings and recommendations
    to guide the data integration process.
    """

    def __init__(self, mainSelf):
        self.__mainSelf = mainSelf

    ###
    def __checkEndpointDistribution(self, endpoint_distribution_results, meandiff_thresh, alpha=0.01):

        """
        Identifies individual datasets whose endpoint distribution significantly differ from the
        others (p-value < alpha: 0.01) (1vsOthers setting). 
        
        In addition, filters out those with a mean difference greater than one standard deviation
        from the mean.

        This type of analysis can only be conducted for regression endpoints.
        """

        # Select the datasets with a p-value bewlo alpha (default: 0.01)
        diff_distribution_df = endpoint_distribution_results.loc[endpoint_distribution_results['p_value'] < alpha]
        # Filter the datasets based on the defined threshold
        diff_distribution_df = diff_distribution_df.loc[diff_distribution_df['meandiff'].abs() >= meandiff_thresh]
        # Sort the filtered datasets by their mean difference
        diff_distribution_df = diff_distribution_df.loc[diff_distribution_df['meandiff'].abs().sort_values(ascending=False).index].reset_index()

        return diff_distribution_df

    ###
    def __checkValueRange(self, data_df):

        """
        Identifies individual datasets where the endpoint value range, expressed in terms of the 
        interquartile range (IQR), is greater than one order of magnitude compared to the IQR calculated
        taking the rest of the data (1vsOthers setting). The use of the IQR instead of the full
        range of values enables the exclusion of outliers. 

        This type of analysis can only be conducted for regression endpoints.
        """

        # Select the datasets with an inconsistent value range
        diff_iqr_results = []
        for ref in self.__mainSelf.sources_list:
            source_data = data_df.loc[data_df[config["NAMES"]["REF"]] == ref]
            other_sources_data = data_df.loc[data_df[config["NAMES"]["REF"]] != ref]
            # Calculate both IQR's
            iqr_source = source_data[config["NAMES"]["VALUE"]].describe()['75%'] - source_data[config["NAMES"]["VALUE"]].describe()['25%']
            iqr_other_sources = other_sources_data[config["NAMES"]["VALUE"]].describe()['75%'] - other_sources_data[config["NAMES"]["VALUE"]].describe()['25%']
            # Check if the source IQR is greater than one order of magnitude compared to other-sources IQR
            if np.floor(np.log10(iqr_source)) != np.floor(np.log10(iqr_other_sources)):
                diff_iqr_results.append([ref, iqr_source, iqr_other_sources])

        # Generate the resulting DataFrame
        diff_iqr_df = pd.DataFrame(data=diff_iqr_results, columns=[config["NAMES"]["REF"], 'iqr_source', 'iqr_entire_dataset'])

        return diff_iqr_df
        
    ###
    def __checkSkewness(self, skewness_results):

        """
        Identifies individual datasets whose endpoint distribution is skewed and classify them 
        into categories of:
        - Severe left skewness --> skewness < -2.0
        - Moderate left skewness --> skewness < -1.0
        - Moderate right skewness --> skewness > 1.0
        - Severe right skewness --> skewness > 2.0

        This type of analysis can only be conducted for regression endpoints.
        """

        # Filter the skewed datasets (skeweness < -1.0 or > 1.0)
        skweness_df = skewness_results.loc[(skewness_results['skewness_value'] < -1.0) | (skewness_results['skewness_value'] > 1.0)]
        # Sort the datasets by their skewness
        skweness_df = skweness_df.sort_values(by='skewness_value', ascending=True).reset_index()
        # Classify them into the corresponding category
        skweness_df['skewness_category'] = skweness_df['skewness_value'].apply(classify_skewness)
        
        return skweness_df   

    ###
    def __checkOutliers(self, outliers_results):

        """
        Identifies the individual datasets with a proportion of outliers (calculated using the 
        Z-score method) greater than 1%.

        This type of analysis can only be conducted for regression endpoints.
        """

        # Filter the datasets with an outlier proportion greater than 1%
        outliers_df = outliers_results.loc[outliers_results['prop_outliers_Zscore'] >= 1]
        # Sort the datasets by their outlier proportion
        outliers_df = outliers_df.sort_values(by='prop_outliers_Zscore', ascending=False).reset_index()
        
        return outliers_df 

    ###
    def __checkAbnormalDataPoints(self, outliers_results):

        """
        Identifies individual datasets containing abnormal data, defined as those data points with 
        a Z-score greater than 5 standard deviations.

        This type of analysis can only be conducted for regression endpoints.    
        """

        # Filter the datasets with abnormal values
        abnormals_df = outliers_results.loc[outliers_results['n_abnormals_Zscore'] > 0]
        # Sort the datasets by their number of abnormal data points
        abnormals_df = abnormals_df.sort_values(by='n_abnormals_Zscore', ascending=False).reset_index()
        
        return abnormals_df 
        
    ###
    def __checkOORDataPoints(self, oor_results):

        """
        Identifies the individual datasets with either upper or lower out-of-range (OOR) values.
        
        This type of analysis can only be conducted for regression endpoints and when either or
        both bounds are defined.
        """

        # Filter the datasets with OOR values
        oor_df = oor_results.loc[(oor_results['n_upper_oor'] > 0) | (oor_results['n_lower_oor'] > 0)]
        
        return oor_df 

    ###
    def __identifyDissimilarDatasets(self, feature_similarity_results):

        """
        Identifies dissimilar datasets, defined as those with an average inter-similarity above
        one standard deviation from the mean. 
        """

        # Create a DataFrame with the results of the feature similarity analysis
        feature_similarity_df = pd.DataFrame(data=[(item[0], sum(item[1]) / len(item[1]), sum(item[2]) / len(item[2])) for item in feature_similarity_results],
                                            columns=[config["NAMES"]["REF"],'average_within_similarity','average_between_similarity'])
        # Calculate the mean and standard deviation of the distribution of between-source similarities
        mean = np.mean([element for item in feature_similarity_results for element in item[2]])
        std = np.std([element for item in feature_similarity_results for element in item[2]])
        # Calculate the similarity threshold
        similarity_cutoff = mean + 1*std
        # Select the dissimilar datasets
        dissimilar_sources_df = feature_similarity_df.loc[feature_similarity_df['average_between_similarity'] >= similarity_cutoff]
        # Sort the datasets by their average inter-similarity
        dissimilar_sources_df = dissimilar_sources_df.sort_values(by='average_between_similarity', ascending=False).reset_index()
        
        return dissimilar_sources_df

    ###
    def __identifyConflictingDatasets(self, discrepancies_df):

        """
        Identifies pairs of confilcting datasets, defined as those sharing more than 50% of their
        molecules and with endpoint value differences exceeding 10%.
        """

        # Endpoint type specifications
        if self.__mainSelf.task == config["TASK_CLASSIFICATION"]:
            diff_col = 'different_values_count_proportion'
        elif self.__mainSelf.task == config["TASK_REGRESSION"]:
            diff_col = 'mean_abs_diff'
        # Select the conflicting dataset pairs
        conflicting_sources_df = discrepancies_df.loc[(discrepancies_df['percent_common_mols'] >= 50) & (discrepancies_df[diff_col] >= 0.1)].reset_index()

        return conflicting_sources_df

    ###
    def __identifyDivergentDatasets(self, data_df):

        """
        Identifies divergent datasets, defined as those sharing less than 10% of their molecules
        with the other sources (1vsOthers setting).
        """

        # Calculate the molecule overlap in 1vsOthers setting
        percent_common_mols_results = []
        for ref in self.__mainSelf.sources_list:
            source_mols = data_df[config["NAMES"]["INCHIKEY"]].loc[data_df[config["NAMES"]["REF"]] == ref].to_list()
            other_sources_mols = data_df[config["NAMES"]["INCHIKEY"]].loc[data_df[config["NAMES"]["REF"]] != ref].to_list()
            # Find common molecules 
            common_mols = set(source_mols).intersection(set(other_sources_mols))
            # Compute the percentage of common molecules
            percent_common_mols = (len(common_mols)/len(source_mols))*100
            percent_common_mols_results.append([ref, percent_common_mols])

        # Generate the resulting DataFrame
        percent_common_mols_df = pd.DataFrame(data=percent_common_mols_results, columns=[config["NAMES"]["REF"], 'percent_common_mols'])
        # Select the divergent datasets
        divergent_sources_df = percent_common_mols_df.loc[percent_common_mols_df['percent_common_mols'] <= 10].reset_index()
        # Sort the datasets by their percentage of common molecules
        divergent_sources_df = divergent_sources_df.sort_values(by='percent_common_mols', ascending=False).reset_index()
        
        return divergent_sources_df

    ###
    def __identifyRedundantDatasets(self, discrepancies_df):

        """
        Identifies pairs of redundant datasets, defined as those sharing more than 80% of their
        molecules.
        """

        # Select the redundant dataset pairs
        redundant_sources_df = discrepancies_df.loc[discrepancies_df['percent_common_mols'] >= 80]
        # Sort the datasets by their percentage of common molecules
        redundant_sources_df = redundant_sources_df.sort_values(by='percent_common_mols', ascending=False).reset_index()

        return redundant_sources_df

    ###
    def __identifyUnderrepresentedDatasets(self, prop_ref_mols_df):

        """
        Identifies underrepresented datasets, defined as those with a proportion of reference molecules
        less than 15%.
        
        This type of analysis can only be conducted when a reference set of molecules is defined.
        """

        # Select the underrepresented datasets
        underrepresented_sources_df = prop_ref_mols_df.loc[prop_ref_mols_df['prop_ref_mols'] <= 15]
        # Sort the datasets by their proportion of reference molecules
        underrepresented_sources_df = underrepresented_sources_df.sort_values(by='prop_ref_mols', ascending=False).reset_index()
        
        return underrepresented_sources_df.rename(columns={'dataset': 'ref'})

    ###
    def __writeToTXT(self, diff_distribution_df, diff_range_df, skewed_distribution_df, outliers_identification_df, 
                    abnormal_identification_df, oor_identification_df, dissimilar_sources_df, conflicting_sources_df, 
                    divergent_sources_df, redundant_sources_df, underrepresented_sources_df):
        
        """
        Generates a TXT file with the insight report results.
        """

        # Endpoint type specifications    
        if self.__mainSelf.task == config["TASK_CLASSIFICATION"]:
            diff_col = 'different_values_count_proportion'
            diff_name = 'Proportion of different values'
        elif self.__mainSelf.task == config["TASK_REGRESSION"]:
            diff_col = 'mean_abs_diff'
            diff_name = 'Mean Absolute Difference'

        with open(self.__mainSelf.directory+'/'+self.__mainSelf.endpoint_name+'/insight_report.txt','w') as f:
            f.write('###########################################################################\n')
            f.write('############################# INSIGHT REPORT ##############################\n')
            f.write('###########################################################################\n')

            if self.__mainSelf.task == config["TASK_REGRESSION"]:
                # WARNING 1: Endpoint Distribution Warning
                f.write(formatWarningTitle(text='Endpoint Distribution'))
                if diff_distribution_df.shape[0] > 0:
                    f.write(f'Data sources with a significantly different {self.__mainSelf.endpoint_name} distribution:\n')
                    for i,row in diff_distribution_df.iterrows():
                        f.write(f'{i+1}) {row["dataset"]} --> Mean difference: {round(row["meandiff"], 3)}\n')
                        if i+1 == diff_distribution_df.shape[0]:
                            f.write(f'\t***Figure reference: {self.__mainSelf.endpoint_name} Distribution Across Data Sources (Violin Plots) [Fig 3.2].\n')  
                else:                    
                    f.write(f'No data source shows a significantly different {self.__mainSelf.endpoint_name} distribution compared to others.\n')

                # WARNING 1.1: Value Range Warning
                f.write(formatWarningTitle(text='Value Range'))
                if diff_range_df.shape[0] > 0:
                    f.write(f'Data sources with an inconsistent {self.__mainSelf.endpoint_name} value range (> 1 order of magnitude):\n')
                    for i,row in diff_range_df.iterrows():
                        f.write(f'{i+1}) {row[config["NAMES"]["REF"]]} --> Source\'s IQR: {round(row["iqr_source"], 2)} vs. Other-sources\' IQR: {round(row["iqr_entire_dataset"], 2)}\n')
                        if i+1 == diff_range_df.shape[0]:
                            f.write(f'\t***Figure reference: {self.__mainSelf.endpoint_name} Distribution Across Data Sources (Violin Plots) [Fig 3.2].\n')  
                else:
                    f.write(f'No data source shows an inconsistent {self.__mainSelf.endpoint_name} value range (> 1 order of magnitude).\n')

                # WARNING 1.2: Skewed Distribution Warning
                f.write(formatWarningTitle(text='Skewed Distribution'))
                if skewed_distribution_df.shape[0] > 0:                
                    f.write(f'Data sources with a skewed {self.__mainSelf.endpoint_name} distribution:\n')
                    for i,row in skewed_distribution_df.iterrows():
                        f.write(f'{i+1}) {row["dataset"]} --> Skewness: {round(row["skewness_value"], 2)} ({row["skewness_category"]})\n')
                        if i+1 == skewed_distribution_df.shape[0]:
                            f.write(f'\t***Figure reference: {self.__mainSelf.endpoint_name} Distribution Across Data Sources (Overlaid Histogram) [Fig 3.1].\n')  
                else:
                    f.write(f'No data source shows a skewed {self.__mainSelf.endpoint_name} distribution.\n')
                
                # WARNING 1.3: Outliers Warning
                f.write(formatWarningTitle(text='Outliers Identification'))
                if outliers_identification_df.shape[0] > 0:                
                    f.write(f'Data sources with a significant proportion of outliers (≥ 1%):\n')
                    for i,row in outliers_identification_df.iterrows():
                        f.write(f'{i+1}) {row["dataset"]} --> Outliers Proportion: {round(row["prop_outliers_Zscore"], 2)}% ({row["n_outliers_Zscore"]} outliers)\n')
                        if i+1 == outliers_identification_df.shape[0]:
                            f.write(f'\t***Figure reference: Visualization of {self.__mainSelf.endpoint_name} Outliers [Fig 1].\n')  
                else:
                    f.write(f'No data source shows a significant proportion of outliers (≥ 1%).\n')

                # WARNING 1.4: Abonormal Data Warning
                f.write(formatWarningTitle(text='Abnormal Data Identification'))
                if abnormal_identification_df.shape[0] > 0:                
                    f.write(f'Data sources with abnormal values:\n')
                    for i,row in abnormal_identification_df.iterrows():
                        f.write(f'{i+1}) {row["dataset"]} --> Abnormals Proportion: {round(row["prop_abnormals_Zscore"], 2)}% ({row["n_abnormals_Zscore"]} data points)\n')
                        if i+1 == abnormal_identification_df.shape[0]:
                            f.write(f'\t***Figure reference: Visualization of {self.__mainSelf.endpoint_name} Outliers [Fig 1].\n')  
                else:
                    f.write(f'No data source contains abnormal values.\n')

                # WARNING 1.5: Out-Of-Range Data Warning
                if (self.__mainSelf.upper_bound is not None) or (self.__mainSelf.lower_bound is not None):
                    f.write(formatWarningTitle(text='Out-Of-Range (OOR) Data Identification'))
                    if oor_identification_df.shape[0] > 0:                
                        f.write(f'Data sources with out-of-range data points:\n')
                        for i,row in oor_identification_df.iterrows():
                            if (self.__mainSelf.upper_bound is not None) and (self.__mainSelf.lower_bound is not None):
                                f.write(f'{i+1}) {row["dataset"]} --> Nº of upper OOR values: {row["n_upper_oor"]} / Nº of lower OOR values: {row["n_lower_oor"]}\n')
                            elif self.__mainSelf.upper_bound is not None:
                                f.write(f'{i+1}) {row["dataset"]} --> Nº of upper OOR values: {row["n_upper_oor"]}\n')
                            elif self.__mainSelf.lower_bound is not None:
                                f.write(f'{i+1}) {row["dataset"]} --> Nº of lower OOR values: {row["n_lower_oor"]}\n')
                    else:
                        f.write(f'No data source shows out-of-range values.\n')

            # WARNING 2: Dissimilar Datasets Warning
            f.write(formatWarningTitle(text='Dissimilar Datasets'))
            if dissimilar_sources_df.shape[0] > 0:    
                f.write(f'Data sources which are dissimilar in terms of the molecular feature profile:\n')
                for i,row in dissimilar_sources_df.iterrows():
                    f.write(f'{i+1}) {row[config["NAMES"]["REF"]]} --> Average between-source distance: {round(row["average_between_similarity"], 2)}\n')
                    if i+1 == dissimilar_sources_df.shape[0]:
                        f.write(f'\t***Figure reference: Feature Space Coverage (KDE plot) [Fig 6].\n') # Pairwise Between-Source Similarity plot could also be used as reference
            else:
                f.write(f'No data source is dissimilar in terms of the molecular feature profile.\n') 

            # WARNING 3: Conflicting Datasets Warning
            f.write(formatWarningTitle(text='Conflicting Datasets'))
            if conflicting_sources_df.shape[0] > 0:
                f.write(f'Data source pairs with conflicting {self.__mainSelf.endpoint_name} value annotations for shared molecules\n(>50% molecules in common + endpoint differences > 10%):\n')
                for i,row in conflicting_sources_df.iterrows():
                    f.write(f'{i+1}) {row["source_1"]} vs. {row["source_2"]} --> Percentage of common molecules: {round(row["percent_common_mols"], 2)}% - {diff_name} : {round(row[diff_col], 2)}\n')
                    if i+1 == conflicting_sources_df.shape[0]:
                        f.write(f'\t***Figure reference: Between-Source Discrepancies [Fig 9].\n') 
            else:
                f.write(f'No data source pair has conflicting {self.__mainSelf.endpoint_name} value annotations for shared molecules.\n')

            # WARNING 4: Divergent Datasets Warning
            f.write(formatWarningTitle(text='Divergent Datasets'))
            if divergent_sources_df.shape[0] > 0:
                f.write(f'Data sources which are divergent due to a low overlap with the rest of sources (<10% shared molecules):\n')
                for i,row in divergent_sources_df.iterrows():
                    f.write(f'{i+1}) {row[config["NAMES"]["REF"]]} --> Percentage of common molecules: {round(row["percent_common_mols"], 2)}%\n')
                    if i+1 == divergent_sources_df.shape[0]:
                        f.write(f'\t***Figure reference: Between-Source Discrepancies [Fig 9].\n')
            else:
                f.write(f'No data source has a low overlap of molecules (less than 10%) with the rest of sources.\n')

            # WARNING 5: Redundant Datasets Warning
            f.write(formatWarningTitle(text='Redundant Datasets'))
            if redundant_sources_df.shape[0] > 0:
                f.write(f'Data source pairs which are redundant due to significant proportion of shared molecules (>80% molecules in common):\n')
                for i,row in redundant_sources_df.iterrows():
                    f.write(f'{i+1}) {row["source_1"]} vs. {row["source_2"]} --> Percentage of common molecules (from first source): {round(row["percent_common_mols"], 2)}%\n')
                    if i+1 == redundant_sources_df.shape[0]:
                        f.write(f'\t***Figure reference: Between-Source Discrepancies [Fig 9].\n') 
            else:
                f.write(f'No data source pair has a significant percentage of shared molecules (greater than 80%).\n')

            # WARNING 6: Reference Molecules Underrepresentation Warning
            if self.__mainSelf.reference_set is not None:
                f.write(formatWarningTitle(text='Reference Molecules Underrepresentation'))
                if underrepresented_sources_df.shape[0] > 0:
                    f.write(f'Data sources with a low proportion of reference molecules (<15%):\n')
                    for i,row in underrepresented_sources_df.iterrows():
                        f.write(f'{i+1}) {row[config["NAMES"]["REF"]]} --> Percentage of reference molecules: {round(row["prop_ref_mols"], 2)}% ({row["num_ref_mols"]} molecules)\n')
                        if i+1 == underrepresented_sources_df.shape[0]:
                            f.write(f'\t***Figure reference: Comparison of {self.__mainSelf.endpoint_name} Distribution for Reference vs. Non-Reference Molecules Across Data Sources [Fig EXTRA].\n')
                else:
                    f.write(f'No data source has a low proportion of reference molecules (less than 15%).\n')

    ###
    def generateInsightReport(self, data_df, skewness_results, outliers_results, oor_results, 
                                feature_similarity_results, endpoint_distribution_results, discrepancies_df, ref_mols_results):
    
        """
        Conducts a final insight report to identify individual datasets:
        - with an endpoint distribution significantly differing from the rest (only for regression endpoints).
        - with an inconsistent endpoint value range (only for regression endpoints).
        - with a skewed endpoint distribution (only for regression endpoints).
        - with a significant outlier proportion (only for regression endpoints).
        - with abnormal data points (only for regression endpoints).
        - with out-of-range data points (only for regression endpoints).
        - which are dissimilar in terms of the molecular feature profile. 
        - which are conflicting due to differing endpoint value annotations for shared molecules. 
        - which are divergent due to a low molecule overlap of molecules with the other sources. 
        - which are redundant for sharing a significant number of molecules. 
        - which are underrepresented in reference molecules (only when a reference set is defined).
        """
        if self.__mainSelf.task == config["TASK_REGRESSION"]:
            # WARNING 1: Endpoint Distribution Warning
            std_dev = data_df[config["NAMES"]["VALUE"]].std() # Calculate standard deviation of the endpoint distribution using the entire dataset
            diff_distribution_df = self.__checkEndpointDistribution(endpoint_distribution_results, meandiff_thresh=std_dev)
            # WARNING 1.1: Value Range Warning
            diff_range_df = self.__checkValueRange(data_df)
            # WARNING 1.2: Skewed Distribution Warning
            skewed_distribution_df = self.__checkSkewness(skewness_results)
            # WARNING 1.3: Outliers Warning
            outliers_df = self.__checkOutliers(outliers_results)
            # WARNING 1.4: Abnormal Data Warning
            abnormal_data_df = self.__checkAbnormalDataPoints(outliers_results)
            # WARNING 1.5: Out-Of-Range Data Warning
            oor_data_df = self.__checkOORDataPoints(oor_results)

        elif self.__mainSelf.task == config["TASK_CLASSIFICATION"]:
            diff_distribution_df = diff_range_df = skewed_distribution_df = outliers_df = abnormal_data_df = oor_data_df = None

        # WARNING 2: Dissimilar Datasets Warning
        dissimilar_sources_df = self.__identifyDissimilarDatasets(feature_similarity_results)

        # WARNING 3: Conflicting Datasets Warning
        conflicting_sources_df = self.__identifyConflictingDatasets(discrepancies_df)

        # WARNING 4: Divergent Datasets Warning
        divergent_sources_df  = self.__identifyDivergentDatasets(data_df)

        # WARNING 5: Redundant Datasets Warning
        redundant_sources_df  = self.__identifyRedundantDatasets(discrepancies_df)

        # WARNING 6: Reference Molecules Underrepresentation Warning
        if self.__mainSelf.reference_set is not None:
            underrepresented_sources_df = self.__identifyUnderrepresentedDatasets(ref_mols_results)
        else:
            underrepresented_sources_df = None

        # Generate a TXT file with results
        self.__writeToTXT(diff_distribution_df, diff_range_df, skewed_distribution_df, outliers_df, abnormal_data_df, oor_data_df,
                        dissimilar_sources_df, conflicting_sources_df, divergent_sources_df, redundant_sources_df, underrepresented_sources_df)
