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

from scipy.stats import skew, skewtest, kurtosis, kurtosistest, chi2_contingency, ks_2samp

from sklearn.metrics import pairwise_distances

from umap import UMAP

from sklearn.preprocessing import MinMaxScaler
from itertools import combinations

from .AI_MoleculeInfo import MoleculeInfo

### Configs
import os
import json

config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'configs.json')

with open(config_path, "r") as file:
    config = json.load(file)

class Calculation():
    """
    Class to compute all parameters and statistics to report.
    """

    def __init__(self, mainSelf):
        self.__mainSelf = mainSelf

    ### 
    def getNmols(self, data_df):

        """
        Gets the total number of molecules in the dataset.
        """

        return [data_df.shape[0]]

    ### 
    def getNmols_perDataset(self, data_df):

        """
        Gets the number of molecules for each individual dataset.
        """

        n_mols_sources = []
        for ref in self.__mainSelf.sources_list:
            n_mols = self.getNmols(data_df.loc[data_df[config["NAMES"]["REF"]] == ref])
            n_mols_sources.append([ref]+n_mols)

        return pd.DataFrame(data=n_mols_sources, columns=['dataset','num_mols'])

    ### 
    def calculatePropRefMols(self, data_df, ref_set):

        """
        Calculates the number and proportion of reference molecules in the dataset.
        """

        # Get the subset of reference molecules in the dataset
        source_mols = set(data_df[config["NAMES"]["INCHIKEY"]].tolist())
        source_ref_mols = ref_set.intersection(source_mols)
        # Compute the number and percentage of reference molecules
        num_mols = len(source_ref_mols)
        prop_mols = (len(source_ref_mols) / len(source_mols)) * 100

        return [num_mols, prop_mols]

    ### 
    def calculatePropRefMols_perDataset(self, data_df, ref_set):

        """
        Calculates the number and proportion of reference molecules for each individual dataset.
        
        This type of analysis can only be conducted when a reference set of molecules is defined.
        """

        ref_mols_sources = []
        for ref in self.__mainSelf.sources_list:
            ref_mols = self.calculatePropRefMols(data_df.loc[data_df[config["NAMES"]["REF"]] == ref], ref_set)
            ref_mols_sources.append([ref]+ref_mols)

        return pd.DataFrame(data=ref_mols_sources, columns=['dataset','num_ref_mols','prop_ref_mols'])

    ###
    def calculateEndpointStatistics(self, data_df):

        """
        Calculates multiple statistics on the endpoint distribution.
        """

        # Compute endpoint statistics
        if self.__mainSelf.task == config["TASK_CLASSIFICATION"]:
            value_counts = ', '.join([f'{value}: {count}' for value, count in sorted(data_df[config["NAMES"]["VALUE"]].value_counts().items())])
            counts = [count for value, count in sorted(data_df[config["NAMES"]["VALUE"]].value_counts().items())]
            ratio = ':'.join(f'{count / np.min(counts):.1f}' for count in counts)
            
            return [value_counts, ratio]
        
        elif self.__mainSelf.task == config["TASK_REGRESSION"]:
            endpoint_description = data_df[config["NAMES"]["VALUE"]].describe()
            
            return endpoint_description.tolist()[1:]
        
    ### 
    def calculateEndpointStatistics_perDataset(self, data_df):

        """
        Calculates multiple statistics on the endpoint distribution for each individual dataset.
        """

        statistics_sources = []
        for ref in self.__mainSelf.sources_list:
            statisticts = self.calculateEndpointStatistics(data_df.loc[data_df[config["NAMES"]["REF"]] == ref])
            statistics_sources.append([ref]+statisticts)

        if self.__mainSelf.task == config["TASK_CLASSIFICATION"]:
            return pd.DataFrame(data=statistics_sources, columns=['dataset','class_counts','class_ratio'])
        elif self.__mainSelf.task == config["TASK_REGRESSION"]:
            return pd.DataFrame(data=statistics_sources, columns=['dataset','mean','standard_deviation','minimum',
                                                                '1st_quartile','median','3rd_quartile','maximum'])
        
    ### 
    def calculateSkewness(self, data_df):

        """
        Calculates the Fisher's Skewness value and tests whether the skew is different from the
        normal distribution.
        
        This type of analysis can only be conducted for regression endpoints.
        """

        # Compute skewness 
        skewness = skew(data_df[config["NAMES"]["VALUE"]].tolist(), axis=0, bias=True)

        # Skewness test
        skewness_test = skewtest(data_df[config["NAMES"]["VALUE"]].tolist(), axis=0)
        statistic = skewness_test[0]
        pvalue = skewness_test[1]

        return [skewness, statistic, pvalue]

    ### 
    def calculateSkewness_perDataset(self, data_df):

        """
        Calculates the Fisher's Skewness for each individual dataset and test whether their skew
        is different from the normal distribution.
        """

        skewness_sources = []
        for ref in self.__mainSelf.sources_list:
            skewness = self.calculateSkewness(data_df.loc[data_df[config["NAMES"]["REF"]] == ref])
            skewness_sources.append([ref]+skewness)

        return pd.DataFrame(data=skewness_sources, columns=['dataset','skewness_value','skewness_statistic','skewness_pvalue'])

    ### 
    def calculateKurtosis(self, data_df):

        """
        Calculates the Fisher's Kurtosis value and tests whether the dataset has normal kurtosis.
    
        This type of analysis can only be conducted for regression endpoints.
        """

        # Calculate kurtosis 
        kurtosis_value = kurtosis(data_df[config["NAMES"]["VALUE"]].tolist(), axis=0, fisher=True, bias=True)

        # Kurtosis test
        kurtosis_test = kurtosistest(data_df[config["NAMES"]["VALUE"]].tolist(), axis=0)
        statistic = kurtosis_test[0]
        pvalue = kurtosis_test[1]

        return [kurtosis_value, statistic, pvalue]

    ### 
    def calculateKurtosis_perDataset(self, data_df):

        """
        Calculates the Fisher's Kurtosis value for each individual dataset and test whether they
        have normal kurtosis. 
        """

        kurtosis_sources = []
        for ref in self.__mainSelf.sources_list:
            kurtosis_info = self.calculateKurtosis(data_df.loc[data_df[config["NAMES"]["REF"]] == ref])
            kurtosis_sources.append([ref]+kurtosis_info)

        return pd.DataFrame(data=kurtosis_sources, columns=['dataset','kurtosis_value','kurtosis_statistic','kurtosis_pvalue'])

    ### 
    def identifyOutliers(self, data_df):

        """
        Finds outliers by using two methods (assuming Normal distribution of the data):

        - Z-score computation: Calculates the z-score for each data point, and determines as
        outliers those with a z-score greater than +3 or less than -3 (more than 3 standard
        deviations away from the mean).

        - 1.5xIQR (interquartile range) rule: Calculates the IQR (Q3-Q1), defines both upper 
        (Q3+(1.5*IQR)) and lower (Q1-(1.5*IQR)) fences, and defines as outliers all data points
        outside these fences.

        In addition, we defined as abnormal those data points with a z-score more than 5 standard
        deviations away from the mean.

        The default method to define the outlier set for further analysis is Z-score calculation.
        
        This type of analysis can only be conducted for regression endpoints.
        """

        # 1st METHOD: Z-score computation
        zscore_df = data_df[[config["NAMES"]["INCHIKEY"], config["NAMES"]["ENDPOINT_ID"], config["NAMES"]["REF"], config["NAMES"]["VALUE"]]]
        zscore_df['zscore'] = zscore_df[config["NAMES"]["VALUE"]].apply(lambda x: (x - zscore_df[config["NAMES"]["VALUE"]].mean()) / zscore_df[config["NAMES"]["VALUE"]].std())
        zscore_outliers = zscore_df.loc[(zscore_df['zscore'] > 3) | (zscore_df['zscore'] < -3)]
        n_outliers_zscore = len(zscore_outliers)
        prop_outliers_zscore = (n_outliers_zscore / data_df.shape[0]) * 100

        # Identify abnormal data points 
        zscore_abnormals = zscore_df.loc[(zscore_df['zscore'] > 5) | (zscore_df['zscore'] < -5)]
        n_abnormals = len(zscore_abnormals)
        prop_abnormals = (n_abnormals / data_df.shape[0]) * 100

        # 2nd METHOD: 1.5xIQR rule
        iqr_df = data_df[[config["NAMES"]["INCHIKEY"], config["NAMES"]["ENDPOINT_ID"], config["NAMES"]["REF"], config["NAMES"]["VALUE"]]]
        iqr = (iqr_df[config["NAMES"]["VALUE"]].describe()['75%'] - iqr_df[config["NAMES"]["VALUE"]].describe()['25%'])
        upper_fence = iqr_df[config["NAMES"]["VALUE"]].describe()['75%'] + (1.5*iqr)
        lower_fence = iqr_df[config["NAMES"]["VALUE"]].describe()['25%'] - (1.5*iqr)
        iqr_outliers = iqr_df.loc[(iqr_df[config["NAMES"]["VALUE"]] > upper_fence) | (iqr_df[config["NAMES"]["VALUE"]] < lower_fence)]
        n_outliers_iqr = len(iqr_outliers)
        prop_outliers_iqr = (n_outliers_iqr / data_df.shape[0]) * 100

        # Define the outliers set
        if self.__mainSelf.outliers_method == 'zscore':
            outliers_set = set(zscore_outliers[config["NAMES"]["INCHIKEY"]].tolist())
        elif self.__mainSelf.outliers_method == 'iqr':
            outliers_set = set(iqr_outliers[config["NAMES"]["INCHIKEY"]].tolist())

        return [n_outliers_zscore, prop_outliers_zscore, n_abnormals, prop_abnormals, n_outliers_iqr, prop_outliers_iqr], outliers_set

    ### 
    def identifyOutliers_perDataset(self, data_df):

        """
        Finds outliers for each individual dataset.
        """

        outliers_sources = []
        outliers_dict = {}
        for ref in self.__mainSelf.sources_list:
            outliers_info, outliers_set = self.identifyOutliers(data_df.loc[data_df[config["NAMES"]["REF"]] == ref])
            outliers_sources.append([ref]+outliers_info)
            outliers_dict[ref] = outliers_set

        outliers_df = pd.DataFrame(data=outliers_sources, columns=['dataset','n_outliers_Zscore','prop_outliers_Zscore','n_abnormals_Zscore',
                                                                'prop_abnormals_Zscore','n_outliers_IQR','prop_outliers_IQR'])

        return outliers_df, outliers_dict

    ### 
    def identifyOORDataPoints(self, data_df):

        """
        Finds out-of-range (OOR) data points based on the predefined excpected bounds.

        This type of analysis can only be conducted for regression endpoints and when either or
        both bounds are defined.
        """

        upper_oor_counts = None
        lower_oor_counts = None

        # Find data points above the upper bound
        if self.__mainSelf.upper_bound is not None:
            upper_oor_data = data_df.loc[data_df[config["NAMES"]["VALUE"]] > self.__mainSelf.upper_bound]
            upper_oor_counts = len(upper_oor_data)

        # Find data points below the lower bound
        if self.__mainSelf.lower_bound is not None:
            lower_oor_data = data_df.loc[data_df[config["NAMES"]["VALUE"]] < self.__mainSelf.lower_bound]
            lower_oor_counts = len(lower_oor_data)

        return [upper_oor_counts, lower_oor_counts]

    ### 
    def identifyOORDataPoints_perDataset(self, data_df):

        """
        Finds out-of-range (OOR) data points for each individual dataset.
        """

        oor_data_sources = []
        for ref in self.__mainSelf.sources_list:
            oor_data = self.identifyOORDataPoints(data_df.loc[data_df[config["NAMES"]["REF"]] == ref])
            oor_data_sources.append([ref]+oor_data)

        return pd.DataFrame(data=oor_data_sources, columns=['dataset','n_upper_oor','n_lower_oor'])


    ### 
    def compareEndpointDistribution_acrossDatasets(self, data_df):

        """
        Conducts statistical tests to compare the endpoint distribution across individual datasets
        in a one-vs-others settings.

        For classification endpoints, the two-sample Chi-Square test is used; while for regression
        endpoints, the two-sample Kolmogorov-Smirnov test was selected.
        """

        # Compares statistically the endpoint distribution of each source vs. the others
        statistical_results_sources = []

        for ref in self.__mainSelf.sources_list:
            if self.__mainSelf.task == config["TASK_CLASSIFICATION"]:
                source_endpoint_value_counts = data_df[config["NAMES"]["VALUE"]].loc[data_df[config["NAMES"]["REF"]] == ref].value_counts().sort_index(ascending=True)
                others_endpoint_value_counts = data_df[config["NAMES"]["VALUE"]].loc[data_df[config["NAMES"]["REF"]] != ref].value_counts().sort_index(ascending=True)
                # Conduct Two-sample Chi-Square test
                chi2_statistic, p_value, dof, expected_freq = chi2_contingency(np.array([source_endpoint_value_counts, others_endpoint_value_counts]), 
                                                                            correction=False)
                statistical_results_sources.append([ref, chi2_statistic, p_value, None])
            
            elif self.__mainSelf.task == config["TASK_REGRESSION"]:
                source_endpoint_values = data_df[config["NAMES"]["VALUE"]].loc[data_df[config["NAMES"]["REF"]] == ref].tolist()
                others_endpoint_values = data_df[config["NAMES"]["VALUE"]].loc[data_df[config["NAMES"]["REF"]] != ref].tolist()
                # Conduct Two-sample Kolmogorov–Smirnov test
                ks_statistic, p_value = ks_2samp(source_endpoint_values, others_endpoint_values)
                # Calculate the mean difference between both distributions
                meandiff = np.mean(source_endpoint_values) - np.mean(others_endpoint_values)
                statistical_results_sources.append([ref, ks_statistic, p_value, meandiff])
        
        return pd.DataFrame(data=statistical_results_sources, columns=['dataset','statistic','p_value', 'meandiff'])

    ### 
    def calculateDistanceMatrix(self, data_df):

        """
        Calculates the molecule distance matrix for the entire dataset. 

        When the feature profile is based on RDKit descriptors, the standardized Euclidean distance
        is computed; while then the molecule profile is based on ECFP4 fingerprints, the Tanimoto 
        similarity coefficient (or Jaccard similarity coefficient) is calculated.
        """

        # Compute the distance matrix
        if self.__mainSelf.feature_type in ['rdkit']: 
            variances = np.var(data_df.loc[:, [f'{self.__mainSelf.features}_{feature}' for feature in MoleculeInfo.AVAILABLE_FEATURES[self.__mainSelf.features]]].to_numpy(), axis=0, ddof=1) # compute feature variances
            variances[variances == 0] = 1e-8 # replace zero variances with a small positive value to avoid division by zero
            distance_matrix = pairwise_distances(data_df.loc[:, [f'{self.__mainSelf.features}_{feature}' for feature in MoleculeInfo.AVAILABLE_FEATURES[self.__mainSelf.features]]].to_numpy(), metric='seuclidean', V=variances)

        elif self.__mainSelf.feature_type in ['ecfp4']: 
            indices = pairwise_distances(data_df.loc[:, [f'{self.__mainSelf.features}_{feature}' for feature in MoleculeInfo.AVAILABLE_FEATURES[self.__mainSelf.features]]].to_numpy(dtype=bool), metric='jaccard') 
            distance_matrix = 1 - indices

        elif self.__mainSelf.feature_type == 'custom': 
            distance_matrix = pairwise_distances(data_df.loc[:, [col for col in data_df.columns if col.startswith(self.__mainSelf.feature_type.upper())]].to_numpy(), metric=self.__mainSelf.distance_metric) 
        return distance_matrix

    ### 
    def calculateFeatureSimilarity(self, distance_matrix, data_df):

        """
        Calculates the distribution of intra- and inter-similarity values across individual datasets
        in a one-vs-others setting. 
        """

        # Compute the distribution of intra- and inter-similarity scores for each data source
        similarity_results = []

        for ref in self.__mainSelf.sources_list:
            source_indices = data_df.loc[data_df[config["NAMES"]["REF"]] == ref].index.to_list()
            other_sources_indices = data_df.loc[data_df[config["NAMES"]["REF"]] != ref].index.to_list()

            # Iterate over each molecule in the current source and compute the minimum within- 
            # and between-source similarities
            within_similarities = []
            between_similarities = []

            for idx in source_indices:
                # INTRA-SIMILARITY
                intra_similarity = np.min([distance_matrix[idx][i] for i in source_indices if i != idx]) if len(source_indices) > 1 else np.nan
                within_similarities.append(intra_similarity)
                # INTER-SIMILARITY
                inter_similarity = np.min([distance_matrix[idx][i] for i in other_sources_indices]) if len(other_sources_indices) > 0 else np.nan
                between_similarities.append(inter_similarity)

            similarity_results.append([ref, within_similarities, between_similarities])         

        return similarity_results 

    ### 
    def calculateFeatureSimilarity_Pairwise(self, distance_matrix, data_df):

        """
        Calculates the mean of the distribution of inter-similarity values for each dataset pair. 
        """

        # Compute the mean of the inter-source similarity scores for each dataset pair
        similarity_results = []

        for ref1 in self.__mainSelf.sources_list:
            for ref2 in self.__mainSelf.sources_list:
                if ref1 != ref2:
                    ref1_indices = data_df.loc[data_df[config["NAMES"]["REF"]] == ref1].index.to_list()
                    ref2_indices = data_df.loc[data_df[config["NAMES"]["REF"]] == ref2].index.to_list()

                    # Iterate over each molecule in the reference source (ref1) and compute the 
                    # minimum between-source similarity
                    similarities = []

                    for idx in ref1_indices:
                        min_similarity = np.min([distance_matrix[idx][i] for i in ref2_indices]) if len(ref2_indices) > 0 else np.nan
                        similarities.append(min_similarity)

                    similarity_results.append([ref1, ref2, np.mean(similarities)])         
        
        return pd.DataFrame(similarity_results, columns=['ref1','ref2','mean_inter_similarity']) 

    ###
    def runUMAP(self, data_df):

        """
        Runs the Uniform Manifold Approximation and Projection (UMAP) dimensionality reduction technique.
        
        When the molecules are defined by RDKit descriptors, the standardized Euclidean distance
        is selected, while when ECFP4 fingerprints are used to define molecules, the Jaccard distance
        is chosen.
        """

        # Run UMAP
        if self.__mainSelf.feature_type in ['rdkit']:
            variances = np.var(data_df.loc[:, [f'{self.__mainSelf.features}_{feature}' for feature in MoleculeInfo.AVAILABLE_FEATURES[self.__mainSelf.features]]].to_numpy(), axis=0, ddof=1) # compute feature variances
            variances[variances == 0] = 1e-8 # replace zero variances with a small positive value to avoid division by zero
            umap = UMAP(n_components=2, init='random', metric='seuclidean', metric_kwds={'V': variances}, random_state=config["RANDOM_SEED"], n_jobs=1)
        
        elif self.__mainSelf.feature_type in ['ecfp4']:
            umap = UMAP(n_components=2, init='random', metric='jaccard', random_state=config["RANDOM_SEED"], n_jobs=1)

        elif self.__mainSelf.feature_type == 'custom':
            umap = UMAP(n_components=2, init='random', metric=self.__mainSelf.distance_metric, random_state=config["RANDOM_SEED"], n_jobs=1)

        if self.__mainSelf.feature_type in ['rdkit','ecfp4']:
            proj = umap.fit_transform(data_df.loc[:, [f'{self.__mainSelf.features}_{feature}' for feature in MoleculeInfo.AVAILABLE_FEATURES[self.__mainSelf.features]]])
        elif self.__mainSelf.feature_type in ['custom']:
            proj = umap.fit_transform(data_df.loc[:, [col for col in data_df.columns if col.startswith(self.__mainSelf.feature_type.upper())]])

        return proj

    ###
    def calculateInterDatasetDiscrepancies(self, data_df):

        """
        Identifies the discrepencies between individual datasets by computing: (1) the number of 
        common molecules, and (2) the differences in endpoint values of these shared molecules.
        
        For classification endpoints, the number of different annotated values is reported, while
        the (standardized) mean absolute difference is computed for regression endpoints.
        """

        # Crate a copy of the DataFrame
        data_df_copy = data_df.copy()

        # Standardize the endpoint values (only for regression endpoints)
        if self.__mainSelf.task == config["TASK_REGRESSION"]:
            scaler = MinMaxScaler()
            data_df_copy[[config["NAMES"]["VALUE"]]] = scaler.fit_transform(data_df_copy[[config["NAMES"]["VALUE"]]])
            
        # Get the molecule set for each dataset
        mols_per_source = {}
        for ref in self.__mainSelf.sources_list:
            mols_per_source[ref] = data_df_copy[[config["NAMES"]["INCHIKEY"], config["NAMES"]["VALUE"]]].loc[data_df_copy[config["NAMES"]["REF"]] == ref]

        # Calculate the number of common molecules and the endpoint value differences across source pairs (pairwise)
        data_info = []
        for ref1, ref2 in combinations(self.__mainSelf.sources_list, 2):
            mols_ref1 = mols_per_source[ref1].set_index(config["NAMES"]["INCHIKEY"])
            mols_ref2 = mols_per_source[ref2].set_index(config["NAMES"]["INCHIKEY"])

            # Find common molecules
            common_mols = mols_ref1.index.intersection(mols_ref2.index)

            if len(common_mols) > 0:
                if self.__mainSelf.task == config["TASK_CLASSIFICATION"]:
                    # Calculate the number of molecules with differing annotated endpoint value
                    diff_metric = (mols_ref1.loc[common_mols][config["NAMES"]["VALUE"]] != mols_ref2.loc[common_mols][config["NAMES"]["VALUE"]]).sum()
                    diff_tokens = [diff_metric, diff_metric/len(common_mols)]  

                elif self.__mainSelf.task == config["TASK_REGRESSION"]:
                    # Calculate the (standardized) endpoint absolute difference 
                    diff_metric = np.abs(mols_ref1.loc[common_mols][config["NAMES"]["VALUE"]] - mols_ref2.loc[common_mols][config["NAMES"]["VALUE"]])
                    diff_tokens = [np.mean(diff_metric), np.median(diff_metric)]
            
            else:
                diff_tokens = [np.nan, np.nan]

            data_info.append([ref1, len(mols_ref1), ref2, len(common_mols), (len(common_mols)/len(mols_ref1))*100] + diff_tokens)
            data_info.append([ref2, len(mols_ref2), ref1, len(common_mols), (len(common_mols)/len(mols_ref2))*100] + diff_tokens)          

        # Convert into a DataFrame
        if self.__mainSelf.task == config["TASK_CLASSIFICATION"]:
            info_df = pd.DataFrame(data=data_info, columns=['source_1', 'num_mols', 'source_2', 'num_common_mols', 'percent_common_mols', 'different_values_count', 'different_values_count_proportion'])
        elif self.__mainSelf.task == config["TASK_REGRESSION"]:
            info_df = pd.DataFrame(data=data_info, columns=['source_1', 'num_mols', 'source_2', 'num_common_mols', 'percent_common_mols', 'mean_abs_diff', 'median_abs_diff'])

        return info_df

    ### 
    def perfromPairwiseKSTest(self, data_df):

        """
        Performs the two-sample Kolmogorov-Smirnov (KS) test to compare the endpoint distribution across 
        all pairs of individual datasets.
        
        This type of analysis can only be conducted for regression endpoints.
        """

        if self.__mainSelf.task == config["TASK_CLASSIFICATION"]:
            return None

        elif self.__mainSelf.task == config["TASK_REGRESSION"]:
            # Perform pairwise comparisons using the two-sample Kolmogorov-Smirnov (KS) test
            ks_results = []
            for ref1, ref2 in combinations(self.__mainSelf.sources_list, 2):  
                        statistic, p_value = ks_2samp(data_df[config["NAMES"]["VALUE"]].loc[data_df[config["NAMES"]["REF"]] == ref1].tolist(), data_df[config["NAMES"]["VALUE"]].loc[data_df[config["NAMES"]["REF"]] == ref2].tolist())
                        # Calculate the absolute mean difference
                        meandiff = np.abs(np.mean(data_df[config["NAMES"]["VALUE"]].loc[data_df[config["NAMES"]["REF"]] == ref1].tolist()) - np.mean(data_df[config["NAMES"]["VALUE"]].loc[data_df[config["NAMES"]["REF"]] == ref2].tolist()))
                        ks_results.append([ref1, ref2, meandiff, p_value])
            # Generate a DataFrame with the KS test results
            ks_results_df = pd.DataFrame(ks_results, columns=['group1', 'group2', 'meandiff', 'p_value'])
            # Add the 'reject' column
            ks_results_df['reject'] = ks_results_df['p_value'].apply(lambda x: True if x < 0.05 else False)

            return ks_results_df
