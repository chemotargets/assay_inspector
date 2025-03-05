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

from DR_DataLoading import *
from DR_Calculation import *
from DR_Visualization import *
from DR_OutputFile import *
from DR_ExpertAssessment import *

from DR_MoleculeData import MoleculeData
from DR_MoleculeInfo import MoleculeInfo

from DR_Utils import logging

import datetime
import seaborn as sns

import os
import json

### Configs
with open("configs.json", "r") as file:
    config = json.load(file)

###
class DataReporting():

    """
    Class to generate the reporting of individual and aggregated datasets build from multiple
    sources. 
    
    NOTE: This class is specific to datasets containing molecules (i.e., molecule-based datasets).
    """
    ### Define constants

    ###
    def __init__(self, data_path, endpoint_name, task, feature_type, outliers_method='zscore',
                 reference_set=None, lower_bound=None, upper_bound=None):

        """
        Class constructor. Requires input data.
        """

        self.data_path = data_path
        self.endpoint_name = endpoint_name
        self.task = task.upper()
        self.feature_type = feature_type

        self.outliers_method = outliers_method

        self.reference_set = reference_set

        self.lower_bound = lower_bound
        self.upper_bound = upper_bound

        features_dict = {'rdkit':MoleculeInfo.FEAT_RDKIT_DESC, 'ecfp4':MoleculeInfo.FEAT_ECFP4}
        self.features = features_dict[self.feature_type]

        self.directory = 'DataReporting_'+datetime.now().strftime("%Y%m%d")

        # Define the color palette 
        color_palette = sns.color_palette('colorblind')
        self.hex_colors = color_palette.as_hex()

    ### 
    def get_individual_reporting(self):

        """
        Generates a report for an individual dataset.
        
        It loads the input data, compute various statistics, and generates multiple plots for
        analysis. It then creates an output file summarizing the key information. 
        
        NOTE: Behavior varies depending on the task type (classification or regression).
        """

        # Load data
        if not os.path.exists(self.data_path):
            logging.error(f"The file {self.data_path} does not exist")
            return 
        data_instance = MoleculeData(source=self.data_path)

        # Generate the DataFrame of molecule feature profiles
        data = self.__getMoleculeProfileDataframe(data_instance, reporting='individual')

        if self.reference_set is not None:
            # Load reference set data
            if not os.path.exists(self.reference_set):
                logging.error(f"The file {self.reference_set} does not exist")
                return 
            ref_instance = MoleculeData(source=self.reference_set)

            # Generate the InChIKey set of reference molecules
            ref_inchikeys = self.__getInChIKeySet(ref_instance)

        # Check data type of the endpoint
        if self.task == config["TASK_CLASSIFICATION"]:
            endpoint_classes = data[config["NAMES"]["VALUE"]].unique()
            if len(endpoint_classes) != 2: 
                logging.error(f"The number of endpoint classes is {len(endpoint_classes)} ({', '.join(endpoint_classes)}), but 2 were expected")

        elif self.task == config["TASK_REGRESSION"]:
            if not pd.api.types.is_numeric_dtype(data[config["NAMES"]["VALUE"]]):
                logging.error(f"The data type of the endpoint value is not numeric but {data[config['NAMES']['VALUE']].dtype}")

        # Count the number of molecules
        n_mols = self.__getNmols(data)

        # Compute endpoint distribution statistics
        statistics = self.__calculateEndpointStatistics(data)

        if self.reference_set is not None:
            # Calculate the percentage of reference molecules
            prop_ref_mols = self.__calculatePropRefMols(data, ref_inchikeys)

        outliers_set = None
        if self.task == config["TASK_REGRESSION"]:
            # Compute skewness and kurtosis
            skewness = self.__calculateSkewness(data)
            kurtosis_df = self.__calculateKurtosis(data)
            # Identify outliers and out-of-range (OOR) data points
            outliers_info, outliers_set = self.__identifyOutliers(data)
            oor_data = self.__identifyOORDataPoints(data)
       
        # Create the main directory and the endpoint subdirectory
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
        if not os.path.exists(os.path.join(self.directory, self.endpoint_name)):
            os.makedirs(os.path.join(self.directory, self.endpoint_name))

        # Generate the final DataFrame and save it into a TSV file
        logging.info(oriMessage= f"Creating report on {self.endpoint_name} individual dataset")

        if self.task == config["TASK_CLASSIFICATION"]:
            summary_df = pd.DataFrame(data=[[self.endpoint_name , 'entire_dataset'] + n_mols + statistics],
                                      columns=[config["NAMES"]["ENDPOINT_ID"],'dataset','num_mols','class_counts','class_ratio'])
        elif self.task == config["TASK_REGRESSION"]:
            summary_df = pd.DataFrame(data=[[self.endpoint_name, 'entire_dataset'] + n_mols + statistics + skewness + kurtosis_df + outliers_info + oor_data],
                                      columns=[config["NAMES"]["ENDPOINT_ID"],'dataset','num_mols','mean','standard_deviation','minimum','1st_quartile','median','3rd_quartile','maximum',
                                              'skewness_value','skewness_statistic','skewness_pvalue','kurtosis_value','kurtosis_statistic','kurtosis_pvalue','n_outliers_Zscore','prop_outliers_Zscore',
                                              'n_abnormals_Zscore','prop_abnormals_Zscore','n_outliers_IQR','prop_outliers_IQR','n_upper_oor','n_lower_oor'])
            
        if self.reference_set is not None:
            summary_df = pd.merge(summary_df, pd.DataFrame(data=[['entire_dataset']+prop_ref_mols], columns=['dataset','num_ref_mols','prop_ref_mols']), on='dataset', how='left')
        
        self.__writeToTSV(summary_df)

       
        if self.task == config["TASK_REGRESSION"]:
            # Visualize outliers
            self.__VisualizeOutliers(data, outliers_set)

        # Plot the endpoint distribution
        self.__plotEndpointDistribution(data, outliers_set)

        # Compute the distance matirx and Plot the similarity distribution
        distance_matrix = self.__calculateDistanceMatrix(data)
        self.__plotSimilarityDistribution(distance_matrix)

        # Run UMAP and Visualize the feature space
        projection = self.__runUMAP(data)
        self.x_range, self.y_range = self.__getAxisRanges(projection)
        self.__plotFeatureSpace_coloredbyEndpoint(projection, data)
        self.__plotFeatureSpace_Hexbin(projection)    

        logging.info(oriMessage= f"The final report and several plots have been saved in the {self.directory+'/'+self.endpoint_name} directory")

    ### 
    def get_comparative_reporting(self):

        """
        Generates a report for an aggregated dataset build from multiple sources.

        It loads the input data, compute various statistics, and generates multiple plots for
        analysis. It then creates an output file summarizing the key information and an expert
        assessment file with alerts to guide data cleaning and preprocessing. 
        
        NOTE: Behavior varies depending on the task type (classification or regression).
        """

        # Load data
        if not os.path.exists(self.data_path):
            logging.error(f"The file {self.data_path} does not exist")
            return 
        data_instance = MoleculeData(source=self.data_path)

        # Generate the DataFrame of molecule feature profiles
        data = self.__getMoleculeProfileDataframe(data_instance, reporting='comparative')
        
        if self.reference_set is not None:
            # Load reference set data
            if not os.path.exists(self.reference_set):
                logging.error(f"The file {self.reference_set} does not exist")
                return 
            ref_instance = MoleculeData(source=self.reference_set)

            # Generate the InChIKey set of reference molecules
            ref_inchikeys = self.__getInChIKeySet(ref_instance)

        # Check data type of the endpoint value
        if self.task == config["TASK_CLASSIFICATION"]:
            endpoint_classes = data[config["NAMES"]["VALUE"]].unique()
            if len(endpoint_classes) != 2:
                logging.error(f"The number of endpoint classes is {len(endpoint_classes)} ({', '.join(endpoint_classes)}), but 2 were expected")

        elif self.task == config["TASK_REGRESSION"]:
            if not pd.api.types.is_numeric_dtype(data[config["NAMES"]["VALUE"]]):
                logging.error(f"The data type of the endpoint value is not numeric but {data[config['NAMES']['VALUE']].dtype}")

        # Get the list of sources and sort them by counts
        sources_list = data[config["NAMES"]["REF"]].unique().tolist()
        sources_counts = {ref: len(data.loc[data[config["NAMES"]["REF"]] == ref]) for ref in sources_list}
        self.sources_list = sorted(sources_list, key=lambda ref: sources_counts[ref], reverse=True)

        # Verify the presence of multiple data sources
        if len(self.sources_list) == 1:
            logging.error(f"Expected multiple data sources, but only one was provided. Use the `.get_individual_reporting()` method instead.")
            return

        # Count the total number of molecules
        n_mols_total = self.__getNmols(data)

        # Count the number of molecules per data source
        n_mols_sources = self.__getNmols_perDataset(data)

        # Compute endpoint distribution statistics on the entire dataset
        statistics_entire_dataset = self.__calculateEndpointStatistics(data)

        # Compute endpoint distribution statistics per data source
        statistics_sources = self.__calculateEndpointStatistics_perDataset(data)

        if self.reference_set is not None:
            # Calculate the total percentage of reference molecules
            prop_ref_mols_total = self.__calculatePropRefMols(data, ref_inchikeys)

            # Calculate the percentage of reference molecules per data source
            prop_ref_mols_sources= self.__calculatePropRefMols_perDataset(data, ref_inchikeys)

        else: 
            prop_ref_mols_total = prop_ref_mols_sources = None

        if self.task == config["TASK_REGRESSION"]:
            # Compute skewness and kurtosis on the entire dataset
            skewness_entire_dataset = self.__calculateSkewness(data)
            kurtosis_entire_dataset = self.__calculateKurtosis(data)
            # Identify outliers and OOR data points on the entire dataset
            outliers_info_entire_dataset, outliers_set_entire_dataset = self.__identifyOutliers(data)
            oor_entire_dataset = self.__identifyOORDataPoints(data)

            # Compute skewness and kurtosis per data sources
            skewness_sources = self.__calculateSkewness_perDataset(data)
            kurtosis_sources = self.__calculateKurtosis_perDataset(data)
            # Identify outliers and OOR data points per data source
            outliers_info_sources, outliers_set_sources = self.__identifyOutliers_perDataset(data)
            oor_sources = self.__identifyOORDataPoints_perDataset(data)

        else: 
            skewness_entire_dataset = kurtosis_entire_dataset = outliers_info_entire_dataset = outliers_set_entire_dataset = oor_entire_dataset = None
            skewness_sources = kurtosis_sources = outliers_info_sources = outliers_set_sources = oor_sources = None

        # Compare statistically the endpoint distribution across data sources
        endpoint_distribution_results = self.__compareEndpointDistribution_acrossDatasets(data)

        # Perform Pairwise KS Test
        ks_results = self.__perfromPairwiseKSTest(data)

        # Compute the distance matirx and get the within- and betwen-source distances per data source
        distance_matrix = self.__calculateDistanceMatrix(data)
        feature_similarity_results = self.__calculateFeatureSimilarity(distance_matrix, data)
        feature_similarity_pairwise_results = self.__calculateFeatureSimilarityPairwise(distance_matrix, data)

        # Create the main directory and the endpoint subdirectory
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
        if not os.path.exists(os.path.join(self.directory, self.endpoint_name)):
            os.makedirs(os.path.join(self.directory, self.endpoint_name))

        # Generate the final DataFrame and save it into a TSV file
        logging.info(f"Creating comparative report on {self.endpoint_name} aggregated dataset")

        summary_df = self.__getSummaryDataFrame(n_mols_total, n_mols_sources, statistics_entire_dataset, statistics_sources, prop_ref_mols_total, prop_ref_mols_sources, 
                                                skewness_entire_dataset, skewness_sources, kurtosis_entire_dataset, kurtosis_sources, outliers_info_entire_dataset, 
                                                outliers_info_sources, oor_entire_dataset, oor_sources, endpoint_distribution_results, feature_similarity_results)
        self.__writeToTSV(summary_df)

        if self.task == config["TASK_REGRESSION"]:
            # Visualize outliers
            self.__VisualizeOutliers(data, outliers_set_entire_dataset) 

        # Plot the intersection across data sources
        self.__plotDatasetsIntersection(data)

        # Plot the endpoint distribution
        self.__plotEndpointDistribution(data, outliers_set_entire_dataset)

        # Plot the endpoint distribution for each data source
        self.__plotEndpointDistributionComparison(data, outliers_set_sources)

        if self.reference_set is not None:
            # Plot the endpoint distribution for reference vs. non-reference molecules for each data source
            self.__plotEndpointDistributionComparison_RefVsNonRefMols(data, ref_inchikeys)

        # Inspect the discrepancies between data sources and Plot the comparative heatmaps
        discrepancies_df = self.__calculateInterDatasetDiscrepancies(data)
        self.__plotComparativeHeatmaps(discrepancies_df)
    
        # Plot KS test results in a heatmap
        self.__plotPairwiseKSTestHeatmap(ks_results)

        # Plot the similarity distribution
        self.__plotSimilarityDistribution(distance_matrix)
        self.__plotFeatureSimilarityHeatmap(feature_similarity_pairwise_results)

        # Run UMAP and Visualize the feature space
        projection = self.__runUMAP(data)
        self.x_range, self.y_range = self.__getAxisRanges(projection)
        self.__plotFeatureSpace_coloredbyEndpoint(projection, data)
        self.__plotFeatureSpace_coloredbyDataset(projection, data)
        self.__plotFeatureSpace_KDEplot(projection, data)
        self.__plotFeatureSpace_Hexbin(projection)    

        # Generate the final expert assessment
        self.__generateExpertAssessment(data, skewness_sources, outliers_info_sources, oor_sources, feature_similarity_results, 
                                        endpoint_distribution_results, discrepancies_df, prop_ref_mols_sources)

        logging.info(f"The final report and several plots have been saved in the {os.path.join(self.directory, self.endpoint_name)} directory")

        self.__generateOutputSummary(summary_df)

    ###
    @property
    def columns(self):
        return self.__columns
    
    @columns.setter
    def columns(self, value):
        self.__columns = value
