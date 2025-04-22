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

from AI_DataLoading import DataLoading
from AI_Calculation import Calculation
from AI_OutputFile import OutputFile
from AI_Visualization import Visualization
from AI_InsightReport import InsightReport

from AI_MoleculeData import MoleculeData
from AI_MoleculeInfo import MoleculeInfo

from AI_Utils import logging

import pandas as pd
import seaborn as sns
from datetime import datetime

### Warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

### Configs
import os
import json

config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'configs.json')

with open(config_path, "r") as file:
    config = json.load(file)

###
class AssayInspector():

    """
    Class to generate the reporting of individual and aggregated datasets build from multiple
    sources. 
    
    NOTE: This class is specific to datasets containing molecules (i.e., molecule-based datasets).
    """
    ### Define constants

    ###
    def __init__(self, data, endpoint_name, task, feature_type, outliers_method='zscore', distance_metric='euclidean', 
                 descriptors_df=None, reference_set=None, lower_bound=None, upper_bound=None):

        """
        Class constructor. Requires input data.
        """

        self.data = data
        self.endpoint_name = endpoint_name
        self.task = task.upper()
        self.feature_type = feature_type

        self.outliers_method = outliers_method
        self.distance_metric = distance_metric

        self.descriptors_df = None if descriptors_df is None else descriptors_df.copy()
        self.reference_set = reference_set

        self.lower_bound = lower_bound
        self.upper_bound = upper_bound

        if self.feature_type in ['rdkit', 'ecfp4']:
            features_dict = {'rdkit':MoleculeInfo.FEAT_RDKIT_DESC, 'ecfp4':MoleculeInfo.FEAT_ECFP4}
            self.features = features_dict[self.feature_type]
        elif self.feature_type == 'custom':
            self.features = self.descriptors_df
        else:
            logging.error(f"Feature type '{self.feature_type}' not implemented")

        self.directory = 'AssayInspector_'+datetime.now().strftime("%Y%m%d")

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

        self.__dataloading = DataLoading(mainSelf=self)
        self.__calculation = Calculation(mainSelf=self)
        self.__outputfile = OutputFile(mainSelf=self)
        self.__visualization = Visualization(mainSelf=self)

        # Load data
        if isinstance(self.data, str) and not os.path.exists(self.data):
            logging.error(f"The file {self.data} does not exist")
            return 
        data_instance = MoleculeData(source=self.data)

        # Generate the DataFrame of molecule feature profiles
        data = self.__dataloading.getMoleculeProfileDataframe(data_instance, reporting='individual')

        if self.reference_set is not None:
            # Load reference set data
            if not os.path.exists(self.reference_set):
                logging.error(f"The file {self.reference_set} does not exist")
                return 
            ref_instance = MoleculeData(source=self.reference_set)

            # Generate the InChIKey set of reference molecules
            ref_inchikeys = self.__dataloading.getInChIKeySet(ref_instance)

        # Check data type of the endpoint
        if self.task == config["TASK_CLASSIFICATION"]:
            endpoint_classes = data[config["NAMES"]["VALUE"]].unique()
            if len(endpoint_classes) != 2: 
                logging.error(f"The number of endpoint classes is {len(endpoint_classes)} ({', '.join(endpoint_classes)}), but 2 were expected")

        elif self.task == config["TASK_REGRESSION"]:
            if not pd.api.types.is_numeric_dtype(data[config["NAMES"]["VALUE"]]):
                logging.error(f"The data type of the endpoint value is not numeric but {data[config['NAMES']['VALUE']].dtype}")

        # Count the number of molecules
        n_mols = self.__calculation.getNmols(data)

        # Compute endpoint distribution statistics
        statistics = self.__calculation.calculateEndpointStatistics(data)

        if self.reference_set is not None:
            # Calculate the percentage of reference molecules
            prop_ref_mols = self.__calculation.calculatePropRefMols(data, ref_inchikeys)

        outliers_set = None
        if self.task == config["TASK_REGRESSION"]:
            # Compute skewness and kurtosis
            skewness = self.__calculation.calculateSkewness(data)
            kurtosis_df = self.__calculation.calculateKurtosis(data)
            # Identify outliers and out-of-range (OOR) data points
            outliers_info, outliers_set = self.__calculation.identifyOutliers(data)
            oor_data = self.__calculation.identifyOORDataPoints(data)
       
        # Create the main directory and the endpoint subdirectory
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
        if not os.path.exists(os.path.join(self.directory, self.endpoint_name)):
            os.makedirs(os.path.join(self.directory, self.endpoint_name))

        # Generate the final DataFrame and save it into a TSV file
        logging.info(f"Creating report on {self.endpoint_name} individual dataset")

        if self.task == config["TASK_CLASSIFICATION"]:
            summary_df = pd.DataFrame(data=[[self.endpoint_name , 'entire_dataset'] + n_mols + statistics],
                                      columns=[config["NAMES"]["ENDPOINT_ID"],'dataset','num_mols','class_counts','class_ratio'])
        elif self.task == config["TASK_REGRESSION"]:
            summary_df = pd.DataFrame(data=[[self.endpoint_name, 'entire_dataset'] + n_mols + statistics + skewness + kurtosis_df + outliers_info + oor_data],
                                      columns=[config["NAMES"]["ENDPOINT_ID"],'dataset','num_mols','mean','standard_deviation','minimum','1st_quartile','median','3rd_quartile','maximum',
                                              'skewness_value','skewness_statistic','skewness_pvalue','kurtosis_value','kurtosis_statistic','kurtosis_pvalue','n_outliers_Zscore','prop_outliers_Zscore',
                                              'n_abnormals_Zscore','prop_abnormals_Zscore','n_outliers_IQR','prop_outliers_IQR','n_upper_oor','n_lower_oor'])
            if self.upper_bound is None:
                summary_df = summary_df.drop(columns=['n_upper_oor'])
            if self.lower_bound is None:
                summary_df = summary_df.drop(columns=['n_lower_oor'])
            
        if self.reference_set is not None:
            summary_df = pd.merge(summary_df, pd.DataFrame(data=[['entire_dataset']+prop_ref_mols], columns=['dataset','num_ref_mols','prop_ref_mols']), on='dataset', how='left')
        
        self.__outputfile.writeToTSV(summary_df)

       
        if self.task == config["TASK_REGRESSION"]:
            # Visualize outliers
            self.__visualization.VisualizeOutliers(data, outliers_set)

        # Plot the endpoint distribution
        self.__visualization.plotEndpointDistribution(data, outliers_set)

        # Compute the distance matirx and Plot the similarity distribution
        distance_matrix = self.__calculation.calculateDistanceMatrix(data)
        self.__visualization.plotSimilarityDistribution(distance_matrix)

        # Run UMAP and Visualize the feature space
        projection = self.__calculation.runUMAP(data)
        self.x_range, self.y_range = self.__visualization.getAxisRanges(projection)
        self.__visualization.plotFeatureSpace_coloredbyEndpoint(projection, data)
        self.__visualization.plotFeatureSpace_Hexbin(projection)    

        logging.info(f"The final report and several plots have been saved in the {self.directory+'/'+self.endpoint_name} directory")

    ### 
    def get_comparative_reporting(self):

        """
        Generates a report for an aggregated dataset build from multiple sources.

        It loads the input data, compute various statistics, and generates multiple plots for
        analysis. It then creates an output file summarizing the key information and an insight
        report file with alerts to guide data cleaning and preprocessing. 
        
        NOTE: Behavior varies depending on the task type (classification or regression).
        """

        self.__dataloading = DataLoading(mainSelf=self)
        self.__calculation = Calculation(mainSelf=self)
        self.__outputfile = OutputFile(mainSelf=self)
        self.__visualization = Visualization(mainSelf=self)
        self.__insightreport = InsightReport(mainSelf=self)

        # Load data
        if isinstance(self.data, str) and not os.path.exists(self.data):
            logging.error(f"The file {self.data} does not exist")
            return 
        data_instance = MoleculeData(source=self.data)

        # Generate the DataFrame of molecule feature profiles
        data = self.__dataloading.getMoleculeProfileDataframe(data_instance, reporting='comparative')

        if self.reference_set is not None:
            # Load reference set data
            if not os.path.exists(self.reference_set):
                logging.error(f"The file {self.reference_set} does not exist")
                return 
            ref_instance = MoleculeData(source=self.reference_set)

            # Generate the InChIKey set of reference molecules
            ref_inchikeys = self.__dataloading.getInChIKeySet(ref_instance)

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
        n_mols_total = self.__calculation.getNmols(data)

        # Count the number of molecules per data source
        n_mols_sources = self.__calculation.getNmols_perDataset(data)

        # Compute endpoint distribution statistics on the entire dataset
        statistics_entire_dataset = self.__calculation.calculateEndpointStatistics(data)

        # Compute endpoint distribution statistics per data source
        statistics_sources = self.__calculation.calculateEndpointStatistics_perDataset(data)

        if self.reference_set is not None:
            # Calculate the total percentage of reference molecules
            prop_ref_mols_total = self.__calculation.calculatePropRefMols(data, ref_inchikeys)

            # Calculate the percentage of reference molecules per data source
            prop_ref_mols_sources= self.__calculation.calculatePropRefMols_perDataset(data, ref_inchikeys)

        else: 
            prop_ref_mols_total = prop_ref_mols_sources = None

        if self.task == config["TASK_REGRESSION"]:
            # Compute skewness and kurtosis on the entire dataset
            skewness_entire_dataset = self.__calculation.calculateSkewness(data)
            kurtosis_entire_dataset = self.__calculation.calculateKurtosis(data)
            # Identify outliers and OOR data points on the entire dataset
            outliers_info_entire_dataset, outliers_set_entire_dataset = self.__calculation.identifyOutliers(data)
            oor_entire_dataset = self.__calculation.identifyOORDataPoints(data)

            # Compute skewness and kurtosis per data sources
            skewness_sources = self.__calculation.calculateSkewness_perDataset(data)
            kurtosis_sources = self.__calculation.calculateKurtosis_perDataset(data)
            # Identify outliers and OOR data points per data source
            outliers_info_sources, outliers_set_sources = self.__calculation.identifyOutliers_perDataset(data)
            oor_sources = self.__calculation.identifyOORDataPoints_perDataset(data)

        else: 
            skewness_entire_dataset = kurtosis_entire_dataset = outliers_info_entire_dataset = outliers_set_entire_dataset = oor_entire_dataset = None
            skewness_sources = kurtosis_sources = outliers_info_sources = outliers_set_sources = oor_sources = None

        # Compare statistically the endpoint distribution across data sources
        endpoint_distribution_results = self.__calculation.compareEndpointDistribution_acrossDatasets(data)

        # Compute the distance matirx and get the within- and betwen-source distances per data source
        distance_matrix = self.__calculation.calculateDistanceMatrix(data)
        feature_similarity_results = self.__calculation.calculateFeatureSimilarity(distance_matrix, data)
        feature_similarity_pairwise_results = self.__calculation.calculateFeatureSimilarity_Pairwise(distance_matrix, data)

        # Create the main directory and the endpoint subdirectory
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
        if not os.path.exists(os.path.join(self.directory, self.endpoint_name)):
            os.makedirs(os.path.join(self.directory, self.endpoint_name))

        # Generate the final DataFrame and save it into a TSV file
        logging.info(f"Creating comparative report on {self.endpoint_name} aggregated dataset")

        summary_df = self.__outputfile.getSummaryDataFrame(n_mols_total, n_mols_sources, statistics_entire_dataset, statistics_sources, prop_ref_mols_total, prop_ref_mols_sources, 
                                                           skewness_entire_dataset, skewness_sources, kurtosis_entire_dataset, kurtosis_sources, outliers_info_entire_dataset, 
                                                           outliers_info_sources, oor_entire_dataset, oor_sources, endpoint_distribution_results, feature_similarity_results)
        self.__outputfile.writeToTSV(summary_df)

        if self.task == config["TASK_REGRESSION"]:
            # Visualize outliers
            self.__visualization.VisualizeOutliers(data, outliers_set_entire_dataset) 

        # Plot the intersection across data sources
        self.__visualization.plotDatasetsIntersection(data)

        # Plot the endpoint distribution
        self.__visualization.plotEndpointDistribution(data, outliers_set_entire_dataset)

        # Plot the endpoint distribution for each data source
        self.__visualization.plotEndpointDistributionComparison(data, outliers_set_sources)

        if self.reference_set is not None:
            # Plot the endpoint distribution for reference vs. non-reference molecules for each data source
            self.__visualization.plotEndpointDistributionComparison_RefVsNonRefMols(data, ref_inchikeys)

        # Inspect the discrepancies between data sources and Plot the comparative heatmaps
        discrepancies_df = self.__calculation.calculateInterDatasetDiscrepancies(data)
        self.__visualization.plotComparativeHeatmaps(discrepancies_df)

        if self.task == config["TASK_REGRESSION"]:
            # Perform Pairwise KS test and Plot the test results in a heatmap
            ks_results = self.__calculation.perfromPairwiseKSTest(data)
            self.__visualization.plotPairwiseKSTestHeatmap(ks_results)

        # Plot the similarity distribution
        self.__visualization.plotSimilarityDistribution(distance_matrix)
        self.__visualization.plotFeatureSimilarityHeatmap(feature_similarity_pairwise_results)

        # Run UMAP and Visualize the feature space
        projection = self.__calculation.runUMAP(data)
        self.x_range, self.y_range = self.__visualization.getAxisRanges(projection)
        self.__visualization.plotFeatureSpace_coloredbyEndpoint(projection, data)
        self.__visualization.plotFeatureSpace_coloredbyDataset(projection, data)
        self.__visualization.plotFeatureSpace_KDEplot(projection, data)
        self.__visualization.plotFeatureSpace_Hexbin(projection)    

        # Generate the final insight report
        self.__insightreport.generateInsightReport(data, skewness_sources, outliers_info_sources, oor_sources, feature_similarity_results, 
                                                         endpoint_distribution_results, discrepancies_df, prop_ref_mols_sources)
    
        logging.info(f"The final report and several plots have been saved in the {os.path.join(self.directory, self.endpoint_name)} directory")

        self.__outputfile.generateOutputSummary(summary_df)

    ###
    @property
    def columns(self):
        return self.__columns
    
    @columns.setter
    def columns(self, value):
        self.__columns = value
