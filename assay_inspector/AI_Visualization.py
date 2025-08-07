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
import os

from itertools import cycle

import plotly.graph_objects as go
import plotly.io as pio
import plotly.express as px
import plotly.figure_factory as ff
from plotly.subplots import make_subplots

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib_venn import venn2, venn3
from upsetplot import UpSet, from_contents

from .AI_Utils import logging

### Configs
import os
import json

config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'configs.json')

with open(config_path, "r") as file:
    config = json.load(file)

class Visualization():
    """
    Class to generate all visualization plots. 
    """

    def __init__(self, mainSelf):
        self.__mainSelf = mainSelf

    ###
    def getAxisRanges(self, proj, expand_factor=0.1):

        """
        Calculates consistent axis ranges for feature space coverage plots to ensure square shape.
        """

        # Find the largest range
        largest_range = max(np.nanmax(proj[:, 0]) - np.nanmin(proj[:, 0]), np.nanmax(proj[:, 1]) - np.nanmin(proj[:, 1]))        
        # Calculate midpoints, expand the largest range, and set both axes to have the same range
        x_mid = (np.nanmax(proj[:, 0]) + np.nanmin(proj[:, 0])) /2
        y_mid = (np.nanmax(proj[:, 1]) + np.nanmin(proj[:, 1])) /2
        expanded_range = largest_range * (1 + expand_factor)
        x_range = (x_mid - (expanded_range / 2), x_mid + (expanded_range / 2))
        y_range = (y_mid - (expanded_range / 2), y_mid + (expanded_range / 2))

        return x_range, y_range

    ###
    def VisualizeOutliers(self, data_df, outliers_set):

        """
        Generates a scatter plot that enables the visualization of outliers on the endpoint sorted
        distribution. 
        """

        # Sort data by value
        data_df = data_df.sort_values(by=config["NAMES"]["VALUE"])
        # Flag molecules as outliers or non-outliers
        colors = data_df['inchikey'].apply(lambda x: 'red' if x in outliers_set else 'black')

        # Create the figure
        fig = go.Figure(data=[go.Scatter(y=data_df[config["NAMES"]["VALUE"]], mode='markers', marker=dict(color=colors))])
        fig.update_layout(
            title={'text':f'[1] Visualization of {self.__mainSelf.endpoint_name} Outliers', 'x':0.5, 'xanchor':'center', 'font':{'size':24}},
            yaxis_title={'text':f'{self.__mainSelf.endpoint_name} value', 'font':{'size':16}},
            xaxis_title={'text':f'Sorted molecules', 'font':{'size':16}},
            height=500, width=500,
            margin=dict(l=50, r=50, t=50, b=50),
        )

        # Save the figure
        pio.write_image(fig, os.path.join(self.__mainSelf.directory, self.__mainSelf.endpoint_name, '1_outlier_visualization.svg'), format='svg')

    ###
    def plotEndpointDistribution(self, data_df, outliers_set):

        """
        Generate the plot(s) which shows the endpoint distribution for the entire dataset.
        """

        # Define y-axis range (only for regression endpoints)
        yaxis_range = [self.__mainSelf.lower_bound, self.__mainSelf.upper_bound]
        if (self.__mainSelf.lower_bound is not None) and (data_df[config["NAMES"]["VALUE"]].min() < self.__mainSelf.lower_bound):
            yaxis_range[0] = data_df[config["NAMES"]["VALUE"]].min()
        if (self.__mainSelf.upper_bound is not None) and (data_df[config["NAMES"]["VALUE"]].max() > self.__mainSelf.upper_bound):
            yaxis_range[1] = data_df[config["NAMES"]["VALUE"]].max()

        # Bar chart for classification endpoints
        if self.__mainSelf.task == config["TASK_CLASSIFICATION"]:
            value_counts = data_df[config["NAMES"]["VALUE"]].value_counts().sort_index(ascending=True)

            # Create the figure
            fig = go.Figure(data=[go.Bar(x=value_counts.index, y=value_counts.values,
                                        text=[f"{(value/data_df.shape[0])*100:.2f}" for value in value_counts.values], textposition='auto',
                                        marker_color=self.__mainSelf.hex_colors[:len(value_counts.index)])])
            fig.update_layout(
                title={'text':f'[2] Distribution of {self.__mainSelf.endpoint_name} Classes', 'x':0.5, 'xanchor':'center', 'font':{'size':24}},
                xaxis_title={'text':f'{self.__mainSelf.endpoint_name} classes', 'font':{'size':16}},
                yaxis_title={'text':'Number of molecules', 'font':{'size':16}},
                height=400, width=600,
                margin=dict(l=50, r=50, t=50, b=50),
                xaxis=dict(tickmode='array', tickvals=value_counts.index.tolist(), ticktext=value_counts.index.tolist())
            )

            # Save the figure
            pio.write_image(fig, os.path.join(self.__mainSelf.directory, self.__mainSelf.endpoint_name, '2_endpoint_distribution_visualization.svg'), format='svg')

        # Histogram and Violin plot for regression endpoints
        elif self.__mainSelf.task == config["TASK_REGRESSION"]:

            # Create the histogram
            fig = go.Figure(data=[go.Histogram(x=data_df[config["NAMES"]["VALUE"]], marker_color=self.__mainSelf.hex_colors[0])])
            fig.update_layout(
                title={'text':f'[2.1] Distribution of {self.__mainSelf.endpoint_name} Values', 'x':0.5, 'xanchor':'center', 'font':{'size':24}},
                xaxis_title={'text':f'{self.__mainSelf.endpoint_name}', 'font':{'size':16}},
                yaxis_title={'text':'Number of molecules', 'font':{'size':16}},
                height=400, width=600,
                margin=dict(l=50, r=50, t=50, b=50),
            )

            # Save the histogram
            pio.write_image(fig, os.path.join(self.__mainSelf.directory, self.__mainSelf.endpoint_name, '2_1_endpoint_distribution_visualization_histogram.svg'), format='svg')

            # Create the violin plot
            fig = go.Figure(data=[go.Violin(y=data_df[config["NAMES"]["VALUE"]], box_visible=True, line_color='black', x0=' ',
                                            meanline_visible=True, meanline={'color':'red'}, 
                                            spanmode='manual', span=[self.__mainSelf.lower_bound, self.__mainSelf.upper_bound])])
            fig.update_layout(
                title={'text':f'[2.2] Distribution of {self.__mainSelf.endpoint_name} Values', 'x':0.5, 'xanchor':'center', 'font':{'size':24}},
                yaxis_title={'text':f'{self.__mainSelf.endpoint_name}', 'font':{'size':16}},
                height=600, width=400,
                margin=dict(l=50, r=50, t=50, b=50),
                yaxis=dict(range=yaxis_range)
            )
            
            # Save the violin plot
            pio.write_image(fig, os.path.join(self.__mainSelf.directory, self.__mainSelf.endpoint_name, '2_2_endpoint_distribution_visualization_violinplot.svg'), format='svg')

        # Additional plot to visualize endpoint distribution when including or excluding
        # outliers (only for regression endpoints)
        if self.__mainSelf.task == config["TASK_REGRESSION"]:
            # Get endpoint distribution excluding outliers
            data_without_outliers = data_df.loc[~data_df[config["NAMES"]["INCHIKEY"]].isin(outliers_set)]

            # Create the figure
            with_outliers = go.Violin(y=data_df[config["NAMES"]["VALUE"]], name='Outliers included', legendgroup='Outliers included', 
                                    side='negative', line_color=self.__mainSelf.hex_colors[0], x0=' ', meanline_visible=True, showlegend=True,
                                    spanmode='manual', span=[self.__mainSelf.lower_bound, self.__mainSelf.upper_bound], 
                                    marker=dict(outliercolor=self.__mainSelf.hex_colors[0]))
            without_outliers = go.Violin(y=data_without_outliers[config["NAMES"]["VALUE"]], name='Outliers excluded', legendgroup='Outliers excluded', 
                                        side='positive', line_color=self.__mainSelf.hex_colors[1], x0=' ', meanline_visible=True, showlegend=True,
                                        spanmode='manual', span=[self.__mainSelf.lower_bound, self.__mainSelf.upper_bound],
                                        marker=dict(outliercolor=self.__mainSelf.hex_colors[1], symbol='x', opacity=0.7))
            fig = go.Figure(data=[with_outliers] + [without_outliers])
            fig.update_layout(
                title={'text':f'[2.3] Distribution of {self.__mainSelf.endpoint_name} Values', 'x':0.5, 'xanchor':'center', 'font':{'size':24}},
                yaxis_title={'text':f'{self.__mainSelf.endpoint_name} value', 'font':{'size':16}},
                height=600, width=500,
                margin=dict(l=50, r=50, t=50, b=50),
                showlegend=True,
                xaxis=dict(range=[-0.4,0.4]),
                yaxis=dict(range=yaxis_range)
            )
        
            # Save the figure
            pio.write_image(fig, os.path.join(self.__mainSelf.directory, self.__mainSelf.endpoint_name, '2_3_endpoint_distribution_visualization_outliers.svg'), format='svg')

    ###
    def plotEndpointDistributionComparison(self, data_df, outliers_dict):

        """
        Generate the plot(s) comparing the endpoint distribution across individual datasets.
        """

        # Define y-axis range (only for regression endpoints)
        yaxis_range = [self.__mainSelf.lower_bound, self.__mainSelf.upper_bound]
        if (self.__mainSelf.lower_bound is not None) and (data_df[config["NAMES"]["VALUE"]].min() < self.__mainSelf.lower_bound):
            yaxis_range[0] = data_df[config["NAMES"]["VALUE"]].min()
        if (self.__mainSelf.upper_bound is not None) and (data_df[config["NAMES"]["VALUE"]].max() > self.__mainSelf.upper_bound):
            yaxis_range[1] = data_df[config["NAMES"]["VALUE"]].max()

        # Stacked Bar charts for classification endpoints
        if self.__mainSelf.task == config["TASK_CLASSIFICATION"]:
            # Calculate class counts and percentages for each source 
            source_class_counts = {ref: data_df[config["NAMES"]["VALUE"]].loc[data_df[config["NAMES"]["REF"]] == ref].value_counts(ascending=True).sort_index() for ref in self.__mainSelf.sources_list}
            source_class_percentages = {ref: (counts / counts.sum()) * 100 for ref, counts in source_class_counts.items()}

            # Define figure data
            bars_nmols = []
            bars_percentage = []
            bars_counts = []

            for i, ref in enumerate(self.__mainSelf.sources_list):
                # Get source data
                n_mols = data_df.loc[data_df[config["NAMES"]["REF"]] == ref].shape[0]
                counts = source_class_counts[ref].sort_index(ascending=True)
                percentages = source_class_percentages[ref].sort_index(ascending=True)
                # Number-of-Molecules Bar Chart
                bars_nmols.append(go.Bar(name=ref, x=[ref], y=[n_mols], marker_color='grey', showlegend=False, opacity=0.7,
                                        text=f"{n_mols}", textposition='outside')) 
                # Class-Percentages Bar Chart
                for j in range(len(percentages)):
                    bars_percentage.append(go.Bar(name=str(percentages.index[j]), x=[ref], y=[percentages.get(j, 0)], 
                                                marker_color=self.__mainSelf.hex_colors[j], opacity=0.7, showlegend=False))    
                # Class-Counts Bar Chart
                for k in range(len(counts)):
                    bars_counts.append(go.Bar(name=str(counts.index[k]), x=[ref], y=[counts.get(k, 0)], 
                                            marker_color=self.__mainSelf.hex_colors[k], opacity=0.7, legendgroup=str(counts.index[j]), showlegend=(i==0),
                                            text=f"{counts.get(k, 0)}", textposition='auto'))          

            # Create the combined figure
            fig = make_subplots(rows=3, cols=1, shared_xaxes=True, row_heights=[0.2, 0.2, 0.6],  vertical_spacing=0.05)

            # Add the number-of-molecules bars to the first subplot
            for bar in bars_nmols:
                fig.add_trace(bar, row=1, col=1)
            # Add the class-percentages bars to the second subplot
            for bar in bars_percentage:
                fig.add_trace(bar, row=2, col=1)
            # Add the class-counts bars to the third subplot
            for bar in bars_counts:
                fig.add_trace(bar, row=3, col=1)

            fig.update_layout(
                barmode='stack',
                title={'text': f'[3] Distribution of {self.__mainSelf.endpoint_name} Classes Across Data Sources', 'x': 0.5, 'xanchor': 'center', 'font': {'size': 24}},
                xaxis3_title={'text': 'Data Sources', 'font': {'size': 16}},  
                yaxis=dict(title={'text': 'Nº of Mols', 'font': {'size': 10}}),  
                yaxis2=dict(title={'text': 'Proportion (%)', 'font': {'size': 10}}, range=[0, 100], minor=dict(nticks=4, ticklen=5, tickcolor="black", showgrid=True)),
                yaxis3=dict(title={'text': 'Number of Molecules', 'font': {'size': 16}}),
                legend_title={'text':f'{self.__mainSelf.endpoint_name} classes', 'font':{'size':16}},
                height=600, width=800,
                margin=dict(l=50,r=50, t=50, b=50),
            )
            
            # Save the figure
            pio.write_image(fig, os.path.join(self.__mainSelf.directory, self.__mainSelf.endpoint_name, '3_endpoint_distribution_comparative_visualization.svg'), format='svg')

        # Overlaid Histogram and Multiple Violin plots for regression endpoints
        elif self.__mainSelf.task == config["TASK_REGRESSION"]:

            # Create the histogram
            fig = go.Figure()
            for i, ref in enumerate(self.__mainSelf.sources_list):
                source_data = data_df.loc[data_df[config["NAMES"]["REF"]] == ref]
                fig.add_trace(go.Histogram(x=source_data[config["NAMES"]["VALUE"]], name=ref, opacity=0.5, marker_color=self.__mainSelf.hex_colors[i % len(self.__mainSelf.hex_colors)]))
            fig.update_layout(
                barmode='overlay',
                title={'text':f'[3.1] Distribution of {self.__mainSelf.endpoint_name} Values Across Data Sources', 'x':0.5, 'xanchor':'center', 'font':{'size':24}},
                xaxis_title={'text':f'{self.__mainSelf.endpoint_name}', 'font':{'size':16}},
                yaxis_title={'text':'Number of molecules', 'font':{'size':16}},
                legend_title={'text':'Data sources', 'font':{'size':16}},
                height=400, width=800,
                margin=dict(l=50, r=50, t=50, b=50),
            )

            # Save the histogram
            pio.write_image(fig, os.path.join(self.__mainSelf.directory, self.__mainSelf.endpoint_name, '3_1_endpoint_distribution_comparative_visualization_histogram.svg'), format='svg')

            # Get the set of unique molecules from the entire dataset
            unique_mols = set(data_df[config["NAMES"]["INCHIKEY"]].tolist())
            
            # Define figure data
            bars_nmols = []
            bars_uniquemols = []
            violins = []

            for ref in self.__mainSelf.sources_list:
                # Get source data
                source_data = data_df.loc[data_df[config["NAMES"]["REF"]] == ref]
                n_mols = source_data.shape[0]
                prop_unique_mols = (len(unique_mols.intersection(set(source_data[config["NAMES"]["INCHIKEY"]].tolist()))) / len(unique_mols)) * 100
                # Number-of-Molecules Bar Chart
                bars_nmols.append(go.Bar(name=ref, x=[ref], y=[n_mols], marker_color=self.__mainSelf.hex_colors[0], showlegend=False, opacity=0.7,
                                        text=f"{n_mols}", textposition='outside')) 
                # Percentage-of-Unique-Molecules Bar Chart
                bars_uniquemols.append(go.Bar(name=ref, x=[ref], y=[prop_unique_mols], marker_color=self.__mainSelf.hex_colors[1], showlegend=False, opacity=0.7,
                                            text=f"{prop_unique_mols}", texttemplate='%{text:.2f}', textposition='outside'))
                # Violin Plot
                violins.append(go.Violin(y=source_data[config["NAMES"]["VALUE"]], name=ref, box_visible=True, line_color='black',
                                        meanline_visible=True, meanline={'color':'red'}, showlegend=False,
                                        spanmode='manual', span=[self.__mainSelf.lower_bound, self.__mainSelf.upper_bound])) 

            # Create the combined figure
            fig = make_subplots(rows=3, cols=1, shared_xaxes=True, row_heights=[0.15, 0.15, 0.7],  vertical_spacing=0.05)

            # Add the number-of-molecules bars to the first subplot
            for bar in bars_nmols:
                fig.add_trace(bar, row=1, col=1)
            # Add the percentage-of-unique-molecules bars to the second subplot
            for bar in bars_uniquemols:
                fig.add_trace(bar, row=2, col=1)
            # Add the violin plots to the third subplot
            for violin in violins:
                fig.add_trace(violin, row=3, col=1)    

            fig.update_layout(
                title={'text':f'[3.2] Distribution of {self.__mainSelf.endpoint_name} Values Across Data Sources', 'x':0.5, 'xanchor':'center', 'font':{'size':24}},
                xaxis3_title={'text': 'Data Sources', 'font': {'size': 16}},  
                yaxis=dict(title={'text': 'Nº of Mols', 'font': {'size': 10}}),  
                yaxis2=dict(title={'text': '% of Unique Mols', 'font': {'size': 10}}, range=[0, 100], minor=dict(nticks=4, ticklen=5, tickcolor="black", showgrid=True)), 
                yaxis3=dict(title={'text': f'{self.__mainSelf.endpoint_name} value', 'font': {'size': 16}}, range=yaxis_range), 
                height=700, width=1000,
                margin=dict(l=50, r=50, t=50, b=50),
                showlegend=False
            )

            # Save the violin plots
            pio.write_image(fig, os.path.join(self.__mainSelf.directory, self.__mainSelf.endpoint_name, '3_2_endpoint_distribution_comparative_visualization_violinplots.svg'), format='svg')

        # Additional plot to visualize endpoint distribution across datasets when including or 
        # excluding outliers (only for regression endpoints)
        if self.__mainSelf.task == config["TASK_REGRESSION"]:
            
            # Create the figure
            violins = []

            # Add dummy traces for legend entries
            violins.append(go.Violin(y=[None], name='Outliers included', legendgroup='Outliers included', line_color=self.__mainSelf.hex_colors[0], showlegend=True))
            violins.append(go.Violin(y=[None], name='Outliers excluded', legendgroup='Outliers excluded', line_color=self.__mainSelf.hex_colors[1], showlegend=True))

            for ref in self.__mainSelf.sources_list:
                # Select source data including and excluding outliers
                source_data = data_df.loc[data_df[config["NAMES"]["REF"]] == ref]
                source_data_no_outliers = source_data.loc[~source_data[config["NAMES"]["INCHIKEY"]].isin(outliers_dict[ref])]
                # Append Violin Plot traces for both endpoint distributions
                violins.append(go.Violin(y=source_data[config["NAMES"]["VALUE"]], name=ref, legendgroup='Outliers included', 
                                        side='negative', line_color=self.__mainSelf.hex_colors[0], meanline_visible=True, showlegend=False,
                                        spanmode='manual', span=[self.__mainSelf.lower_bound, self.__mainSelf.upper_bound],
                                        marker=dict(outliercolor=self.__mainSelf.hex_colors[0])))
                violins.append(go.Violin(y=source_data_no_outliers[config["NAMES"]["VALUE"]], name=ref,  legendgroup='Outliers excluded', 
                                        side='positive', line_color=self.__mainSelf.hex_colors[1], meanline_visible=True, showlegend=False,
                                        spanmode='manual', span=[self.__mainSelf.lower_bound, self.__mainSelf.upper_bound],
                                        marker=dict(outliercolor=self.__mainSelf.hex_colors[1], symbol='x', opacity=0.7)))

            fig = go.Figure(data=violins)
            fig.update_layout(
                title={'text':f'[3.3] Distribution of {self.__mainSelf.endpoint_name} Values Across Data Sources', 'x':0.5, 'xanchor':'center', 'font':{'size':24}},
                yaxis_title={'text':f'{self.__mainSelf.endpoint_name} value', 'font':{'size':16}},
                height=500, width=1100,
                margin=dict(l=50, r=50, t=50, b=50),
                showlegend=True,
                xaxis=dict(range=[1.6, len(self.__mainSelf.sources_list)+1.4]),
                yaxis=dict(range=yaxis_range)
            )

            # Save the figure
            pio.write_image(fig, os.path.join(self.__mainSelf.directory, self.__mainSelf.endpoint_name, '3_3_endpoint_distribution_comparative_visualization_outliers.svg'), format='svg')

    ###
    def plotEndpointDistributionComparison_RefVsNonRefMols(self, data_df, ref_set):

        """
        Generate the plot comparing the endpoint distribution of reference vs. non-reference
        molecules across individual datasets.

        This plot can only be created when a reference set of molecules is defined.
        """

        # Define y-axis range (only for regression endpoints)
        yaxis_range = [self.__mainSelf.lower_bound, self.__mainSelf.upper_bound]
        if (self.__mainSelf.lower_bound is not None) and (data_df[config["NAMES"]["VALUE"]].min() < self.__mainSelf.lower_bound):
            yaxis_range[0] = data_df[config["NAMES"]["VALUE"]].min()
        if (self.__mainSelf.upper_bound is not None) and (data_df[config["NAMES"]["VALUE"]].max() > self.__mainSelf.upper_bound):
            yaxis_range[1] = data_df[config["NAMES"]["VALUE"]].max()

        # Grouped and Stacked Bar charts for classification endpoints
        if self.__mainSelf.task == config["TASK_CLASSIFICATION"]:

            # Calculate class counts and percentages for each source 
            endpoint_classes = sorted(data_df[config["NAMES"]["VALUE"]].unique().tolist())
            source_class_counts_ref = {ref: data_df[config["NAMES"]["VALUE"]].loc[(data_df[config["NAMES"]["REF"]] == ref) & (data_df[config["NAMES"]["INCHIKEY"]].isin(ref_set))].value_counts().reindex(endpoint_classes, fill_value=0).sort_index(ascending=True) for ref in self.__mainSelf.sources_list}
            source_class_counts_nonref = {ref: data_df[config["NAMES"]["VALUE"]].loc[(data_df[config["NAMES"]["REF"]] == ref) & (~data_df[config["NAMES"]["INCHIKEY"]].isin(ref_set))].value_counts().reindex(endpoint_classes, fill_value=0).sort_index(ascending=True) for ref in self.__mainSelf.sources_list}
            source_class_percentages_ref = {ref: (counts / counts.sum()) * 100 for ref, counts in source_class_counts_ref.items()}
            source_class_percentages_nonref = {ref: (counts / counts.sum()) * 100 for ref, counts in source_class_counts_nonref.items()}

            # Define figure data
            bars_nmols = []
            bars_percentages_refmols = []
            bars_counts_refmols = []

            for i, ref in enumerate(self.__mainSelf.sources_list):
                # Get source data
                source_data = data_df.loc[data_df[config["NAMES"]["REF"]] == ref]

                n_mols = source_data.shape[0]
                prop_ref_mols = (len(ref_set.intersection(set(source_data[config["NAMES"]["INCHIKEY"]].tolist()))) / len(set(source_data[config["NAMES"]["INCHIKEY"]].tolist()))) * 100
                prop_nonref_mols = 100 - prop_ref_mols

                percentages_ref_mols = source_class_percentages_ref[ref].sort_index(ascending=True)
                percentages_nonref_mols = source_class_percentages_nonref[ref].sort_index(ascending=True)

                counts_ref_mols = source_class_counts_ref[ref].sort_index(ascending=True)
                counts_nonref_mols = source_class_counts_nonref[ref].sort_index(ascending=True)

                # Number-of-Molecules Bar Chart
                bars_nmols.append(go.Bar(name=f'{ref} (Total molecules)', x=[ref], y=[n_mols], marker_color='grey', showlegend=False, opacity=0.7,
                                        text=f"{n_mols}", textposition='outside')) 
                # Class-Percentages Bar Chart
                for j in range(len(percentages_ref_mols)):
                    if j == 0:
                        bars_percentages_refmols.append(go.Bar(name=f'{str(percentages_ref_mols.index[j])} (Reference molecules)', x=[ref], y=[percentages_ref_mols.get(j, 0)], 
                                                            marker_color=self.__mainSelf.hex_colors[j], showlegend=(i==0), offsetgroup='ref', opacity=0.7)) 
                        bars_percentages_refmols.append(go.Bar(name=f'{str(percentages_nonref_mols.index[j])} (Non-Reference molecules)', x=[ref], y=[percentages_nonref_mols.get(j, 0)], 
                                                            marker_color=self.__mainSelf.hex_colors[j+2], showlegend=(i==0), offsetgroup='nonref', opacity=0.7)) 
                    else:
                        bars_percentages_refmols.append(go.Bar(name=f'{str(percentages_ref_mols.index[j])} (Reference molecules)', x=[ref], y=[percentages_ref_mols.get(j, 0)], 
                                                            marker_color=self.__mainSelf.hex_colors[j], showlegend=(i==0), base=[percentages_ref_mols.get(j-1, 0)], offsetgroup='ref', opacity=0.7)) 
                        bars_percentages_refmols.append(go.Bar(name=f'{str(percentages_nonref_mols.index[j])} (Non-Reference molecules)', x=[ref], y=[percentages_nonref_mols.get(j, 0)], 
                                                            marker_color=self.__mainSelf.hex_colors[j+2], showlegend=(i==0), base=[percentages_nonref_mols.get(j-1, 0)], offsetgroup='nonref', opacity=0.7)) 

                # Class-Counts Bar Chart
                for k in range(len(counts_ref_mols)):
                    if k == 0:
                        bars_counts_refmols.append(go.Bar(name=f'{str(counts_ref_mols.index[k])} (Reference molecules)', x=[ref], y=[counts_ref_mols.get(k, 0)], 
                                                        marker_color=self.__mainSelf.hex_colors[k], showlegend=False, offsetgroup='ref', opacity=0.7)) 
                        bars_counts_refmols.append(go.Bar(name=f'{str(counts_nonref_mols.index[k])} (Non-Reference molecules)', x=[ref], y=[counts_nonref_mols.get(k, 0)], 
                                                        marker_color=self.__mainSelf.hex_colors[k+2], showlegend=False, offsetgroup='nonref', opacity=0.7)) 
                    else:
                        bars_counts_refmols.append(go.Bar(name=f'{str(counts_ref_mols.index[k])} (Reference molecules)', x=[ref], y=[counts_ref_mols.get(k, 0)], 
                                                        marker_color=self.__mainSelf.hex_colors[k], showlegend=False, base=[counts_ref_mols.get(k-1, 0)], offsetgroup='ref', opacity=0.7,
                                                        text=f"{prop_ref_mols}", texttemplate='%{text:.2f}%', textposition='outside')) 
                        bars_counts_refmols.append(go.Bar(name=f'{str(counts_nonref_mols.index[k])} (Non-Reference molecules)', x=[ref], y=[counts_nonref_mols.get(k, 0)], 
                                                        marker_color=self.__mainSelf.hex_colors[k+2], showlegend=False, base=[counts_nonref_mols.get(k-1, 0)], offsetgroup='nonref', opacity=0.7,
                                                        text=f"{prop_nonref_mols}", texttemplate='%{text:.2f}%', textposition='outside',))  
        
            # Create the combined figure
            fig = make_subplots(rows=3, cols=1, shared_xaxes=False, row_heights=[0.2, 0.2, 0.6],  vertical_spacing=0.05)

            # Add the number-of-molecules bars to the first subplot
            for bar in bars_nmols:
                fig.add_trace(bar, row=1, col=1)
            # Add the class-percentages bars to the second subplot
            for bar in bars_percentages_refmols:
                fig.add_trace(bar, row=2, col=1)
            # Add the class-counts bars to the third subplot
            for bar in bars_counts_refmols:
                fig.add_trace(bar, row=3, col=1)

            fig.update_layout(
                barmode='group',
                title={'text': f'[ADD.] Comparison of {self.__mainSelf.endpoint_name} Classes for Reference vs. Non-Reference Molecules Across Data Sources', 'x': 0.5, 'xanchor': 'center', 'font': {'size': 20}},
                xaxis1=dict(showticklabels=False),
                xaxis2=dict(showticklabels=False),
                xaxis3_title={'text': 'Data Sources', 'font': {'size': 16}},  
                yaxis1=dict(title={'text': 'Nº of Mols', 'font': {'size': 10}}),  
                yaxis2=dict(title={'text': 'Proportion (%)', 'font': {'size': 10}}, range=[0, 100], minor=dict(nticks=4, ticklen=5, tickcolor="black", showgrid=True)),
                yaxis3=dict(title={'text': 'Number of Molecules', 'font': {'size': 16}}),
                legend_title={'text':f'{self.__mainSelf.endpoint_name} classes', 'font':{'size':16}},
                legend=dict(x=1, y=0.70, xanchor='left', yanchor='middle'),
                height=700, width=1000,
                margin=dict(l=50, r=50, t=50, b=50),
            )

        # Mutliple Split Violin plots for regression endpoints
        elif self.__mainSelf.task == config["TASK_REGRESSION"]:
            
            # Define figure data
            bars_nmols = []
            bars_refmols = []
            violins = []

            # Add dummy traces for legend entries
            violins.append(go.Violin(y=[None], name='Reference molecules', legendgroup='Reference molecules', line_color=self.__mainSelf.hex_colors[1], showlegend=True))
            violins.append(go.Violin(y=[None], name='Non-reference molecules', legendgroup='Non-reference molecules', line_color=self.__mainSelf.hex_colors[2], showlegend=True))

            for ref in self.__mainSelf.sources_list:
                # Get source data
                source_data = data_df.loc[data_df[config["NAMES"]["REF"]] == ref]
                # Calculate statistics
                n_mols = source_data.shape[0]
                prop_ref_mols = (len(ref_set.intersection(set(source_data[config["NAMES"]["INCHIKEY"]].tolist()))) / len(set(source_data[config["NAMES"]["INCHIKEY"]].tolist()))) * 100
                prop_nonref_mols = 100 - prop_ref_mols
                # Number-of-Molecules Bar Chart
                bars_nmols.append(go.Bar(name=ref, x=[ref], y=[n_mols], marker_color=self.__mainSelf.hex_colors[0], showlegend=False, opacity=0.7,
                                        text=f"{n_mols}", textposition='outside')) 
                # Percentage-of-Reference-Molecules Bar Chart
                bars_refmols.append(go.Bar(name='Reference molecules', x=[ref], y=[prop_ref_mols], marker_color=self.__mainSelf.hex_colors[1], showlegend=False, opacity=0.7))
                bars_refmols.append(go.Bar(name='Non-reference molecules', x=[ref], y=[prop_nonref_mols], marker_color=self.__mainSelf.hex_colors[2], showlegend=False, opacity=0.7))
                
                # Select source data for reference and non-reference molecules
                source_data_ref_mols = source_data.loc[source_data[config["NAMES"]["INCHIKEY"]].isin(ref_set)]
                source_data_nonref_mols = source_data.loc[~source_data[config["NAMES"]["INCHIKEY"]].isin(ref_set)]
                # Split Violin PlotS
                if len(source_data_ref_mols) > 0:
                    violins.append(go.Violin(y=source_data_ref_mols[config["NAMES"]["VALUE"]], name=ref, legendgroup='Reference molecules', 
                                            side='negative', line_color=self.__mainSelf.hex_colors[1], meanline_visible=True, showlegend=False,
                                            spanmode='manual', span=[self.__mainSelf.lower_bound, self.__mainSelf.upper_bound],
                                            marker=dict(outliercolor=self.__mainSelf.hex_colors[1])))
                if len(source_data_nonref_mols) > 0:
                    violins.append(go.Violin(y=source_data_nonref_mols[config["NAMES"]["VALUE"]], name=ref,  legendgroup='Non-reference molecules', 
                                            side='positive', line_color=self.__mainSelf.hex_colors[2], meanline_visible=True, showlegend=False,
                                            spanmode='manual', span=[self.__mainSelf.lower_bound, self.__mainSelf.upper_bound],
                                            marker=dict(outliercolor=self.__mainSelf.hex_colors[2], symbol='x', opacity=0.7)))
                
            # Create the combined figure
            fig = make_subplots(rows=3, cols=1, shared_xaxes=True, row_heights=[0.15, 0.15, 0.7],  vertical_spacing=0.05)

            # Add the number-of-molecules bars to the first subplot
            for bar in bars_nmols:
                fig.add_trace(bar, row=1, col=1)
            # Add the percentage-of-reference-molecules bars to the second subplot
            for bar in bars_refmols:
                fig.add_trace(bar, row=2, col=1)
            # Add the split violin plots to the third subplot
            for violin in violins:
                fig.add_trace(violin, row=3, col=1)
            
            fig.update_layout(
                barmode='stack',
                title={'text':f'[EXTRA] Comparison of {self.__mainSelf.endpoint_name} Values for Reference vs. Non-Reference Molecules Across Data Sources', 'x':0.5, 'xanchor':'center', 'font':{'size':20}},
                xaxis3_title={'text': 'Data Sources', 'font': {'size': 16}},  
                yaxis=dict(title={'text': 'Nº of Mols', 'font': {'size': 10}}),  
                yaxis2=dict(title={'text': 'Proportion (%)', 'font': {'size': 10}}, range=[0, 100], minor=dict(nticks=4, ticklen=5, tickcolor="black", showgrid=True)), 
                yaxis3=dict(title={'text': f'{self.__mainSelf.endpoint_name} value', 'font': {'size': 16}}, range=yaxis_range), 
                legend=dict(x=1, y=0.75, xanchor='left', yanchor='middle'),
                height=700, width=1000,
                margin=dict(l=50, r=50, t=50, b=50),
            )

        # Save the figure
        pio.write_image(fig, os.path.join(self.__mainSelf.directory, self.__mainSelf.endpoint_name, 'EXTRA_endpoint_distribution_comparative_visualization_refmols.svg'), format='svg')

    ###
    def plotFeatureSpace_coloredbyEndpoint(self, proj, data_df):

        """
        Generates the plot showing the projection of the molecule feature space colored by 
        endpoint value.
        """

        # UMAP projection colored by endpoint value 
        if self.__mainSelf.task == config["TASK_CLASSIFICATION"]:
            data_df[config["NAMES"]["VALUE"]] = data_df[config["NAMES"]["VALUE"]].astype('category')
            fig = px.scatter(proj, x=0, y=1,
                            color=data_df[config["NAMES"]["VALUE"]], 
                            color_discrete_sequence=self.__mainSelf.hex_colors, category_orders={'value': data_df[config["NAMES"]["VALUE"]].cat.categories.tolist()},
                            labels={'color':f'{self.__mainSelf.endpoint_name} class'}, opacity=0.5)
            fig.update_layout(
                title={'text':f'[4] Feature Space Coverage', 'x':0.5, 'xanchor':'center', 'font':{'size':24}},
                xaxis_title={'text':f'UMAP1', 'font':{'size':16}},
                yaxis_title={'text':f'UMAP2', 'font':{'size':16}},
                legend={'font': {'size': 16}},
                height=800, width=800,
                margin=dict(l=50, r=50, t=50, b=50),
                xaxis=dict(range=self.__mainSelf.x_range), yaxis=dict(range=self.__mainSelf.y_range)            
            )

        elif self.__mainSelf.task == config["TASK_REGRESSION"]:
            fig = px.scatter(proj, x=0, y=1,
                            color=data_df[config["NAMES"]["VALUE"]], 
                            color_continuous_scale=px.colors.sequential.Viridis, color_continuous_midpoint=data_df[config["NAMES"]["VALUE"]].mean(),
                            labels={'color':f'{self.__mainSelf.endpoint_name}'}, opacity=0.5)
            fig.update_layout(
                title={'text':f'[4] Feature Space Coverage', 'x':0.5, 'xanchor':'center', 'font':{'size':24}},
                xaxis_title={'text':'UMAP1', 'font':{'size':16}},
                yaxis_title={'text':'UMAP2', 'font':{'size':16}},
                legend={'font': {'size': 16}},
                height=800, width=800,
                margin=dict(l=50, r=50, t=50, b=50),
                xaxis=dict(range=self.__mainSelf.x_range), yaxis=dict(range=self.__mainSelf.y_range)            
            )

        # Save the figure
        fig.write_image(os.path.join(self.__mainSelf.directory, self.__mainSelf.endpoint_name, '4_UMAP_projection_coloredby_endpoint.svg'), format='svg')


    ###
    def plotFeatureSpace_coloredbyDataset(self, proj, data_df):
        """
        Generates the plot showing the projection of the molecular feature space colored by 
        data source.
        """

        # UMAP projection colored by data source
        fig = px.scatter(proj, x=0, y=1,
                        color=data_df[config["NAMES"]["REF"]], 
                        color_discrete_sequence=self.__mainSelf.hex_colors, category_orders={'color':self.__mainSelf.sources_list},
                        labels={'color':'Data source'}, opacity=0.5)
        fig.update_layout(
            title={'text':'[5] Feature Space Coverage', 'x':0.5, 'xanchor':'center', 'font':{'size':24}},
            xaxis_title={'text':'UMAP1', 'font':{'size':16}},
            yaxis_title={'text':'UMAP2', 'font':{'size':16}},
            legend={'font': {'size': 16}},
            height=800, width=1000,
            margin=dict(l=50, r=50, t=50, b=50),
            xaxis=dict(range=self.__mainSelf.x_range), yaxis=dict(range=self.__mainSelf.y_range)            
        )

        # Save the figure
        fig.write_image(os.path.join(self.__mainSelf.directory, self.__mainSelf.endpoint_name, '5_UMAP_projection_coloredby_dataset.svg'), format='svg')

    ###
    def plotFeatureSpace_KDEplot(self, proj, data_df):

        """
        Generates the kernel density estimate (KDE) plot showing the projection of the molecular
        feature space across individual datasets. 
        """

        # Define the palette based on the number of sources
        if len(self.__mainSelf.hex_colors) < len(self.__mainSelf.sources_list):
            palette = [color for color, _ in zip(cycle(self.__mainSelf.hex_colors), range(len(self.__mainSelf.sources_list)))]
        else:
            palette = self.__mainSelf.hex_colors[:len(self.__mainSelf.sources_list)]

        # KDE plot
        plt.figure(figsize=(9, 9))
        kde = sns.kdeplot(data=pd.DataFrame(proj, columns=['0', '1']), x='0', y='1', 
                        hue=data_df[config["NAMES"]["REF"]], hue_order=self.__mainSelf.sources_list, palette=palette)

        plt.title('[6] Feature Space Coverage (KDE plot)', fontsize=20)
        plt.xlabel('UMAP1', fontsize=14)
        plt.ylabel('UMAP2', fontsize=14)
        plt.xlim(self.__mainSelf.x_range)
        plt.ylim(self.__mainSelf.y_range)
        plt.grid(True)
        plt.tight_layout()

        # Save the figure
        plt.savefig(os.path.join(self.__mainSelf.directory, self.__mainSelf.endpoint_name, '6_UMAP_projection_KDEplot.svg'), format='svg')
        plt.close()

    ###
    def plotFeatureSpace_Hexbin(self, proj):

        """
        Generates the plot showing the density of the molecular feature space in a Hexagonal 
        Binning plot.
        """

        # Density plot
        fig = go.Figure(go.Histogram2d(x=proj[:, 0], y=proj[:, 1], histnorm='density', 
                                    autobinx=False, xbins=dict(start=self.__mainSelf.x_range[0], end=self.__mainSelf.x_range[1], size=(self.__mainSelf.x_range[1] - self.__mainSelf.x_range[0])/100),
                                    autobiny=False, ybins=dict(start=self.__mainSelf.y_range[0], end=self.__mainSelf.y_range[1], size=(self.__mainSelf.y_range[1] - self.__mainSelf.y_range[0])/100),
                                    colorbar=dict(title='Point count'), colorscale='Electric'))
        fig.update_layout(
            title={'text':'[7] Feature Space Density (Hexbin)', 'x':0.5, 'xanchor':'center', 'font':{'size':24}},
            xaxis_title={'text':'UMAP1', 'font':{'size':16}},
            yaxis_title={'text':'UMAP2', 'font': {'size': 16}},
            legend={'font': {'size': 16}},
            height=800, width=800,
            margin=dict(l=50, r=50, t=50, b=50),
            xaxis=dict(range=self.__mainSelf.x_range), yaxis=dict(range=self.__mainSelf.y_range)
        )

        # Save the figure
        pio.write_image(fig, os.path.join(self.__mainSelf.directory, self.__mainSelf.endpoint_name, '7_UMAP_projection_Hexbin.svg'), format='svg')

    ###
    def plotDatasetsIntersection(self, data_df):

        """
        Generates a plot showing the intersection of molecules across individual datasets.

        The type of plot is selected based on the number of sources:
        - A Venn diagram is used when there are two or three datasets.
        - An UpSet plot is used when there are more than three datasets.
        """

        # Get the molecule set for each dataset
        mols_per_source = {}
        for ref in self.__mainSelf.sources_list:
            mols_per_source[ref] = data_df[config["NAMES"]["INCHIKEY"]].loc[data_df[config["NAMES"]["REF"]] == ref].unique().tolist()

        # 2 DATASETS
        if len(mols_per_source) == 2: 
            # Retrieve the InChIKey set for each dataset
            set1 = set(mols_per_source[self.__mainSelf.sources_list[0]])
            set2 = set(mols_per_source[self.__mainSelf.sources_list[1]])
            # Plot the Venn diagram
            plt.figure(figsize=(8, 8))
            venn2([set1, set2], (self.__mainSelf.sources_list[0], self.__mainSelf.sources_list[1]), set_colors=(self.__mainSelf.hex_colors[0], self.__mainSelf.hex_colors[1]))
            plt.title('[8] Venn Diagram', fontsize=24)
            plt.tight_layout()

        # 3 DATASETS
        elif len(mols_per_source) == 3: 
            # Retrieve the InChIKey set for each dataset
            set1 = set(mols_per_source[self.__mainSelf.sources_list[0]])
            set2 = set(mols_per_source[self.__mainSelf.sources_list[1]])
            set3 = set(mols_per_source[self.__mainSelf.sources_list[2]])
            # Plot the Venn diagram
            plt.figure(figsize=(8, 8))
            venn3([set1, set2, set3], (self.__mainSelf.sources_list[0], self.__mainSelf.sources_list[1], self.__mainSelf.sources_list[2]), set_colors=(self.__mainSelf.hex_colors[0], self.__mainSelf.hex_colors[1], self.__mainSelf.hex_colors[2]))
            plt.title('[8] Venn Diagram', fontsize=24)
            plt.tight_layout()
        
        # MORE THAN 3 DATASETS
        elif len(mols_per_source) > 3:
            # Pass the dictionary to the from_contents() function
            upset_data = from_contents(mols_per_source)
            # Plot the UpSet plot
            matplotlib.rcParams["font.size"] = 10

            if len(self.__mainSelf.sources_list) == 4: # 4 data sources
                upset = UpSet(upset_data, subset_size='count', sort_by='cardinality', show_counts="{:,}", facecolor=self.__mainSelf.hex_colors[0], element_size=40) 
            else:
                try: 
                    upset = UpSet(upset_data, min_subset_size=data_df.shape[0]*0.01, subset_size='count', # only show subsets with at least 1% of the total data points
                                sort_by='cardinality', show_counts="{:,}", facecolor=self.__mainSelf.hex_colors[0], element_size=40) 
                except:
                    logging.warning(oriMessage="The UpSet plot has not been generated because no subset contains at least 1% of the total data points.")
                    return
            
            plot = upset.plot(fig=plt.figure(figsize=(8, 8)))
            plot["intersections"].set_ylabel("Subset size")
            plot["totals"].set_xlabel("Dataset size")
            plt.suptitle('[8] Upset Plot', fontsize=24)
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)    
        
        # Save the image
        plt.savefig(os.path.join(self.__mainSelf.directory, self.__mainSelf.endpoint_name, '8_intersection_across_datasets.svg'), format='svg')
        plt.close()

    ###
    def plotComparativeHeatmaps(self, info_df):

        """
        Generates the plot which integrates two complementary pairwise heatmaps showing the discrepencies
        between data sources based on the proportion of common molecules and the differences in endpoint
        values for these shared compounds.
        """

        # Endpoint type specifications
        if self.__mainSelf.task == config["TASK_CLASSIFICATION"]:
            diff_col = 'different_values_count_proportion' 
            colorbar_title = 'Proportion of different<br>values'
            text_template = '%{text}'
        elif self.__mainSelf.task == config["TASK_REGRESSION"]:
            diff_col = 'mean_abs_diff' 
            colorbar_title = 'Mean Absolute Difference<br>(standardized)'
            text_template = '%{text:.3f}'

        # Create the pivot tables for the first heatmap
        num_common_mols_df = info_df.pivot(index='source_1', columns='source_2', values='num_common_mols')
        num_common_mols_df = num_common_mols_df.loc[self.__mainSelf.sources_list[::-1], self.__mainSelf.sources_list] # reorder heatmap matrix
        percent_common_mols_df = info_df.pivot(index='source_1', columns='source_2', values='percent_common_mols')
        percent_common_mols_df = percent_common_mols_df.loc[self.__mainSelf.sources_list[::-1], self.__mainSelf.sources_list] # reorder heatmap matrix

        # Create the pivot table for the second heatmap 
        difference_df = info_df.pivot(index='source_1', columns='source_2', values=diff_col)
        difference_df = difference_df.loc[self.__mainSelf.sources_list, self.__mainSelf.sources_list] # reorder heatmap matrix
        difference_values_df = difference_df.copy()

        if self.__mainSelf.task == config["TASK_CLASSIFICATION"]:
            difference_values_df = info_df.pivot(index='source_1', columns='source_2', values='different_values_count')

        # Create the figure
        fig = make_subplots(rows=1, cols=2, subplot_titles=('Common Molecules Counts', f'{self.__mainSelf.endpoint_name} Value Difference'), shared_yaxes=True)

        # Common-Molecules Heatmap 
        fig.add_trace(go.Heatmap(z=percent_common_mols_df.values, x=percent_common_mols_df.columns, y=percent_common_mols_df.index,
                                text=num_common_mols_df.replace(np.nan, "").values, texttemplate='%{text}',
                                colorbar=dict(title='Percentage<br>(from left source)', len=1, x=0.45), colorscale='Reds', zmin=0, zmax=100), row=1, col=1) 
        # Endpoint-Differences Heatmap 
        fig.add_trace(go.Heatmap(z=difference_df.values, x=difference_df.columns, y=difference_df.index,
                                text=difference_values_df.replace(np.nan, "").values, texttemplate=text_template, 
                                colorbar=dict(title=colorbar_title, len=1, x=1), colorscale='Reds', zmin=0, zmax=1), row=1, col=2) 
        
        fig.update_layout(
            title={'text':'[9] Pairwise Between-Source Discrepancies', 'x':0.5, 'xanchor':'center', 'font':{'size':24}},
            height=600, width=1500,
            margin=dict(l=50, r=50, t=50, b=50)
        )

        # Save the figure
        pio.write_image(fig, os.path.join(self.__mainSelf.directory, self.__mainSelf.endpoint_name, '9_comparative_heatmaps.svg'), format='svg')

    ###
    def plotSimilarityDistribution(self, distance_matrix):

        """
        Generates the plot showing the distribution of intermolecular distances in the entire
        dataset.
        """
        
        # Feature type specifications
        if self.__mainSelf.feature_type in ['rdkit']:
            distance_metric = 'Standardized Euclidean distance'
            xrange = [distance_matrix.min(), distance_matrix.max()]
        elif self.__mainSelf.feature_type in ['ecfp4']:
            distance_metric = 'Tanimoto distance'
            xrange = [0, 1]
        elif self.__mainSelf.feature_type in ['custom']:
            distance_metric = f'{self.__mainSelf.distance_metric.capitalize()} distance'
            xrange = [distance_matrix.min(), distance_matrix.max()]

        # Randomly sample 10,000 molecular distances
        if (distance_matrix.shape[0] * (distance_matrix.shape[0] - 1) // 2) > 10000:
            random_subset = np.random.choice(distance_matrix[np.triu_indices_from(distance_matrix, k=1)], 10000, replace=False)
        else:
            random_subset = distance_matrix[np.triu_indices_from(distance_matrix, k=1)]

        # Create the figure
        fig = go.Figure(data=[go.Histogram(x=random_subset, marker_color=self.__mainSelf.hex_colors[0])])
        fig.update_layout(
            title={'text':f'[10] Distribution of {distance_metric} values', 'x':0.5, 'xanchor':'center', 'font':{'size':24}},
            xaxis_title={'text':f'{distance_metric}', 'font':{'size':16}},
            yaxis_title={'text':''},
            height=400, width=600,
            margin=dict(l=50, r=50, t=50, b=50),
            yaxis=dict(showticklabels=False),
            xaxis=dict(range=xrange)
        )

        # Save the figure
        pio.write_image(fig, os.path.join(self.__mainSelf.directory, self.__mainSelf.endpoint_name, '10_intermolecular_distances_visualization.svg'), format='svg')

    ###
    def plotFeatureSimilarityHeatmap(self, similarity_results):

        """
        Generates the plot showing the feature profile similarities across individual datasets
        by displaying the average between-source distances. 
        """

        # Feature type specifications
        if self.__mainSelf.feature_type in ['rdkit','custom']:
            zmax = similarity_results['mean_inter_similarity'].max()
        elif self.__mainSelf.feature_type in ['ecfp4']:
            zmax = 1

        # Create the pivot table for the heatmap 
        similarity_table = similarity_results.pivot(index='ref1', columns='ref2', values='mean_inter_similarity') 
        similarity_table = similarity_table.loc[self.__mainSelf.sources_list, self.__mainSelf.sources_list] # reorder heatmap matrix

        # Create the figure
        fig = go.Figure(data=[go.Heatmap(z=similarity_table.values, x=similarity_table.columns, y=similarity_table.index,
                                        text=similarity_table.replace(np.nan, "").values, texttemplate='%{text:.2f}',
                                        colorbar=dict(title=f'Average Between-Distance<br>(from left source)', len=1, x=1), colorscale='amp', 
                                        zmin=0, zmax=zmax)])
        fig.update_layout(
            title={'text':f'[11] Pairwise Between-Source Distances', 'x':0.5, 'xanchor':'center', 'font':{'size':24}},
            height=600, width=800,
            margin=dict(l=50, r=50, t=50, b=50),
        )
        fig.update_annotations(font_size=24)
        fig.update_yaxes(autorange="reversed")

        # Save the figure
        pio.write_image(fig, os.path.join(self.__mainSelf.directory, self.__mainSelf.endpoint_name, '11_interdataset_distances_heatmap.svg'), format='svg')

    def plotPairwiseKSTestHeatmap(self, ks_results):

        """
        Generates the heatmap displaying the pairwise results of the two-sample KS test. Mean 
        differences are shown only for statistically significant pairwise comparisons (p-value < 
        alpha: 0.05). The overall plot includes a side dendrogram illustrating the clustering 
        of the datasets.

        This plot can only be generated for regression endpoints.
        """

        # Create the pivot tables for the heatmap 
        ks_df_meandiff = ks_results.pivot(index='group1', columns='group2', values='meandiff') 
        ks_df_reject = ks_results.pivot(index='group1', columns='group2', values='reject')

        # Mask non-significant pairs by setting their mean difference value to NaN
        ks_df_meandiff_masked = ks_df_meandiff.where(ks_df_reject, np.nan)
        ks_df_meandiff_masked = ks_df_meandiff_masked.combine_first(ks_df_meandiff_masked.T)
        ks_df_meandiff_masked = ks_df_meandiff_masked.replace(np.nan, 0)
        np.fill_diagonal(ks_df_meandiff_masked.values, np.nan)

        # Create the upper dendrogram and get the reordered indices
        dendro_side = ff.create_dendrogram(ks_df_meandiff_masked.replace(np.nan, 0).values, orientation='top', labels=ks_df_meandiff_masked.index)
        dendro_order = dendro_side['layout']['xaxis']['ticktext']
        for trace in dendro_side['data']:
            trace['xaxis'] = 'x2'
            trace['yaxis'] = 'y2' 

        # Sort the heatmap data based on the dendrogram order
        ks_df_meandiff_masked = ks_df_meandiff_masked.loc[dendro_order, dendro_order]
        ks_df_meandiff_masked = ks_df_meandiff_masked.iloc[::-1, :]

        # Define a custom colorscale
        colorscale = px.colors.sequential.ice
        colorscale[0] = '#D3D3D3'

        # Create the heatmap
        heatmap = go.Heatmap(z=ks_df_meandiff_masked.values, x=ks_df_meandiff_masked.columns, y=ks_df_meandiff_masked.index, 
                            text=ks_df_meandiff_masked.replace([np.nan, 0], "").values, texttemplate='%{text:.3f}',
                            colorbar=dict(title='Mean difference', len=1, x=1), colorscale=colorscale,
                            zmin=0, zmax=ks_results['meandiff'].max(), showscale=True)

        # Combine the dendrogram and heatmap
        fig = go.Figure(data=list(dendro_side['data']) + [heatmap])

        # Update layout
        fig.update_layout(
            title={'text':f'[12] Two-sample KS Test Results', 'x':0.5, 'xanchor':'center', 'font':{'size':24}},
            height=700, width=850,
            margin=dict(l=50, r=50, t=50, b=50),
            showlegend=False, 
            xaxis=dict(domain=[0, 1]),
            xaxis2=dict(domain=[0,1], showticklabels=False, autorange='reversed'),
            yaxis=dict(domain=[0, 0.775]),
            yaxis2=dict(domain=[0.775, 1], showticklabels=False, autorange='reversed')
        )
        fig.update_annotations(font_size=24)

        # Save the figure
        pio.write_image(fig, os.path.join(self.__mainSelf.directory, self.__mainSelf.endpoint_name, '12_kstest_heatmap.svg'), format='svg')
