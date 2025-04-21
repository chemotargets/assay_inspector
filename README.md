<div align="center">
  <h1>
    Data consistency assessment facilitates transfer learning in ADME modeling
  </h1>
  <p><i>AssayInspector: A Python package for diagnostic assessment of data consistency in molecular datasets</i></p>

  ![Python](https://img.shields.io/badge/python-3.9+-blue)
  ![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
  
</div>

<div align="center">
	<div style="background-color: white; display: inline-block; padding: 10px; border-radius: 10px; width:60vw">
		<img src="https://github.com/chemotargets/assay_inspector/blob/master/AssayInspector.png?raw=true" alt="AssayInspector">
	</div>
</div>

&nbsp;

Data heterogeneity and distributional misalignments pose critical challenges for machine learning models, often compromising predictive accuracy. These challenges are exemplified in preclinical safety modeling, a crucial step in early-stage drug discovery where limited data and experimental constraints exacerbate integration issues. Analyzing public ADME datasets, we uncovered significant misalignments between benchmark and gold-standard sources that degrade model performance. Our analyses further revealed that dataset discrepancies arise from differences in various factors, from experimental conditions in data collection to chemical space coverage. This highlights the importance of **rigorous data consistency assessment (DCA) prior to modeling**. To facilitate a systematic DCA across diverse datasets, we developed **AssayInspector**, a **model-agnostic package** that leverages *statistics*, *visualizations*, and *diagnostic summaries* to identify *outliers*, *batch effects*, and *discrepancies*. Beyond preclinical safety, DCA can play a crucial role in federated learning scenarios, enabling effective transfer learning across heterogeneous data sources and supporting reliable integration across diverse scientific domains.

**Keywords:** data reporting, molecular property, ADME, physicochemical, machine learning, data aggregation, predictive accuracy, benchmark

## Installation

To install the latest release and use the package, run the following command to download it from GitHub:  
```bash 
git clone https://github.com/chemotargets/assay_inspector.git
cd assay_inspector
```

Then, create the environment to run the package as follows: 
```bash
conda env create -f AssayInspector_env.yml
```

## Getting Started

To run `AssayInspector`, you first need to prepare your input data. The file should be in `.tsv` or `.csv` format and include the following required columns:
* `smiles`: The SMILES string representation of each molecule in the dataset.
* `value`: The annotated value for each molecule — use a numerical value for regression tasks or a binary label (0 or 1) for classification tasks.
* `ref`: The reference source name from which each value-molecule annotation was obtained.
* `endpoint`: The name of the endpoint to analyze.

You can find two example input files for the [half-life] (test_scripts/data_testing/logHL_aggregated_dataset.tsv) and [clearance] (test_scripts/data_testing/logCL_aggregated_dataset.tsv) datasets.

## Usage

Once the input data file has been prepared, you can run `AssayInspector` in the following way:

```python
from AI_Main import AssayInspector

# Prepare AssayInspector report
report = AssayInspector(
	data_path='path/to/dataset/file.tsv',
	endpoint_name='endpoint',
	task='regression',
	feature_type='ecfp4',
	reference_set='path/to/reference_set.tsv' # optional
)

# Run AssayInspector report
report.get_individual_reporting()
report.get_comparative_reporting()
```

#### AssayInspector arguments

| Argument | Type | Description |
| --- | --- | --- |
| `data_path` | `str` | Path to the input dataset file (`.csv` or `.tsv` format). |
| `endpoint_name` | `str` | Name of the endpoint to analyze. |
| `task` | `str` | Type of task: either `'regression'` or `'classification'`. |
| `feature_type` | `str` | Type of features to use: one of `'ecfp4'`, `'rdkit'`, or `'custom'`. |
| `outliers_method` | `str` | *(Optional)* Method to detect outliers: `'zscore'` *(default)* or `'iqr'`. |
| `distance_metric` | `str` | *(Optional)* Distance metric for custom descriptors: `'euclidean'` *(default)*. |
| `descriptors_df` | `pd.DataFrame` | *(Optional)* DataFrame containing molecular descriptors for dataset molecules (required when `feature_type='custom'`). |
| `reference_set` | `str` | *(Optional)* Path to an additional dataset used for comparative analysis. |
| `lower_bound` | `int` or `float` | *(Optional)* Lower bound to define the endpoint applicability domain. |
| `upper_bound` | `int` or `float` | *(Optional)* Upper bound to define the endpoint applicability domain. |

The resulting output will be saved in a folder named `AssayInspector_YYYYMMDD`, which will contain:
- A tabular file that summarizes key descriptive parameters for each data source.
- A comprehensive set of visualization plots that facilitate the detection of inconsistencies across data sources.
- An insight report containing multiple alerts and recommendations to guide data cleaning and preprocessing.

## Examples

Below are a few sample outputs generated by `AssayInspector`.

| Endpoint | Outlier Visualization | Endpoint Distribution Comparative Visualization |
|-------------------------------|-------------------------------|----------------------------|
| Half-life | ![Outlier Visualization](examples/1_outlier_visualization_logHL.svg) | ![Endpoint Distribution Comparative Visualization](examples/3_2_endpoint_distribution_comparative_visualization_violinplots_logHL.svg) |
| Clearance | ![Outlier Visualization](examples/1_outlier_visualization_logCL.svg) | ![Endpoint Distribution Comparative Visualization](examples/3_2_endpoint_distribution_comparative_visualization_violinplots_logCL.svg) |

## License

`AssayInspector` is licensed under the MIT License. See the [LICENSE](./LICENSE.md) file.

<!--
## Cite us
Please cite [our paper](url) if you use *AssayInspector* in your own work:

```
@article {TAG,
         title = {Data consistency assessment facilitates transfer learning in ADME modeling},
         author = {Parrondo-Pizarro, Raquel and Menestrina, Luca and Garcia-Serna, Ricard and Fernández-Torras, Adrià and Mestres, Jordi},
         journal = {Journal},
         volume = {Vol},
         year = {Year},
         doi = {doi},
         URL = {url},
         publisher = {Publisher},
}
```
-->