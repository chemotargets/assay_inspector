<div align="center">
  <h2>
    Data consistency assessment facilitates transfer learning in ADME modeling
  </h2>
</div>

<div align="center">
  <img src="https://raw.githubusercontent.com/chemotargets/assay_inspector/main/AssayInspector.png", alt="AssayInspector", width=60%>
</div>

&nbsp;

Data heterogeneity and distributional misalignments pose critical challenges for machine learning models, often compromising predictive accuracy. These challenges are exemplified in preclinical safety modeling, a crucial step in early-stage drug discovery where limited data and experimental constraints exacerbate integration issues. Analyzing public ADME datasets, we uncovered significant misalignments between benchmark and gold-standard sources that degrade model performance. Our analyses further revealed that dataset discrepancies arise from differences in various factors, from experimental conditions in data collection to chemical space coverage. This highlights the importance of **rigorous data consistency assessment (DCA) prior to modeling**. To facilitate a systematic DCA across diverse datasets, we developed AssayInspector, a **model-agnostic package** that leverages *statistics*, *visualizations*, and *diagnostic summaries* to identify *outliers*, *batch effects*, and *discrepancies*. Beyond preclinical safety, DCA can play a crucial role in federated learning scenarios, enabling effective transfer learning across heterogeneous data sources and supporting reliable integration across diverse scientific domains.

**Keywords:** data reporting, molecular property, ADME, physicochemical, machine learning, data aggregation, predictive accuracy, benchmark

### Usage
**Download** from GitHub:  
```
wget https://github.com/chemotargets/assay_inspector/archive/refs/heads/master.zip
unzip master.zip
cd assay_inspector-master/
```

**Create Environment**:  
```
conda env create -f AssayInspector_env.yml
```

**Perform Analysis**

```
from AI_Main import AssayInspector

report = AssayInspector(
	data_path='path/to/dataset/file.tsv',
	endpoint_name='endpoint',
	task='regression',
	feature_type='ecfp4,
	reference_set='path/to/reference_set.tsv'
)

report.get_individual_reporting()
report.get_comparative_reporting()
```


<!--
#### Citation Note
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