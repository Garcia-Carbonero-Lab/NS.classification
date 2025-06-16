# NS.Classification
This repository contains the code and some of the data used to generate and analyze neuroendocrine molecular subtypes (NS) from diverse tumor origins. The main analysis paper is available [here](URL to paper).

## Data Preparation

The first step is to download data that is not included in the repository and move it to the corresponding folders.  
To develop the molecular classification, including custom data and downstream analysis, you need to download the following files:

1. **Transcriptomic data**: Download raw transcriptomic data (.CEL files) from GEO using the [GeoXXX](URL) accession. Once downloaded, create the folder `data/transcriptome/raw/files/` and move all files there.
2. **Methylomic data**: Download raw methylomic data (.IDAT files) from GEO using the [GeoXXX](URL) accession. After downloading, create the folder `data/methylome/raw/IDAT/` and move the files into it.
3. **Infinium Methylation Annotation File**: Download the methylation annotation file from the Illumina website [here](https://webdata.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip) and place it in the `data/methylome/annotation/` folder.

## Config Preparation

Fill out the configuration template available in the `config/` folder. In the template, specify the path to the working directory (`wkdir`) where the results will be saved and call it config.tsv 

## Environment Setup

Each script requires a specific environment, indicated in the script as a comment at the top. The composition of each environment is detailed in a `.yml` file located in the `envs/` folder. To manage environments, you should [install Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

To create an environment, open a terminal and run the following command, replacing the file path as needed:


```
  conda env create -f envs/preprocess.yml

```
Before running a script, activate the required environment using `conda activate`, as shown below:
```
  conda activate preprocess
```


## Neuroendocrine Subtype Generation and Analysis Using Scripts

Once the appropriate environment is activated, scripts can be run to reproduce the analysis using the R `-f` flag. For example:

```
R -f scripts/preprocess/expression.matrix.R

```

### Data Preprocessing

The first step is to preprocess the data and prepare files for the Multi-Omics Factor Analysis (MOFA) and subsequent analysis. Follow the order below, as some scripts generate outputs that are required by later scripts:

1. **Create normalized expression matrix**: `scripts/preprocess/expression.matrix.R` (arrays_m environment)
2. **Create normalized methylation matrix**: `scripts/preprocess/methylation.matrix.R` (methylation_m environment)
3. **Apply quality correction to both omics**: `scripts/preprocess/quality.correction.R` (preprocess environment)
4. **Feature selection for MOFA model creation**: `scripts/preprocess/feature.selection.R` (preprocess environment)
5. **CIMP calculation**: `scripts/preprocess/CIMP.calculation.R` (preprocess environment)
6. **Immune infiltration deconvolution**: Run over the corrected transcriptomic data to obtain immune infiltration scores, using `preprocess/immune.infiltration.matrix.R` (immune environment)
7. **Gene set filtering**: `scripts/genesets/filter.genesets.R` (geneset environment)
8. **Gene set annotation**: `scripts/genesets/anotate.genesets.R` (geneset environment)
9. **Master regulator network generation**: `scripts/master_regulators/aracne.R` (master_regulator environment)
10. **Master regulator activation inference**: `scripts/master_regulators/viper.R` (master_regulator environment)

### MOFA Analysis
After preprocessing, the next step is to integrate expression and methylation data using MOFA to generate a MOFA model. This model comprises eight factors derived separately from each available primary tumor type (lung, stomach, pancreas, small intestine, and colorectal). Each latent factor represents a driving source of variation across data modalities. These factors should then be examined to determine their relationship with clinical and biological features and the variability they capture.

1. **MOFA model creation**: `scripts/MOFA/MOFA.creation.R` (MOFA environment)
2. **MOFA model downstream analysis**: `scripts/MOFA/MOFA.analysis.R` (MOFA environment)

### Clustering Process

We applied an ensemble classification strategy using the MOFA model, excluding factors 3 and 7, as they represent unwanted variability related to expression quality and hepatocyte contamination, respectively. The ensemble was built using the consensus of eight classification algorithms (hierarchical clustering, K-means, DIANA, PAM, mixed Gaussian models, BLOCK, C-means, and HDBSCAN) via the CSPA method.

1. **Ensemble Classification Generation**: `clustering/classification.R` (clustering.ensemble environment)

### Downstream Clustering Analysis

Once the molecular subtypes were identified, we analyzed the molecular and clinical characteristics of each subtype.

1. **Clinical differences**: `scripts/classification.analysis/clinical.R` (any environment with R)
2. **Immune score differences**: `scripts/classification.analysis/immune.R` (deg environment)
3. **Differentially expressed genes**: `scripts/classification.analysis/deg.R` (deg environment)
4. **Differentially methylated probes**: `scripts/classification.analysis/dep.R` (deg environment)
5. **Differentially inferred activated master regulators**: `scripts/classification.analysis/master.regulators.R` (deg environment)
6. **Gene set enrichment analysis**: `scripts/classification.analysis/gsea.R` (gsea environment)
7. **Survival analysis**: `scripts/classification.analysis/survival.NS.R` (survival environment)
8. **Predictive evaluation of Cox models with and without neuroendocrine subtypes**: `scripts/classification.analysis/survival.eval.R` (surveval environment) (Note is surveval and not survival environment)
9. **Prognostic impact of differentially expressed genes and methylation probes**: `scripts/classification.analysis/prognostic.features.R` (survival environment)
10. **Correlation between gene expression and CpG methylation**: `scripts/classification.analysis/corr.met.exp.R` (correlation.exp.met environment)

# Validation Analysis

To corroborate the robustness of the results obtained, the same analysis approach has been applied to an external cohort of patients with transcriptomic and methylation data available. This external cohort consists of samples from four different studies and three different primary tumors:

- **[Alcala et al.](https://pubmed.ncbi.nlm.nih.gov/31431620/)**: This cohort consists of 34 lung NETs. Transcriptomic data is available in the EGA under the code [EGAD00010001719](https://ega-archive.org/datasets/EGAD00010001719), and can also be downloaded along with clinical data from the Fernández-Cuesta lab GitHub repository [IARCbioinfo/DRMetrics](https://github.com/IARCbioinfo/DRMetrics). Methylation data is available in the EGA under the code [EGAD00010001720](https://ega-archive.org/datasets/EGAD00010001720). For this analysis, the folder `data/external_data/Alcala_LUNET/` should be present. The expression data file, `read_counts_all.txt`, must be moved to this folder, and the methylation data should be placed in a folder called `IDAT` inside the `data/external_data/Alcala_LUNET/` folder.

- **[Laddha et al.](https://pubmed.ncbi.nlm.nih.gov/31300474/)**: This cohort consists of 18 lung NETs. Transcriptomic data is available in GEO under the code [GSE118131](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118131) and in the same matrix as Alcala et al. data in the Fernández-Cuesta lab GitHub repository [IARCbioinfo/DRMetrics](https://github.com/IARCbioinfo/DRMetrics). Methylation data is available under the code [GSE118132](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118132). For this analysis, the folder `data/external_data/Laddha_LUNET/` should be present. The expression data file with genes as rows and samples as columns, `counts.tsv`, must be moved to this folder, and the methylation data should be placed in a folder called `IDAT` inside the `data/external_data/Laddha_LUNET/` folder.

- **[Lakis et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC7859232/)**: This cohort consists of 11 pancreatic NETs. Transcriptomic, methylomic, and clinical data were available in the ICGC portal until June 2024. To obtain this data, please contact us at ccarretero.imas12@h12o.es. All files should be moved to the `data/external_data/Lakis_PANET/` folder.

- **[Patte et al.](https://pubmed.ncbi.nlm.nih.gov/40038310/)**: This cohort consists of 117 small intestine NETs. Methylation data is available under the code [GSE229203](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229203), while transcriptomic and clinical data can be accessed through EGA under dataset [EGAD50000000906] (https://ega-archive.org/datasets/EGAD50000000906), In this analysis, we use the gene expression matrix as provided by the authors. To proceed, ensure the folder `data/external_data/SiNETs_Patte/` exists. The expression matrix  `expression_data_genesymbol.rds`, containing genes as rows and samples as columns, should be placed directly in this folder. The raw methylation data should be placed in a subfolder named `IDAT` within the same directory.

## Preprocess Public Data

First, apply the preprocessing steps to the expression and methylation data of each dataset separately, using the scripts in the `scripts/validation/preprocess/` folder:

* **Alcala Folder**:
    * **Preprocess Transcriptomic Data of Alcala**: `scripts/validation/preprocess/Alcala/Alcala.expression.R` (geo environment)
    * **Preprocess Methylation Data of Alcala**: `scripts/validation/preprocess/Alcala/Alcala.methylation.R` (met.arrays environment)

* **Laddha Folder**:
    * **Preprocess Transcriptomic Data of Laddha**: `scripts/validation/preprocess/Laddha/Laddha.expression.R` (geo environment)
    * **Preprocess Methylation Data of Laddha**: `scripts/validation/preprocess/Laddha/Laddha.methylation.R` (met.arrays environment)

* **Lakis Folder**: 
    
    * **Generation clinical data of Lakis**: `scripts/validation/preprocess/Lakis.base.generation.R` (any R installed)
    * **Preprocess Transcriptomic Data of Lakis**: `scripts/validation/preprocess/Lakis/Lakis.expression.R` (geo environment)
    * **Preprocess Methylation Data of Lakis**: `scripts/validation/preprocess/Lakis/Lakis.methylation.R` (geo environment)
    * **Filtering clinical data of Lakis**: `scripts/validation/preprocess/Lakis.base.filter.R` (geo environmet)

* **Patte Folder**: 
    * **Preprocess Transcriptomic Data of Patte**: `scripts/validation/preprocess/Patte/Patte.expression.R` (geo environment)
    * **Preprocess Methylation Data of Patte**: `scripts/validation/preprocess/Patte/Patte.methylation.R` (met.arrays environment)

Once each dataset was preprocessed we combine and filter all datasets runing:
 
1. `scripts/validation/preprocess/combine.expression.R` (preprocess environment)
2. `scripts/validation/preprocess/combine.expression.R` (preprocess environment)
3. `scripts/validation/preprocess/Base.R` (geo environment)

For finally preprocess steps we obtained:

1. **CIMP calculation**: `scripts/validation/preprocess/calculate.CIMP.R` (preprocess environment)
2. **Immune infiltration deconvolution**: `scripts/validation/preprocess/calculate.CIMP.R` (preprocess environment)

## MOFA creation


We developed two MOFA models:

1. **Discovery Model** – Trained exclusively on discovery samples.  
2. **Combination Model** – Trained on both discovery and validation samples.

The relevant analysis script is located at:  
`scripts/validation/MOFA/MOFA.comparations.R`  
(MOFA environment)

### Step 1: Factor Comparison

We begin by comparing the **Discovery MOFA model** with the **Combination MOFA model** to identify which latent factors from the discovery dataset correspond to those in the combination model.  
Script:  
`scripts/validation/MOFA/MOFA.comparations.R`
(MOFA environment)

### Step 2: Sample Classification into Clusters

Next, we classify samples from both the **Combination** and **Validation** MOFA models into three clusters.  
This classification is based on the latent factors corresponding to original **Factor 1** and **Factor 2**.

Scripts:
- `scripts/validation/cluster.process/classification.combination.R` (clustering.ensemble environment)
- `scripts/validation/cluster.process/classification.validation.R`
(clustering.ensemble environment)

### Step 3: Subtype Assignment (NSs)

Finally, we assign samples to **NSs (novel subtypes)** based on the following process:

1. Each validation sample is assigned to a subtype based on the majority subtype of the discovery samples that cluster with it in the **Combination MOFA model**.

2. The **NSs in the Validation MOFA model** are then annotated based on the majority subtype labels transferred from the Combination model.

Script:  
`scripts/validation/cluster.process/clustering.comparation.R` (any with R installed environment)

## SNs molecular and clinical associations validation ##

We analyze the associations between NSs and clinical and molecular data using the following scripts:

1. **Clinical differences**: `scripts/validation/cluster.analysis/clinical.R` (any environment with R)
2. **Immune score differences**: `scripts/validation/cluster.analysis/immune.R` (deg environment)
3. **Differentially expressed genes**: `scripts/validation/cluster.analysis/deg.R` (deg environment)
4. **Differentially methylated probes**: `sscripts/validation/cluster.analysis/dmp.R` (deg environment)
5. **Gene set enrichment analysis**: `scripts/validation/cluster.analysis/gsea.R` (gsea environment)
6. **Survival analysis**: `scripts/validation/cluster.analysis/survival.NS.R` (survival environment)
7. **CNVs enrichtment**: `scripts/validation/cnv.R` (clinical environment)
