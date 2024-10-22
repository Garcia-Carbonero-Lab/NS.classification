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

1. **Clinical differences**: `scripts/classification/analysis/clinical.R` (any environment with R)
2. **Immune score differences**: `scripts/classification/analysis/immune.R` (deg environment)
3. **Differentially expressed genes**: `scripts/classification/analysis/deg.R` (deg environment)
4. **Differentially methylated probes**: `scripts/classification/analysis/dep.R` (deg environment)
5. **Differentially inferred activated master regulators**: `scripts/classification/analysis/master.regulators.R` (deg environment)
6. **Gene set enrichment analysis**: `scripts/classification/analysis/gsea.R` (gsea environment)
7. **Survival analysis**: `scripts/classification/analysis/survival.NS.R` (survival environment)
8. **Predictive evaluation of Cox models with and without neuroendocrine subtypes**: `scripts/classification/analysis/survival.eval.R` (surveval environment) (Note is surveval and not survival environment)
9. **Prognostic impact of differentially expressed genes and methylation probes**: `scripts/classification/analysis/prognostic.features.R` (survival environment)
10. **Correlation between gene expression and CpG methylation**: `scripts/classification/analysis/corr.met.exp.R` (correlation.exp.met environment)

# Validation Analysis

To corroborate the robustness of the results obtained, the same analysis approach has been applied to an external cohort of patients with transcriptomic and methylation data available. This external cohort consists of samples from four different studies and three different primary tumors:

- **[Alcala et al.](https://pubmed.ncbi.nlm.nih.gov/31431620/)**: This cohort consists of 34 lung NETs. Transcriptomic data is available in the EGA under the code [EGAD00010001719](https://ega-archive.org/datasets/EGAD00010001719), and can also be downloaded along with clinical data from the Fernández-Cuesta lab GitHub repository [IARCbioinfo/DRMetrics](https://github.com/IARCbioinfo/DRMetrics). Methylation data is available in the EGA under the code [EGAD00010001720](https://ega-archive.org/datasets/EGAD00010001720). For this analysis, the folder `data/external_data/Alcala_LUNET/` should be present. The expression data file, `read_counts_all.txt`, must be moved to this folder, and the methylation data should be placed in a folder called `IDAT` inside the `data/external_data/Alcala_LUNET/` folder.

- **[Laddha et al.](https://pubmed.ncbi.nlm.nih.gov/31300474/)**: This cohort consists of 18 lung NETs. Transcriptomic data is available in GEO under the code [GSE118131](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118131) and in the same matrix as Alcala et al. data in the Fernández-Cuesta lab GitHub repository [IARCbioinfo/DRMetrics](https://github.com/IARCbioinfo/DRMetrics). Methylation data is available under the code [GSE118132](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118132). For this analysis, the folder `data/external_data/Laddha_LUNET/` should be present. The expression data file with genes as rows and samples as columns, `counts.tsv`, must be moved to this folder, and the methylation data should be placed in a folder called `IDAT` inside the `data/external_data/Laddha_LUNET/` folder.

- **[Lakis et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC7859232/)**: This cohort consists of 11 pancreatic NETs. Transcriptomic, methylomic, and clinical data were available in the ICGC portal until June 2024. To obtain this data, please contact us at ccarretero.imas12@h12o.es. All files should be moved to the `data/external_data/Lakis_PANET/` folder.

- **[Karpathakis et al.](https://pubmed.ncbi.nlm.nih.gov/26169971/)**: This cohort consists of 27 small intestine NETs. Methylation data is available under the code [GSE73832](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73832). However, transcriptomic and clinical data are available only upon request from the authors. The methylation matrix, `matrix.txt`, and the requested expression data should be moved to the `data/external_data/Karpathakis_SINET/` folder.

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
    * **Filtering clinical data of Lakis**: `scripts/validation/preprocess/Lakis.base.filter.R` (geo)

* **Karpathakis Folder**: 
    * **Preprocess Transcriptomic Data of Karpathakis**: `scripts/validation/preprocess/Karpathakis/Karpathakis.expression.R` (geo environment)
    * **Preprocess Methylation Data of Lakis**: `scripts/validation/preprocess/Karpathakis/Karpathakis.methylation.R` (geo environment)

Once each dataset was preprocessed we combine and filter all datasets runing:
 
1. `scripts/validation/preprocess/combine.expression.R` (preprocess environment)
2. `scripts/validation/preprocess/combine.expression.R` (preprocess environment)
3. `scripts/validation/preprocess/Base.R` (geo environment)

For finally preprocess steps we obtained:

1. **CIMP calculation**: `scripts/validation/preprocess/calculate.CIMP.R` (preprocess environment)
2. **Immune infiltration deconvolution**: `scripts/validation/preprocess/calculate.CIMP.R` (preprocess environment)

## MOFA creation

We develop two models. The model using the discovery samples and the validation samples called the combination model and the discovery model using only validation samples.


