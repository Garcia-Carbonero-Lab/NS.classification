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
R -f scripts/preprocess/exprrssion.matrix.R

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
