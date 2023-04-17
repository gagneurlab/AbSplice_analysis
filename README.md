# AbSplice analysis

Repository for the analyses done in the paper "Aberrant splicing across human tissues"

# Folder Structure
The project is setup as a snakemake pipeline. 

* The `workflow` folder contains the config files and the analysis scripts for each dataset.
    * `config` contains yaml files for the GTEx, mito and ALS datasets, as well as a yaml that defines the file structure used in all analysis. 
    * `scripts` contains all processing and analysis scripts. The `common` folder contains all scripts shared acorss datasets. `gtex_v8`, `mito` and `als` contain analysis scripts specific to the respective dataset. Each of the dataset folders contains a Snakefile that runs the analysis of the respective dataset. As mito and als depend on gtex_v8 results (e.g. AbSplice model, SpliceMaps from GTEx tissues), the gtex_v8 Snakefile needs to be executed first.

# Overview of the analysis pipeline
In `scripts/common` most of the scripts of the pipeline are included.
The main tasks of the pipeline are:
* `vcf_annotation` filters the provided vcf file for rare, high quality variants 
* `junction_annotation` generates SpliceMaps based on RNA-seq split read counts from FRASER
* `outlier_ground_truth` generates FRASER outliers on gene and junction level, filtered for significance and rare variants in the vicinity
* `splicing_pred` generates splicing predictions for SpliceAI, MMSplice, MTSplice, CADD-Splice, SQUIRLS and AbSplice
* `splicing_result` provides maximum predictions on a gene and variant level
* `benchmark` creates the benchmark dataframes from predictions and outliers. Where each gene, sample, tissue combination in the benchmark contains at least one rare variant on the gene

# Figure generation
All figures from the manuscript can be generated from the the provided scripts in [here](https://github.com/gagneurlab/AbSplice_analysis/tree/master/figures_R). The minimal dataset to produce the figures is available [here](https://zenodo.org/record/7628916). To run the scripts, first download and unzip the provided data folder to the root directory of this repository.








