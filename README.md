# AbSplice analysis

Repository for the analyses done in the [publication](https://www.nature.com/articles/s41588-023-01373-3) "Aberrant splicing prediction across human tissues".

# Folder Structure
The project is setup as a snakemake pipeline. 

* The `workflow` folder contains the config files and the analysis scripts for each dataset.
    * `config` contains yaml files for the GTEx, mito and ALS datasets, as well as a yaml that defines the file structure used in all analysis. 
    * `scripts` contains all processing and analysis scripts. The `common` folder contains all scripts shared across datasets. `gtex_v8`, `mito` and `als` contain analysis scripts specific to the respective dataset. Each of the dataset folders contains a Snakefile that runs the analysis of the respective dataset. As mito and als depend on gtex_v8 results (e.g. AbSplice model, SpliceMaps from GTEx tissues), the gtex_v8 Snakefile needs to be executed first.

# Overview of the analysis pipeline
In `scripts/common` most of the scripts of the pipeline are included.
The main tasks of the pipeline are:
* `vcf_annotation` filters the provided vcf file for rare, high quality variants 
* `junction_annotation` generates SpliceMaps based on RNA-seq split read counts from FRASER
* `outlier_ground_truth` generates FRASER outliers on gene and junction level, filtered for significance and rare variants in the vicinity
* `splicing_pred` generates splicing predictions for SpliceAI, MMSplice, MTSplice, CADD-Splice, SQUIRLS and AbSplice
* `splicing_result` provides maximum predictions on a gene and variant level
* `benchmark` creates the benchmark dataframes from predictions and outliers. Where each gene, sample, tissue combination in the benchmark contains at least one rare variant on the gene

# Necessary files to provide by the user
In order for the pipeline to work users need to provide some files:

- outlier calls from DROP. Create a symlink to the DROP results in the folder `./workflow/data/resources/{dataset}/DROP/`, where dataset can be `gtex_v8`, `mito` and `als`.  
We used DROP v.1.1.2. For GTEx v8 the necessary config file to run DROP is provided in [here](https://github.com/gagneurlab/AbSplice_analysis/blob/master/workflow/data/resources/gtex_v8/DROP_config.yaml): `./workflow/data/resources/gtex_v8/DROP_config.yaml`, with the corresponding sample annotation for DROP in [here](https://github.com/gagneurlab/AbSplice_analysis/blob/master/workflow/data/resources/gtex_v8/DROP_sample_annotation.tsv): `./workflow/data/resources/gtex_v8/DROP_sample_annotation.tsv`.

- rocksdb database with precomputed SpliceAI scores. Create a symlink in `./workflow/data/resources/common/{genome}/SpliceAI/spliceai.db`, where genome can be `hg19` and `hg38`. These databases can be created from: https://github.com/gagneurlab/spliceai_rocksdb.
- rocksdb database with gnomAD minor allele frequencies. Create a symlink in `./workflow/data/resources/common/{genome}/gnomAD_maf_db/rocksdb/maf.db`, where genome can be `hg19` and `hg38`. These databases can be created from: https://github.com/gagneurlab/gnomad_rocksdb.
- vcf files from the dataset. Users need to normalize the vcf files and store them in `./workflow/data/resources/{dataset}/vcf_normalized/{vcf_id}.vcf.gz`, where dataset can be `gtex_v8`, `mito` and `als`. `vcf_id` can be anything (e.g. sampleID or split huge vcf by chromosome to parallelize computations).

# Run the pipeline
Follow these steps to run the pipeline:
- Install the absplice environment. Follow the steps indicated in the [absplice GitHub](https://github.com/gagneurlab/absplice/tree/master). You can use the docker image or conda environment.
- Activate the absplice environment.
- From the absplice environment, run the snakemake workflow for gtex_v8 (below just with 1 job; increase this depending on your available resources):
    ```
    cd workflow/scripts/gtex_v8
    python -m snakemake -j 1 --use-conda
    ```
- Similarly run the snakemake workflow for [ALS](https://github.com/gagneurlab/AbSplice_analysis/tree/master/workflow/scripts/als). 
- The mito dataset is not publicly available.

# Figure generation
All figures from the manuscript can be generated from the the provided scripts in [here](https://github.com/gagneurlab/AbSplice_analysis/tree/master/figures_R). The minimal dataset to produce the figures is available [here](https://zenodo.org/record/7628916). To run the scripts, first download and unzip the provided data folder to the root directory of this repository.








