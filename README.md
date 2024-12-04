# PhyloTwin-Preprocessor - GBIF Occurrence Data Processing Pipeline

A Nextflow pipeline for processing GBIF species occurrence data with spatial outlier detection and H3 grid binning. 

The pipeline offers two processing modes:
- "Atomic" mode: Processes each species independently
- "Batched" mode: Processes multiple species in batches (optimized for HPC environments)


# Dependencies

## Primary dependencies

These are the tools that the user must install and configure on their system to run the pipeline:

- [Nextflow](https://www.nextflow.io/) >= 24.10 (requires [Java](https://www.java.com/en/) >= 17 & <= 23)
- [Singularity](https://sylabs.io/)/[Apptainer](https://apptainer.org/) or [Docker](https://www.docker.com/)

## Secondary dependencies

These are the containerized tools or packages required by the pipeline, which will be automatically handled within the containers:

- [DuckDB](https://duckdb.org/)
- [ELKI](https://elki-project.github.io/)
- [R](https://www.r-project.org/)
- [Python](https://www.python.org/)
- [GNU Parallel](https://www.gnu.org/software/parallel/)
- [aria2](https://aria2.github.io/)
- [jq](https://github.com/jqlang/jq)

R packages:
- [data.table](https://rdatatable.gitlab.io/data.table/)
- [ape](https://github.com/emmanuelparadis/ape)


