# PhyloTwin-Preprocessor - GBIF Occurrence Data Processing Pipeline

A Nextflow pipeline for processing GBIF species occurrence data with spatial outlier detection and H3 grid binning. 

The pipeline offers two processing modes:
- "Atomic" mode: Processes each species independently
- "Batched" mode: Processes multiple species in batches (optimized for HPC environments)


