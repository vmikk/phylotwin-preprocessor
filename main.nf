#!/usr/bin/env nextflow

// Usage:
// nextflow run main.nf \
//   -resume \
//   --occurrences "gbif_occurrences/" \
//   --outdir "gbif_occurrences_processed/" \
//   --basis_of_record "PRESERVED_SPECIMEN,MATERIAL_CITATION,MACHINE_OBSERVATION"
//   --t1 "./original_trees/" \
//   --t2 "./processed_trees/" \
//   --outlier_h3_resolution 6 \
//   --outlier_occurrence_threshold 100 \
//   --dbscan_minpts 3 \
//   --dbscan_eps 1000000 \
//   --h3_resolution 4 \
//   --userid "user"


// List of phylogenetic trees
// NB. taxon names (first element of tuple) should be unique
phylo_trees = [

    // Curated trees
    tuple("Plants",      "${params.t1}/Tietje_2023_Seed_plants_TACT.tre.gz",           "${params.t2}/Tietje_2023_Seed_plants_TACT.nex.gz"),
    tuple("Ferns",       "${params.t1}/FTOL_1-7-0_sanger_con_dated.tre.gz",            "${params.t2}/Ferns_FTOL_1-7-0.nex.gz"),
    tuple("Mushrooms",   "${params.t1}/Varga_2019_Mushrooms.tre.gz",                   "${params.t2}/Varga_2019_Mushrooms.nex.gz"),
    tuple("Mammals",     "${params.t1}/VertLife_Mammals.tre.gz",                       "${params.t2}/VertLife_Mammals.nex.gz"),
    tuple("Birds",       "${params.t1}/VertLife_Birds.tre.gz",                         "${params.t2}/VertLife_Birds.nex.gz"),
    tuple("Squamates",   "${params.t1}/VertLife_Squamates.tre.gz",                     "${params.t2}/VertLife_Squamates.nex.gz"),
    tuple("Amphibians",  "${params.t1}/VertLife_Amphibians.tre.gz",                    "${params.t2}/VertLife_Amphibians.nex.gz"),
    tuple("Fish",        "${params.t1}/FishTree_actinopt_12k_treePL.tre.gz",           "${params.t2}/FishTree.nex.gz"),
    tuple("Butterflies", "${params.t1}/Kawahara_2023_Butterflies_Newick.tre.gz",       "${params.t2}/Kawahara_2023_Butterflies.nex.gz"),
    tuple("Bees",        "${params.t1}/Henriquez-Piskulich_2024_BeetreeOfLife.tre.gz", "${params.t2}/Henriquez-Piskulich_2024_BeetreeOfLife.nex.gz"),

    // Synthetic trees from Open Tree of Life
    tuple("OTT_Aves",       "${params.t1}/OTL_Aves-1.3.nwk.gz",          "${params.t2}/OTL_Aves-1.3.nex.gz"),
    // NB. For the following synthetic trees, there are no initial trees, only processed trees
    tuple("OTT_Bryophyta",  "${params.t2}/OTT_v15.1__Bryophyta.nex.gz",  "${params.t2}/OTT_v15.1__Bryophyta.nex.gz"),
    tuple("OTT_Arthropoda", "${params.t2}/OTT_v15.1__Arthropoda.nex.gz", "${params.t2}/OTT_v15.1__Arthropoda.nex.gz"),
    tuple("OTT_Chordata",   "${params.t2}/OTT_v15.1__Chordata.nex.gz",   "${params.t2}/OTT_v15.1__Chordata.nex.gz"),
    tuple("OTT_Annelida",   "${params.t2}/OTT_v15.1__Annelida.nex.gz",   "${params.t2}/OTT_v15.1__Annelida.nex.gz"),
    tuple("OTT_Mollusca",   "${params.t2}/OTT_v15.1__Mollusca.nex.gz",   "${params.t2}/OTT_v15.1__Mollusca.nex.gz"),
    tuple("OTT_Nematoda",   "${params.t2}/OTT_v15.1__Nematoda.nex.gz",   "${params.t2}/OTT_v15.1__Nematoda.nex.gz"),
    tuple("OTT_Cnidaria",   "${params.t2}/OTT_v15.1__Cnidaria.nex.gz",   "${params.t2}/OTT_v15.1__Cnidaria.nex.gz"),
    tuple("OTT_Fungi",      "${params.t2}/OTT_v15.1__Fungi.nex.gz",      "${params.t2}/OTT_v15.1__Fungi.nex.gz"),

]

// GTDB_ar53_r220.tre.gz
// GTDB_bac120_r220.tre.gz

params.task_batching = false
params.batchsize = 40   // .buffer(size: params.batchsize, remainder: true)
params.cleanupwd = true


// Directory for publishing outputs
OUTDIR = params.userid ? params.outdir + "/" + params.userid : params.outdir
OUT_0_TRR = OUTDIR + "/Tree_stats"         // Number of GBIF records represented and not represented in phylogenetic trees
OUT_1_PRQ = OUTDIR + "/Filtered_parquet"   // Filtered and binned species occurrences
OUT_2_OUT = OUTDIR + "/Outliers"           // Pre-binned species occurrences with outlier scores

// Count number of occurrences per species
process count_occurrences {

    input:
      path occurrences

    output:
      path "Occurrences_counts.parquet", emit: occ_counts_parquet
      path "Occurrences_counts.csv.gz",  emit: occ_counts_csv

    script:
    def tempDirArg = task.ext.tempDir ? "-x ${task.ext.tempDir}" : ""
    def memoryArg  = task.memory  ? "-m ${task.memory.toMega()}.MB" : ""
    def basisOfRecordArg = params.basis_of_record ? "-b ${params.basis_of_record}" : ""
    """
    echo -e "Counting species occurrences\n"

    echo "Input directory: " ${occurrences}
    echo "CPUs: " ${task.cpus}
    if [ -n "${task.memory}"  ];     then echo "Memory: ${memoryArg}";          fi
    if [ -n "${task.ext.tempDir}" ]; then echo "Temp directory: ${tempDirArg}"; fi

    gbif_records_count.sh \
      -i ${occurrences} \
      -o "Occurrences_counts.parquet" \
      -t ${task.cpus} \
      ${memoryArg} ${tempDirArg} \
      ${basisOfRecordArg}

    if [ ${params.cleanupwd} == true ]; then
      echo "Cleaning up"
      rm -r ${occurrences}
    fi

    echo "..Done"
    """
}


// Define species set (given phylogenetic trees available)
process define_species_set {

    publishDir "${OUT_0_TRR}", pattern: '*_stats.txt', mode: "${params.publish_dir_mode}", overwrite: true 

    tag "${taxon}"

    input:
      tuple val(taxon), path(phylotree_initial, stageAs: "t1/*"), path(phylotree_processed, stageAs: "t2/*"), path(occ_counts_csv)

    output:
      path "${taxon}_Occurrences_large.txt",  emit: occ_large
      path "${taxon}_Occurrences_small.txt",  emit: occ_small
      path "${taxon}_stats.txt",              emit: stats

    script:
    """
    echo -e "Defining species set\n"
    echo "Occurrences counts: "   ${occ_counts_csv}
    echo "Phylogenetic tree: "    ${taxon}
    echo "Processed tree file: "  ${phylotree_processed}
    echo "Initial tree file: "    ${phylotree_initial}
    echo "Occurrence threshold: " ${params.outlier_occurrence_threshold}

    echo -e "\nCounting the number of GBIF records represented and not represented in a phylogenetic tree\n"

    gbif_records_in_tree.R \
      --gbif                 ${occ_counts_csv} \
      --tree                 ${phylotree_processed} \
      --sourcetree           ${phylotree_initial} \
      --occurrence_threshold ${params.outlier_occurrence_threshold} \
      --outputstats          ${taxon}_stats.txt \
      --output_large_occ     ${taxon}_Occurrences_large.txt \
      --output_low_occ       ${taxon}_Occurrences_small.txt

    if [ ${params.cleanupwd} == true ]; then
      echo "Cleaning up"
      rm -r ${phylotree_initial} ${phylotree_processed}
    fi
    """
}

// Define species set (batched)
process define_species_set_batched {

    publishDir "${OUT_0_TRR}", pattern: 'results/*_stats.txt', saveAs: { filename -> file(filename).name }, mode: "${params.publish_dir_mode}", overwrite: true

    input:
      path(occ_counts_csv)
      path(phylotrees_table)
      path(phylotree_initial,   stageAs: "t1/*")
      path(phylotree_processed, stageAs: "t2/*")

    output:
      path "results/*_Occurrences_large.txt",  emit: occ_large
      path "results/*_Occurrences_small.txt",  emit: occ_small
      path "results/*_stats.txt",              emit: stats

    script:
    """
    echo -e "Defining species set [BATCHED]\n"
    echo "Occurrences counts: "   ${occ_counts_csv}
    echo "Tree table (by taxon): "${phylotrees_table}
    echo "Processed trees: "      ${phylotree_processed}
    echo "Initial trees: "        ${phylotree_initial}
    echo "Occurrence threshold: " ${params.outlier_occurrence_threshold}

    ## Remove paths from the tree table
    awk -F'\t' 'BEGIN { OFS=FS } { sub(".*/", "", \$2); sub(".*/", "", \$3); print }' \
      ${phylotrees_table} \
      > tree_table_no_paths.txt

    echo -e "\nCounting the number of GBIF records represented and not represented in a phylogenetic tree\n"

    mkdir results

    parallel -j1 -a tree_table_no_paths.txt --colsep '\t' \
      "echo {1} && \
      gbif_records_in_tree.R \
        --gbif                 ${occ_counts_csv} \
        --sourcetree           t1/{2} \
        --tree                 t2/{3} \
        --occurrence_threshold ${params.outlier_occurrence_threshold} \
        --outputstats          results/{1}_stats.txt \
        --output_large_occ     results/{1}_Occurrences_large.txt \
        --output_low_occ       results/{1}_Occurrences_small.txt"

    if [ ${params.cleanupwd} == true ]; then
      echo "Cleaning up"
      rm -r t1 t2
    fi
    """
}


// Pool species lists from different taxonomic groups
// In different trees, there might the same species (e.g., outgroup species or if trees overlap),
// therefore, we need to keep them only once in the pooled species list
// In addition, we can exclude extinct species
process pool_species_lists {

    input:
      path(occurrence_counts, stageAs: "occ_counts/*")
      path(extinct_taxa)

    output:
      path "Occurrences_large.txt", emit: occ_large
      path "Occurrences_small.txt", emit: occ_small

    script:
    def extinctArg = extinct_taxa ? "--extinct ${extinct_taxa}" : ''
    """
    echo -e "Pooling species lists\n"

    echo "Occurrence counts: " ${occurrence_counts}
    echo "Extinct taxa: "      ${extinct_taxa}

    pool_occurrence_counts.R \
      --inpdir        occ_counts/ \
      --output_large  Occurrences_large.txt \
      --output_small  Occurrences_small.txt \
      ${extinctArg}

    if [ ${params.cleanupwd} == true ]; then
      echo "Cleaning up"
      rm -r occ_counts
    fi
    """
}

// Prepare data for spatial outlier removal (atomic, per species)
process prepare_species {

    tag "${species}"

    input:
      tuple val(species), path(occurrences)

    output:
      tuple val(species), path("*.csv.gz"),  emit: h3_binned_csv
      path("*.parquet"),                     emit: h3_binned

    script:
    def tempDirArg = task.ext.tempDir ? "-x ${task.ext.tempDir}" : ""
    def memoryArg  = task.memory  ? "-m ${task.memory.toMega()}.MB" : ""
    def basisOfRecordArg = params.basis_of_record ? "-b ${params.basis_of_record}" : ""
    def duckdbArg = (workflow.containerEngine == 'singularity') ? '-e "/usr/local/bin/duckdb_ext"' : ''
    def normalized_species = species.replaceAll(/[^a-zA-Z0-9]+/, '_')
    """
    echo -e "Preparing data for spatial outlier removal\n"

    echo "Species occurrences: " ${occurrences}
    echo "Species name: "        ${species}
    echo "Output file: "         ${normalized_species}.parquet
    echo "H3 resolution: "       ${params.outlier_h3_resolution}
    if [ -n "${task.memory}"  ];     then echo "Memory: ${memoryArg}";          fi
    if [ -n "${task.ext.tempDir}" ]; then echo "Temp directory: ${tempDirArg}"; fi
    echo "Container engine: "    ${workflow.containerEngine}

    h3_preprocessor.sh \
      -i ${occurrences}'/*' \
      -s "${species}" \
      -o "${normalized_species}.parquet" \
      -r ${params.outlier_h3_resolution} \
      -t ${task.cpus} \
      ${memoryArg} ${tempDirArg} \
      ${basisOfRecordArg} \
      ${duckdbArg}


    if [ ${params.cleanupwd} == true ]; then
      echo "Cleaning up"
      rm -r ${occurrences}
    fi

    echo "..Done"
    """
}


// Prepare data for spatial outlier removal, in batches
process prepare_species_batched {

    input:
      tuple val(species), path(occurrences)

    output:
      path("prepare_species/*.csv.gz"),  emit: h3_binned_csv
      path("prepare_species/*.parquet"), emit: h3_binned

    script:
    def tempDirArg = task.ext.tempDir ? "-x ${task.ext.tempDir}" : ""
    def memoryArg  = task.memory  ? "-m ${task.memory.toMega()}.MB" : ""
    def basisOfRecordArg = params.basis_of_record ? "-b ${params.basis_of_record}" : ""
    def duckdbArg = (workflow.containerEngine == 'singularity') ? '-e "/usr/local/bin/duckdb_ext"' : ''
    """
echo -e "Preparing data for spatial outlier removal [BATCHED]\n"

echo "Species occurrences: " ${occurrences}
echo "H3 resolution: " ${params.outlier_h3_resolution}

echo -e "\n"
echo "Species: ${species}"
echo -e "\n"

if [ -n "${task.memory}"  ];     then echo "Memory: ${memoryArg}";          fi
if [ -n "${task.ext.tempDir}" ]; then echo "Temp directory: ${tempDirArg}"; fi
echo "Container engine: " ${workflow.containerEngine}

echo "..Writing species to file"
cat << EOF > ids
${species.join("\n")}
EOF

mkdir -p prepare_species

echo "..Processing species"
parallel -j1 -a ids --halt now,fail=1 \
  --rpl '{normalized_species} s/[^[:alnum:]]/_/g' \
  "echo {} && \
  h3_preprocessor.sh \
    -i ${occurrences}'/*' \
    -s {} \
    -o prepare_species/{normalized_species}.parquet \
    -r ${params.outlier_h3_resolution} \
    -t ${task.cpus} \
    ${memoryArg} ${tempDirArg} \
    ${basisOfRecordArg} \
    ${duckdbArg}"

if [ ${params.cleanupwd} == true ]; then
  echo "Cleaning up"
  rm ids
  rm -r ${occurrences}
  rm prepare_species/*.sql
fi

echo -e "\n..All done"
    """
}


// Detect spatial outliers (per species) using DBSCAN
process dbscan {
    // cpus 2

    tag "${species}"

    input:
      tuple val(species), path(h3_binned_csv)

    output:
      tuple val(species), path("*.DBSCAN-outliers.txt.gz"), emit: dbscan_scores

    script:
    """
    echo -e "Detecting spatial outliers using DBSCAN\n"
    echo "Species: "    ${species}

    ## Normalize species name for file naming
    normalized_species=\$(echo "${species}" | sed 's/[^[:alnum:]]/_/g')

    elki_outlier.sh \
      --input  ${h3_binned_csv} \
      --output "\${normalized_species}_DBSCAN-scores.txt.gz" \
      --method DBSCANOutlierDetection \
      --geomodel WGS84SpheroidEarthModel \
      --indextype RStarTree \
      --k ${params.dbscan_minpts} \
      --d ${params.dbscan_eps}

    echo "ELKI finished"
    echo "Merging H3 cell IDs and scores"

    ## Merge H3 cell IDs and scores
    paste \
      <(zcat "${h3_binned_csv}" | cut -d',' -f3-4 | sed 's/,/\t/g') \
      <(zcat "\${normalized_species}_DBSCAN-scores.txt.gz" | cut -d' ' -f2) \
      | gzip -3 \
      > "\${normalized_species}.DBSCAN-outliers.txt.gz"

    if [ ${params.cleanupwd} == true ]; then
      echo "Cleaning up"
      rm "\${normalized_species}_DBSCAN-scores.txt.gz"
      rm "${h3_binned_csv}"
    fi

    echo "..Done"
    """
}


// Detect spatial outliers using DBSCAN (batched)
process dbscan_batched {
    // cpus 2

    input:
      path(h3_binned_csv, stageAs: "inp/*")

    output:
      path("mrg/*.DBSCAN-outliers.txt.gz"), emit: dbscan_scores

    script:
    """
echo -e "Detecting spatial outliers using DBSCAN [BATCHED]\n"
echo "Number of input files detected: " \$(find inp -name '*.csv.gz' | wc -l)

echo "..Processing species occurrences"
mkdir -p dbs
find inp -name '*.csv.gz' \
  | parallel -j1 --halt now,fail=1 \
    --rpl '{fname} s:^.*/:: ; s/.csv.gz//' \
    "echo {} && \
    elki_outlier.sh \
      --input  {} \
      --output dbs/{fname}_DBSCAN-scores.txt.gz \
      --method DBSCANOutlierDetection \
      --geomodel WGS84SpheroidEarthModel \
      --indextype RStarTree \
      --k ${params.dbscan_minpts} \
      --d ${params.dbscan_eps}"

echo "ELKI finished"
echo "Merging H3 cell IDs and scores"

## Merge H3 cell ID, species name, and outlier scores
mkdir -p mrg
find inp -name '*.csv.gz' \
  | parallel -j1 --halt now,fail=1 \
    --rpl '{fname} s:^.*/:: ; s/.csv.gz//' \
    "echo {} && \
    paste \
      <(zcat "{}" | cut -d',' -f3-4 | sed 's/,/\t/g') \
      <(zcat "dbs/{fname}_DBSCAN-scores.txt.gz" | cut -d' ' -f2) \
      | gzip -3 \
      > mrg/{fname}.DBSCAN-outliers.txt.gz"

## Clean up
if [ ${params.cleanupwd} == true ]; then
  echo "Cleaning up"
  rm -r inp
  rm -r dbs
fi

echo -e "\n..All done"
    """
}


// Count number of grid cells identified as outliers
process count_outliers {

  publishDir "${OUT_2_OUT}", mode: "${params.publish_dir_mode}", overwrite: true

  input:
    path(outlier_scores, stageAs: "scores/*")

  output:
    path("outlier_scores.txt.gz"), emit: outlier_scores
    path("outlier_scores_summary.txt")

  script:
  def tempDirArg = task.ext.tempDir ? "-x ${task.ext.tempDir}" : ""
  def memoryArg  = task.memory  ? "-m ${task.memory.toMega()}.MB" : ""
  """
  echo -e "Outlier removal summary\n"
  if [ -n "${task.memory}"  ];     then echo "Memory: ${memoryArg}";          fi
  if [ -n "${task.ext.tempDir}" ]; then echo "Temp directory: ${tempDirArg}"; fi

  outlier_removal_summary.sh \
    -i \$(pwd)/scores \
    -o outlier_scores.txt.gz \
    -q 0.5 -t ${task.cpus} \
    ${memoryArg} ${tempDirArg}

  if [ ${params.cleanupwd} == true ]; then
    echo "Cleaning up"
    rm -r scores
  fi

  """
}


// Process DBSCAN results (remove grids marked as outliers)
// NB. Scores from `DBSCANOutlierDetection` are binary!
process process_dbscan {

    tag "${species}"

    input:
      tuple val(species), path(dbscan_scores), path(occ)

    output:
      tuple val("${species}"), path("*.parquet"), path("*.txt"), emit: nooutliers

    script:
    def tempDirArg = task.ext.tempDir ? "-x ${task.ext.tempDir}" : ""
    def memoryArg  = task.memory  ? "-m ${task.memory.toMega()}.MB" : ""
    def basisOfRecordArg = params.basis_of_record ? "-b ${params.basis_of_record}" : ""
    def duckdbArg = (workflow.containerEngine == 'singularity') ? '-e "/usr/local/bin/duckdb_ext"' : ''
    def normalized_species = species.replaceAll(/[^a-zA-Z0-9]+/, '_')
    """
    echo -e "Processing DBSCAN results\n"
    echo "Species: " ${species}
    if [ -n "${task.memory}"  ];     then echo "Memory: ${memoryArg}";          fi
    if [ -n "${task.ext.tempDir}" ]; then echo "Temp directory: ${tempDirArg}"; fi
    echo "Container engine: " ${workflow.containerEngine}

    remove_outliers_from_parquet.sh \
      -i ${occ}'/*' \
      -w ${dbscan_scores} \
      -r ${params.outlier_h3_resolution} \
      -o "${normalized_species}.parquet" \
      -s "${species}" \
      -t ${task.cpus} \
      ${memoryArg} ${tempDirArg} \
      ${basisOfRecordArg} \
      ${duckdbArg}

    ## For compatibility with low-occurrence species, we need a file with species names
    echo ${species} > "${normalized_species}.txt"

    if [ ${params.cleanupwd} == true ]; then
      echo "Cleaning up"
      rm ${dbscan_scores}
      if [ -d ${occ} ]; then
        rm -r ${occ}
      else
        rm ${occ}
      fi
    fi
    """
}


// Process DBSCAN results (remove grids marked as outliers, batched)
process process_dbscan_batched {

    input:
      tuple path(dbscan_scores, stageAs: "scores/*"), path(occ)

    output:
      path("flt/*.parquet"), emit: nooutliers

    script:
    def tempDirArg = task.ext.tempDir ? "-x ${task.ext.tempDir}" : ""
    def memoryArg  = task.memory  ? "-m ${task.memory.toMega()}.MB" : ""
    def basisOfRecordArg = params.basis_of_record ? "-b ${params.basis_of_record}" : ""
    def duckdbArg = (workflow.containerEngine == 'singularity') ? '-e "/usr/local/bin/duckdb_ext"' : ''
    """
echo -e "Processing DBSCAN results [BATCHED]\n"
if [ -n "${task.memory}"  ];     then echo "Memory: ${memoryArg}";          fi
if [ -n "${task.ext.tempDir}" ]; then echo "Temp directory: ${tempDirArg}"; fi
echo "Container engine: "  ${workflow.containerEngine}

echo "Number of input files detected: " \$(find scores -name '*.txt.gz' | wc -l)

echo -e "\n\n..Processing species occurrences"
mkdir -p flt

# Notes: species parameter is not provided to the `remove_outliers_from_parquet.sh`,
# so the first species in the outlier scores file will be used

find scores -name '*.txt.gz' \
  | parallel -j1 --halt now,fail=1 \
    --rpl '{fname} s:^.*/:: ; s/.txt.gz//' \
    "echo {} && \
    remove_outliers_from_parquet.sh \
      -i ${occ}'/*' \
      -w {} \
      -r ${params.outlier_h3_resolution} \
      -o flt/{fname}.parquet \
      -t ${task.cpus} \
      ${memoryArg} ${tempDirArg} \
      ${basisOfRecordArg} \
      ${duckdbArg} \
    && echo -e '\n\n'"


## Clean up
if [ ${params.cleanupwd} == true ]; then
  echo "Cleaning up"
  rm -r scores
  rm -r flt/*.sql
  rm -r ${occ}
fi
    """
}


// Filter and bin GBIF records to the H3 grid cells of final resolution
process filter_and_bin {

    tag "${name}"
    // publishDir "${OUT_1_PRQ}", mode: "${params.publish_dir_mode}"

    input:
      tuple val(name), path(occ), path(spp)

    output:
      path "*_filtered.parquet", emit: aggregated

    script:
    // Check if input is a directory (raw data for low-occurrence species)
    // or a single parquet file (outlier-removed data for large-occurrence species)
    def inp = occ.getName().endsWith('.parquet') ? occ : "${occ}'/*'"
    
    def tempDirArg = task.ext.tempDir ? "-x ${task.ext.tempDir}" : ""
    def memoryArg  = task.memory  ? "-m ${task.memory.toMega()}.MB" : ""
    def basisOfRecordArg = params.basis_of_record ? "-b ${params.basis_of_record}" : ""
    def duckdbArg = (workflow.containerEngine == 'singularity') ? '-e "/usr/local/bin/duckdb_ext"' : ''
    def normalized_species = name.replaceAll(/[^a-zA-Z0-9]+/, '_')
    """
    echo -e "Filtering and binning GBIF records\n"

    echo "Species occurrences: " ${inp}
    echo "Species: "             ${name}
    echo "Output file: "         ${normalized_species}_filtered.parquet
    echo "H3 resolution: "       ${params.h3_resolution}
    if [ -n "${task.memory}"  ];     then echo "Memory: ${memoryArg}";          fi
    if [ -n "${task.ext.tempDir}" ]; then echo "Temp directory: ${tempDirArg}"; fi
    echo "Container engine: "    ${workflow.containerEngine}

    filter_and_bin.sh \
      -i ${inp} \
      -s ${spp} \
      -o ${normalized_species}_filtered.parquet \
      -r ${params.h3_resolution} \
      -t ${task.cpus} \
      ${memoryArg} ${tempDirArg} \
      ${basisOfRecordArg} \
      ${duckdbArg}

    echo -e "\nFiltering and binning finished. Cleaning up\n"

    ## Clean up
    if [ ${params.cleanupwd} == true ]; then
      echo "Cleaning up"
      if [ -e ${occ} ]; then
        if [ -d ${occ} ]; then
          rm -r ${occ}
        else
          rm ${occ}
        fi
      fi
    fi
    """
}


// Filter and bin GBIF records to the H3 grid cells of final resolution (batched)
process filter_and_bin_batched {

    input:
      path(occ, stageAs: "flt/*")

    output:
      path "binned/*.parquet", emit: aggregated

    script:
    def tempDirArg = task.ext.tempDir ? "-x ${task.ext.tempDir}" : ""
    def memoryArg  = task.memory  ? "-m ${task.memory.toMega()}.MB" : ""
    def basisOfRecordArg = params.basis_of_record ? "-b ${params.basis_of_record}" : ""
    def duckdbArg = (workflow.containerEngine == 'singularity') ? '-e "/usr/local/bin/duckdb_ext"' : ''
    """
    echo -e "Filtering and binning GBIF records [BATCHED]\n"

    echo "Filtered species occurrences: " ${occ}
    echo "H3 resolution: " ${params.h3_resolution}
    if [ -n "${task.memory}"  ];     then echo "Memory: ${memoryArg}";          fi
    if [ -n "${task.ext.tempDir}" ]; then echo "Temp directory: ${tempDirArg}"; fi
    echo "Container engine: "    ${workflow.containerEngine}

    ## Create files with species keys (for compatibility with low-occurrence species filtering)
    mkdir -p spkeys
    find flt -name '*.parquet' \
      | parallel -j1 --halt now,fail=1 \
        --rpl '{/:} s:(.*/)?([^/.]+)(\\.[^/]+)*\$:\$2:' \
        "echo {/:} > spkeys/{/:}.txt"

    ## Filter and bin occurrences
    mkdir -p binned
    find flt -name '*.parquet' \
      | parallel -j1 --halt now,fail=1 \
        --rpl '{/:} s:(.*/)?([^/.]+)(\\.[^/]+)*\$:\$2:' \
        "filter_and_bin.sh \
          -i {} \
          -s spkeys/{/:}.txt \
          -o binned/{/:}.parquet \
          -r ${params.h3_resolution} \
          -t ${task.cpus} \
          ${memoryArg} ${tempDirArg} \
          ${basisOfRecordArg} \
          ${duckdbArg}"

    ## Clean up
    if [ ${params.cleanupwd} == true ]; then
      echo -e "\nFiltering and binning finished. Cleaning up\n"
      rm -r spkeys
      rm -r flt
      rm -r binned/*.sql
    fi

    echo -e "\n..All done\n"
    """
}


// Merge multiple parquet files (for large-occurrence species) into a single file
process pool_parquets {

  publishDir "${OUTDIR}/", mode: "${params.publish_dir_mode}", overwrite: true

  input:
    path(parquets, stageAs: "parquets_staged/*")
  
  output:
    path("Filtered_parquet/*.parquet"), emit: prq

  script:
  def tempDirArg = task.ext.tempDir ? "-x ${task.ext.tempDir}" : ""
  def memoryArg  = task.memory  ? "-m ${task.memory.toMega()}.MB" : ""
  """
  echo "Merging parquet files"

  if [ -n "${task.memory}"  ];     then echo "Memory: ${memoryArg}";          fi
  if [ -n "${task.ext.tempDir}" ]; then echo "Temp directory: ${tempDirArg}"; fi

  pool_parquets.sh \
    -i parquets_staged \
    -o Filtered_parquet \
    -t ${task.cpus} \
    ${memoryArg} ${tempDirArg}

  echo "Pooling parquet files finished"

  ## Clean up
  if [ ${params.cleanupwd} == true ]; then
    echo "Cleaning up"
    rm -r parquets_staged
  fi
  """
}



workflow {
  if (params.task_batching == true) {
    task_batching()
  } else {
    atomic_tasks()
  }
}

// Workflow with atomic tasks
workflow atomic_tasks {

  // Input directory with GBIF species occurrences (parquet files)
  ch_occurrence_dir = Channel.fromPath(params.occurrences)

  // Count species occurrences
  count_occurrences(ch_occurrence_dir)

  // Phylogenetic trees
  ch_phylotrees = Channel.fromList(phylo_trees)

  // Filter out missing files
  ch_phylotrees = ch_phylotrees.filter { name, file1, file2 ->  
    def file2_obj = file(file2)
    def exists = file2_obj.exists()
    if (!exists) {
        println "Record removed: $name (missing file: $file2)"
    }
    return exists
  }

  // Combine phylogenetic trees and occurrence counts
  // tuple(taxon, phylotree_initial, phylotree_processed, occ_counts_csv)
  ch_phylo_occ = ch_phylotrees.combine(count_occurrences.out.occ_counts_csv)

  // Define species set, split into large and small occurrence sets
  define_species_set(ch_phylo_occ)

  // Channel with occurrence counts (collected over different taxa)
  ch_occcounts = define_species_set.out.occ_large
    .merge(define_species_set.out.occ_small)
    .collect()
  
  // Extinct taxa (optional, if provided)
  ch_extinct_taxa = params.extinct_taxa ? Channel.fromPath(params.extinct_taxa) : Channel.empty()

  // Pool species lists from different taxonomic groups
  pool_species_lists(ch_occcounts, ch_extinct_taxa)

  // Run outlier removal for "large" species
  if (params.nodbscan == false) {

    // Channel with species for spatial outlier removal
    // NB. results returned by `splitText` operator are always terminated by a `\n` newline character, so we need to trim it
    ch_large_species = pool_species_lists.out.occ_large
      .splitText()
      .map{ it -> it.trim() }
      .combine(ch_occurrence_dir)

    // Prepare data for spatial outlier removal (large species only)
    prepare_species(ch_large_species)

    // Remove spatial outliers using DBSCAN
    dbscan(prepare_species.out.h3_binned_csv)

    // Add raw occurrence path to the channel with DBSCAN scores
    // tuple(species, dbscan_scores, occurrence_path)
    dbscan.out.dbscan_scores
      .combine(ch_occurrence_dir)
      .set { ch_dbscan_scores }

    // Process DBSCAN results (remove grids marked as outliers)
    process_dbscan(ch_dbscan_scores)

    // Outlier removal summary
    dbscan.out.dbscan_scores
      .map { it -> it[1] }
      .collect()
      .set { ch_all_scores }
  
    count_outliers(ch_all_scores)

    // Data for low-occurrence filtering   tuple("low", raw_data, low_occ_species )
    ch_occurrence_dir
      .merge(pool_species_lists.out.occ_small) { occ, spp -> tuple("low", occ, spp) }
      .set { ch_spk_low }

    // Data for large-occurrence filtering   tuple( species, outlier_removed_data, large_occ_species )
    ch_spk_large = process_dbscan.out.nooutliers

    ch_spk_low
      .concat(ch_spk_large)
      .set { ch_spk }

    // end of outlier removal

  } else {

    // Prepare tuples for filtering and binning of species (without outlier removal)
    ch_occurrence_dir
      .merge(pool_species_lists.out.occ_small) { occ, spp -> tuple("low", occ, spp) }
      .set { ch_spk_low }

    ch_occurrence_dir
      .merge(pool_species_lists.out.occ_large) { occ, spp -> tuple("large", occ, spp) }
      .set { ch_spk_large }

    // tuple("low",   raw_data, low_occ_species )
    // tuple("large", raw_data, large_occ_speciess )

    ch_spk_low
      .concat(ch_spk_large)
      .set { ch_spk }

  }

  // Filter and bin occurrences
  filter_and_bin(ch_spk)

  // Merge parquet files into a bigger chunks (300k records per file)
  ch_all_filtered = filter_and_bin.out.aggregated.collect()
  pool_parquets(ch_all_filtered)

}


// Workflow with task batching
workflow task_batching {

  // Input directory with GBIF species occurrences (parquet files)
  ch_occurrence_dir = Channel.fromPath(params.occurrences)

  // Count species occurrences
  count_occurrences(ch_occurrence_dir)

  // Phylogenetic trees
  ch_phylotrees = Channel.fromList(phylo_trees)

  // Collect phylogenetic trees into a table (since trees can be named in non-regular way)
  ch_phylotrees
    .map { tuple -> tuple.join('\t') }
    .collectFile(name: 'table.txt', newLine: true)
    .set { ch_phylotrees_table }

  ch_trees_1 = ch_phylotrees.map { tuple -> tuple[1] }.collect()  // initial trees
  ch_trees_2 = ch_phylotrees.map { tuple -> tuple[2] }.collect()  // processed trees

  // Define species set, split into large and small occurrence sets
  define_species_set_batched(
    count_occurrences.out.occ_counts_csv,
    ch_phylotrees_table,
    ch_trees_1,
    ch_trees_2
  )

  // Channel with occurrence counts (collected over different taxa)
  ch_occcounts = define_species_set_batched.out.occ_large
    .merge(define_species_set_batched.out.occ_small)
    .collect()
  
  // Extinct taxa (optional, if provided)
  ch_extinct_taxa = params.extinct_taxa ? Channel.fromPath(params.extinct_taxa) : Channel.empty()

  // Pool species lists from different taxonomic groups
  pool_species_lists(ch_occcounts, ch_extinct_taxa)

  // Run outlier removal for "large" species
  if (params.nodbscan == false) {

    // Channel with species keys for spatial outlier removal
    // NB. results returned by `splitText` operator are always terminated by a `\n` newline character, so we need to trim it
    // Output: [ [1, 2, 3], glob_occ ]
    ch_large_species = pool_species_lists.out.occ_large
      .splitText()
      .map{ it -> it.trim() }
      .buffer(size: params.batchsize, remainder: true)
      .map { it -> [it] }
      .combine(ch_occurrence_dir)

    // Prepare data for spatial outlier removal (large species only)
    prepare_species_batched(ch_large_species)

    // Remove spatial outliers using DBSCAN
    dbscan_batched(prepare_species_batched.out.h3_binned_csv)

    // Add raw occurrence path to the channel with DBSCAN scores
    // Output: [ [ out1, out2, out3 ], glob_occ ]
    dbscan_batched.out.dbscan_scores
      .map { it -> [it] }
      .combine(ch_occurrence_dir)
      .set { ch_dbscan_scores }

    // Process DBSCAN results (remove grids marked as outliers, batched)
    process_dbscan_batched(ch_dbscan_scores)

    // Outlier removal summary
    dbscan_batched.out.dbscan_scores
      .collect()
      .set { ch_all_scores }

    count_outliers(ch_all_scores)

    // Low-occurrence filtering   tuple("low", raw_data, low_occ_specieskeys )
    ch_occurrence_dir
      .merge(pool_species_lists.out.occ_small) { occ, spp -> tuple("low", occ, spp) }
      .set { ch_spk_low }

    filter_and_bin(ch_spk_low)

    // Large-occurrence filtering
    filter_and_bin_batched(process_dbscan_batched.out.nooutliers)
    
    // Merge parquet files into a bigger chunks (300k records per file)
    filter_and_bin_batched.out.aggregated
      .concat(filter_and_bin.out.aggregated)
      .collect()
      .set { ch_all_filtered }

    pool_parquets(ch_all_filtered)

    // end of outlier removal

  } else {

    // Prepare tuples for filtering and binning of species (without outlier removal)
    ch_occurrence_dir
      .merge(pool_species_lists.out.occ_small) { occ, spp -> tuple("low", occ, spp) }
      .set { ch_spk_low }

    ch_occurrence_dir
      .merge(pool_species_lists.out.occ_large) { occ, spp -> tuple("large", occ, spp) }
      .set { ch_spk_large }

    // tuple("low",   raw_data, low_occ_specieskeys )
    // tuple("large", raw_data, large_occ_specieskeys )

    ch_spk_low
      .concat(ch_spk_large)
      .set { ch_spk }

    // Filter and bin occurrences
    filter_and_bin(ch_spk)

    // Merge parquet files into a bigger chunks (300k records per file)
    ch_all_filtered = filter_and_bin.out.aggregated.collect()
    pool_parquets(ch_all_filtered)

    // end of no-outlier-removal
  }

}


// On completion
workflow.onComplete {
    println "Pipeline completed at : $workflow.complete"
    println "Duration              : ${workflow.duration}"
    println "Execution status      : ${workflow.success ? 'All done!' : 'Failed' }"
}

// On error
workflow.onError {
    println "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

