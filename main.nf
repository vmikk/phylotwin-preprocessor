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
    tuple("Plants",      "${params.t1}/Tietje_2023_Seed_plants_TACT.tre.gz",  "${params.t2}/Tietje_2023_Seed_plants_TACT.nwk"),
    tuple("Ferns",       "${params.t1}/FTOL_1-7-0_sanger_con_dated.tre.gz",   "${params.t2}/Ferns_FTOL_1-7-0.nwk"),
    tuple("Mushrooms",   "${params.t1}/Varga_2019_Mushrooms.tre.gz",          "${params.t2}/Varga_2019_Mushrooms.nwk"),
    tuple("Mammals",     "${params.t1}/VertLife_Mammals.tre.gz",              "${params.t2}/VertLife_Mammals.nwk"),
    tuple("Birds",       "${params.t1}/VertLife_Birds.tre.gz",                "${params.t2}/VertLife_Birds.nwk"),
    tuple("Squamates",   "${params.t1}/VertLife_Squamates.tre.gz",            "${params.t2}/VertLife_Squamates.nwk"),
    tuple("Amphibians",  "${params.t1}/VertLife_Amphibians.tre.gz",           "${params.t2}/VertLife_Amphibians.nwk"),
    tuple("Fish",        "${params.t1}/FishTree_actinopt_12k_treePL.tre.xz",  "${params.t2}/FishTree.nwk"),
    tuple("Butterflies", "${params.t1}/Kawahara_2023_Butterflies_Newick.tre",          "${params.t2}/Kawahara_2023_Butterflies.nwk"),
    tuple("Bees",        "${params.t1}/Henriquez-Piskulich_2024_BeetreeOfLife.tre.gz", "${params.t2}/Henriquez-Piskulich_2024_BeetreeOfLife.nwk"),
]

// GTDB_ar53_r220.tre.gz
// GTDB_bac120_r220.tre.gz


// Directory for publishing outputs
OUTDIR = params.userid ? params.outdir + "/" + params.userid : params.outdir
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

    echo "..Done"
    """
}


// Define species set (given phylogenetic trees available)
process define_species_set {

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
    pool_occurrence_counts.R \
      --inpdir        occ_counts/ \
      --output_large  Occurrences_large.txt \
      --output_small  Occurrences_small.txt \
      ${extinctArg}
    """
}

// Prepare data for spatial outlier removal (atomic, per species)
process prepare_species {

    tag "${specieskey}"

    input:
      tuple val(specieskey), path(occurrences)

    output:
      tuple val(specieskey), path("${specieskey}.csv.gz"),  emit: h3_binned_csv
      path("${specieskey}.parquet"),                        emit: h3_binned

    script:
    def tempDirArg = task.ext.tempDir ? "-x ${task.ext.tempDir}" : ""
    def memoryArg  = task.memory  ? "-m ${task.memory.toMega()}.MB" : ""
    def basisOfRecordArg = params.basis_of_record ? "-b ${params.basis_of_record}" : ""
    def duckdbArg = (workflow.containerEngine == 'singularity') ? '-e "/usr/local/bin/duckdb_ext"' : ''
    """
    echo -e "Preparing data for spatial outlier removal\n"

    echo "Species occurrences: " ${occurrences}
    echo "Species key: "         ${specieskey}
    echo "Output file: "         ${specieskey}.parquet
    echo "H3 resolution: "       ${params.outlier_h3_resolution}
    if [ -n "${task.memory}"  ];     then echo "Memory: ${memoryArg}";          fi
    if [ -n "${task.ext.tempDir}" ]; then echo "Temp directory: ${tempDirArg}"; fi
    echo "Container engine: "    ${workflow.containerEngine}

    h3_preprocessor.sh \
      -i ${occurrences}'/*' \
      -s ${specieskey} \
      -o ${specieskey}.parquet \
      -r ${params.outlier_h3_resolution} \
      -t ${task.cpus} \
      ${memoryArg} ${tempDirArg} \
      ${basisOfRecordArg} \
      ${duckdbArg}

    echo "..Done"
    """
}


// Prepare data for spatial outlier removal, in batches
process prepare_species_batched {

    input:
      tuple val(specieskeys), path(occurrences)

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
echo "Species keys: ${specieskeys}"
echo -e "\n"

if [ -n "${task.memory}"  ];     then echo "Memory: ${memoryArg}";          fi
if [ -n "${task.ext.tempDir}" ]; then echo "Temp directory: ${tempDirArg}"; fi
echo "Container engine: " ${workflow.containerEngine}

echo "..Writing species keys to file"
cat << EOF > ids
${specieskeys.join("\n")}
EOF

mkdir -p prepare_species

echo "..Processing species keys"
parallel -j1 -a ids --halt now,fail=1 \
  "echo {} && \
  h3_preprocessor.sh \
    -i ${occurrences}'/*' \
    -s {} \
    -o prepare_species/{}.parquet \
    -r ${params.outlier_h3_resolution} \
    -t ${task.cpus} \
    ${memoryArg} ${tempDirArg} \
    ${basisOfRecordArg} \
    ${duckdbArg}"

## Clean up
rm ids
rm -r ${occurrences}
rm prepare_species/*.sql

echo -e "\n..All done"

    """
}


// Detect spatial outliers (per species) using DBSCAN
process dbscan {
    // cpus 2

    tag "${specieskey}"

    input:
      tuple val(specieskey), path(h3_binned_csv)

    output:
      tuple val(specieskey), path("${specieskey}.DBSCAN-outliers.txt.gz"), emit: dbscan_scores

    script:
    """
    echo -e "Detecting spatial outliers using DBSCAN\n"
    echo "Species key: "    ${specieskey}

    elki_outlier.sh \
      --input  ${h3_binned_csv} \
      --output ${specieskey}_DBSCAN-scores.txt.gz \
      --method DBSCANOutlierDetection \
      --geomodel WGS84SpheroidEarthModel \
      --indextype RStarTree \
      --k ${params.dbscan_minpts} \
      --d ${params.dbscan_eps}

    echo "ELKI finished"
    echo "Merging H3 cell IDs and scores"

    ## Merge H3 cell IDs and scores
    paste \
      <(zcat "${h3_binned_csv}" | cut -d',' -f3) \
      <(zcat "${specieskey}_DBSCAN-scores.txt.gz" | cut -d' ' -f2) \
      | gzip -3 \
      > "${specieskey}.DBSCAN-outliers.txt.gz"

    rm "${specieskey}_DBSCAN-scores.txt.gz"

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
    --rpl '{/:} s:(.*/)?([^/.]+)(\\.[^/]+)*\$:\$2:' \
    "echo {} && \
    elki_outlier.sh \
      --input  {} \
      --output dbs/{/:}_DBSCAN-scores.txt.gz \
      --method DBSCANOutlierDetection \
      --geomodel WGS84SpheroidEarthModel \
      --indextype RStarTree \
      --k ${params.dbscan_minpts} \
      --d ${params.dbscan_eps}"

echo "ELKI finished"
echo "Merging H3 cell IDs and scores"

## Merge SpeciesKey, H3 cell IDs, and scores
mkdir -p mrg
find inp -name '*.csv.gz' \
  | parallel -j1 --halt now,fail=1 \
    --rpl '{/:} s:(.*/)?([^/.]+)(\\.[^/]+)*\$:\$2:' \
    "echo {} && \
    paste \
      <(zcat "{}" | cut -d',' -f3) \
      <(zcat "dbs/{/:}_DBSCAN-scores.txt.gz" | cut -d' ' -f2) \
      | gzip -3 \
      > mrg/{/:}.DBSCAN-outliers.txt.gz"

## Clean up
rm -r inp
rm -r dbs/*DBSCAN-scores.txt.gz

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
    -q 0.5 -t 2

  """
}


// Process DBSCAN results (remove grids marked as outliers)
// NB. Scores from `DBSCANOutlierDetection` are binary!
process process_dbscan {

    tag "${specieskey}"

    input:
      tuple val(specieskey), path(dbscan_scores), path(occ)

    output:
      tuple val("${specieskey}"), path("${specieskey}.parquet"), path("${specieskey}.txt"), emit: nooutliers

    script:
    def tempDirArg = task.ext.tempDir ? "-x ${task.ext.tempDir}" : ""
    def memoryArg  = task.memory  ? "-m ${task.memory.toMega()}.MB" : ""
    def basisOfRecordArg = params.basis_of_record ? "-b ${params.basis_of_record}" : ""
    def duckdbArg = (workflow.containerEngine == 'singularity') ? '-e "/usr/local/bin/duckdb_ext"' : ''
    """
    echo -e "Processing DBSCAN results\n"
    echo "Species key: " ${specieskey}
    if [ -n "${task.memory}"  ];     then echo "Memory: ${memoryArg}";          fi
    if [ -n "${task.ext.tempDir}" ]; then echo "Temp directory: ${tempDirArg}"; fi
    echo "Container engine: " ${workflow.containerEngine}

    remove_outliers_from_parquet.sh \
      -i ${occ}'/*' \
      -w ${dbscan_scores} \
      -r ${params.outlier_h3_resolution} \
      -o ${specieskey}.parquet \
      -s ${specieskey} \
      -t ${task.cpus} \
      ${memoryArg} ${tempDirArg} \
      ${basisOfRecordArg} \
      ${duckdbArg}

    ## For compatibility with low-occurrence species, we need a file with species keys
    echo ${specieskey} > "${specieskey}.txt"

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
find scores -name '*.txt.gz' \
  | parallel -j1 --halt now,fail=1 \
    --rpl '{/:} s:(.*/)?([^/.]+)(\\.[^/]+)*\$:\$2:' \
    "echo {} && \
    remove_outliers_from_parquet.sh \
      -i ${occ}'/*' \
      -w {} \
      -r ${params.outlier_h3_resolution} \
      -o flt/{/:}.parquet \
      -s {/:} \
      -t ${task.cpus} \
      ${memoryArg} ${tempDirArg} \
      ${basisOfRecordArg} \
      ${duckdbArg} \
    && echo -e '\n\n'"

    """
}


// Filter and bin GBIF records to the H3 grid cells of final resolution
process filter_and_bin {

    tag "${name}"
    // publishDir "${OUT_1_PRQ}", mode: "${params.publish_dir_mode}"

    input:
      tuple val(name), path(occ), path(spp)

    output:
      path "${name}_filtered.parquet", emit: aggregated

    script:
    // Check if input is a directory (raw data for low-occurrence species)
    // or a single parquet file (outlier-removed data for large-occurrence species)
    def inp = occ.getName().endsWith('.parquet') ? occ : "${occ}'/*'"
    
    def tempDirArg = task.ext.tempDir ? "-x ${task.ext.tempDir}" : ""
    def memoryArg  = task.memory  ? "-m ${task.memory.toMega()}.MB" : ""
    def basisOfRecordArg = params.basis_of_record ? "-b ${params.basis_of_record}" : ""
    def duckdbArg = (workflow.containerEngine == 'singularity') ? '-e "/usr/local/bin/duckdb_ext"' : ''
    """
    echo -e "Filtering and binning GBIF records\n"

    echo "Species occurrences: " ${inp}
    echo "Species keys: "        ${name}
    echo "Output file: "         ${name}_filtered.parquet
    echo "H3 resolution: "       ${params.h3_resolution}
    if [ -n "${task.memory}"  ];     then echo "Memory: ${memoryArg}";          fi
    if [ -n "${task.ext.tempDir}" ]; then echo "Temp directory: ${tempDirArg}"; fi
    echo "Container engine: "    ${workflow.containerEngine}

    filter_and_bin.sh \
      -i ${inp} \
      -s ${spp} \
      -o ${name}_filtered.parquet \
      -r ${params.h3_resolution} \
      -t ${task.cpus} \
      ${memoryArg} ${tempDirArg} \
      ${basisOfRecordArg} \
      ${duckdbArg}


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
  """
  echo "Merging parquet files"

  pool_parquets.sh \
    -i parquets_staged \
    -o Filtered_parquet

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

  // Channel with species keys for spatial outlier removal
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
  // tuple(specieskey, dbscan_scores, occurrence_path)
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

  // Data for low-occurrence filtering   tuple( raw_data, low_occ_specieskeys )
  ch_occurrence_dir
    .merge(pool_species_lists.out.occ_small) { occ, spp -> tuple("low", occ, spp) }
    .set { ch_spk_low }

  // Data for large-occurrence filtering   tuple( outlier_removed_data, large_occ_specieskey )
  ch_spk_large = process_dbscan.out.nooutliers
  
  ch_spk_low
    .concat(ch_spk_large)
    .set { ch_spk }

  // Filter and bin occurrences
  filter_and_bin(ch_spk)

  // Merge parquet files into a bigger chunks (300k records per file)
  ch_all_filtered = filter_and_bin.out.aggregated.collect()
  pool_parquets(ch_all_filtered)

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

