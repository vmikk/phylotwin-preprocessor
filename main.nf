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


// Count number of occurrences per species
process count_occurrences {

    input:
      path occurrences

    output:
      path "Occurrences_counts.parquet", emit: occ_counts_parquet
      path "Occurrences_counts.csv.gz",  emit: occ_counts_csv

    script:
    def tempDirArg = task.tempDir ? "-x ${task.tempDir}" : ""
    def memoryArg  = task.memory  ? "-m ${task.memory}"  : ""
    def basisOfRecordArg = params.basis_of_record ? "-b ${params.basis_of_record}" : ""
    """
    echo -e "Counting species occurrences\n"

    echo "Input directory: " ${occurrences}
    echo "CPUs: " ${task.cpus}
    if [ ! -z ${memoryArg}  ]; then echo "Memory: ${memoryArg}";          fi
    if [ ! -z ${tempDirArg} ]; then echo "Temp directory: ${tempDirArg}"; fi

    gbif_records_count.sh \
      -i ${occurrences} \
      -o "Occurrences_counts.parquet" \
      -t ${task.cpus} \
      "${memoryArg}" "${tempDirArg}" \
      "${basisOfRecordArg}"

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
process pool_species_lists {

    input:
      path(occurrence_counts, stageAs: "occ_counts/*")

    output:
      path "Occurrences_large.txt", emit: occ_large
      path "Occurrences_small.txt", emit: occ_small

    script:
    """
    pool_occurrence_counts.R \
      --inpdir        occ_counts/ \
      --output_large  Occurrences_large.txt \
      --output_small  Occurrences_small.txt
    """
}

// Prepare data for spatial outlier removal (per species)
process prepare_species {

    tag "${specieskey}"

    input:
      tuple val(specieskey), path(occurrences)

    output:
      tuple val(specieskey), path("${specieskey}.csv.gz"),  emit: h3_binned_csv
      path("${specieskey}.parquet"),                        emit: h3_binned

    script:
    def tempDirArg = task.tempDir ? "-x ${task.tempDir}" : ""
    def memoryArg  = task.memory  ? "-m ${task.memory}"  : ""
    def basisOfRecordArg = params.basis_of_record ? "-b ${params.basis_of_record}" : ""
    """
    echo -e "Preparing data for spatial outlier removal\n"

    echo "Species occurrences: " ${occurrences}
    echo "Species key: "         ${specieskey}
    echo "Output file: "         ${specieskey}.parquet
    echo "H3 resolution: "       ${params.outlier_h3_resolution}
    if [ ! -z ${memoryArg} ];  then echo "Memory: ${memoryArg}";          fi
    if [ ! -z ${tempDirArg} ]; then echo "Temp directory: ${tempDirArg}"; fi

    h3_preprocessor.sh \
      -i ${occurrences}'/*' \
      -s ${specieskey} \
      -o ${specieskey}.parquet \
      -r ${params.outlier_h3_resolution} \
      -t ${task.cpus} \
      "${memoryArg}" "${tempDirArg}" \
      "${basisOfRecordArg}"

    echo "..Done"
    """
}


// Detect spatial outliers (per species) using DBSCAN
process dbscan {
    // cpus 2

    tag "${specieskey}"

    input:
      tuple val(specieskey), path(h3_binned_csv)

    output:
      tuple val(specieskey), path("${specieskey}_DBSCAN-outliers.txt.gz"), emit: dbscan_scores

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
      > "${specieskey}_DBSCAN-outliers.txt.gz"

    rm "${specieskey}_DBSCAN-scores.txt.gz"

    echo "..Done"
    """
}


// Workflow
workflow {

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
  
  // Pool species lists from different taxonomic groups
  pool_species_lists(ch_occcounts)

  // Channel with species keys for spatial outlier removal
  // NB. results returned by `splitText` operator are always terminated by a `\n` newline character, so we need to trim it
  ch_large_specieskey = pool_species_lists.out.occ_large.splitText().map{it -> it.trim()}

  // Prepare data for spatial outlier removal
  prepare_species(ch_occurrence_dir, ch_large_specieskey)



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

