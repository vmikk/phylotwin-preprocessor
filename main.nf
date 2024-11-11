#!/usr/bin/env nextflow

params.occurrences  = './gbif_occurrences/'
params.outdir       = './gbif_occurrences_processed/'

// Phylogenetic trees
params.t1   = './original_trees/'   // Original trees
params.t2   = './processed_trees/'  // Trees with tips matched to GBIF specieskeys

phylo_trees = [
    tuple("Plants",      "${params.t1}/Tietje_2023_Seed_plants_TACT.nwk", "${params.t2}/Tietje_2023_Seed_plants_TACT.nwk"),
    tuple("Mammals",     "${params.t1}/VertLife_Mammals.nwk",             "${params.t2}/VertLife_Mammals.nwk"),
    tuple("Butterflies", "${params.t1}/Kawahara_2023_Butterflies.nwk",    "${params.t2}/Kawahara_2023_Butterflies.nwk")
]


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
      "${memoryArg}" "${tempDirArg}"

    echo "..Done"
    """
}


// Define species set (given phylogenetic trees available)
process define_species_set {

    tag "${taxon}"

    input:
      path occ_counts_csv
      tuple val(taxon), path(phylotree_initial), path(phylotree_processed)

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

}


// Workflow
workflow {

  // Input directory with GBIF species occurrences (parquet files)
  ch_occurrence_dir = Channel.fromPath(params.occurrences)

  // Count species occurrences
  count_occurrences(ch_occurrence_dir)

  // Phylogenetic trees
  ch_phylotrees = Channel.fromList(phylo_trees)

  // Define species set, split into large and small occurrence sets
  define_species_set(
    count_occurrences.out.occ_counts_csv,
    ch_phylotrees
  )


// On completion
workflow.onComplete {
    println "Pipeline completed at : $workflow.complete"
    println "Duration              : ${workflow.duration}"
    println "Execution status      : ${workflow.success ? 'All done!' : 'Failed' }"
}


