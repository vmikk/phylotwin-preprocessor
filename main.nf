#!/usr/bin/env nextflow

params.occurrences = './gbif_occurrences/'
params.outdir      = './gbif_occurrences_processed/'

// Count number of occurrences per species
process count_occurrences {

    input:
      path occurrences

    output:
      path "Occurrences_counts.parquet", emit: occ_counts_parquet
      path "Occurrences_counts.csv",     emit: occ_counts_csv

    script:
    """
    echo "Counting species occurrences"
    echo "Input directory: " ${input}

    extra_args = ""

    if [ -n ${task.memory} ]; then
        extra_args += " -m ${task.memory}"
    fi
    if [ -n ${task.tempDir} ]; then
        extra_args += " -x ${task.tempDir}"
    fi

    gbif_records_count.sh \
      -i ${occurrences} \
      -o "Occurrences_counts.parquet" \
      -t ${task.cpus} \
      "$extra_args"

    echo "..Done"
    """
}


// Define species set (given phylogenetic trees available)
process define_species_set {

    tag "${phylotree}"

    input:
      path occ_counts_csv
      tuple val(phylotree), path(phylotree_processed), path(phylotree_initial)

    output:
      path "${phylotree}_Occurrences_large.txt",  emit: occ_large
      path "${phylotree}_Occurrences_small.txt",  emit: occ_small
      path "${phylotree}_stats.txt",              emit: stats

    script:
    """
    echo -e "Defining species set\n"
    echo "Occurrences counts: "   ${occ_counts_csv}
    echo "Phylogenetic tree: "    ${phylotree}
    echo "Processed tree file: "  ${phylotree_processed}
    echo "Initial tree file: "    ${phylotree_initial}
    echo "Occurrence threshold: " ${params.outlier_occurrence_threshold}

    echo -e "\nCounting the number of GBIF records represented and not represented in a phylogenetic tree\n"

    gbif_records_in_tree.R \
      --gbif                 ${occ_counts_csv} \
      --tree                 ${phylotree_processed} \
      --sourcetree           ${phylotree_initial} \
      --occurrence_threshold ${params.outlier_occurrence_threshold} \
      --outputstats          ${phylotree}_stats.txt \
      --output_large_occ     ${phylotree}_Occurrences_large.txt \
      --output_low_occ       ${phylotree}_Occurrences_small.txt

    """
}


// Workflow
workflow {

  // Input directory with GBIF species occurrences (parquet files)
  ch_occurrence_dir = Channel.fromPath(params.occurrences)

  // Count species occurrences
  count_occurrences(ch_occurrence_dir)


}



