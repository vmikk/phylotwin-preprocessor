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



