#!/usr/bin/env nextflow
include { alignReads } from "./aligning"
include { call_hotspots } from "./hotspots_calling"
include { trimReadsFromFile; trimReads } from "./trimming"
include { get_container; set_key_for_group_tuple } from "./helpers"


process symlink_or_download {
    publishDir "${params.outdir}/${sample_id}/stats", pattern: "${metadata}"
    publishDir "${params.outdir}/${sample_id}", pattern: "${srr}/*.fastq.gz"
    cpus params.threads
    tag "${srr}"
    maxForks 4
    scratch true
    errorStrategy 'ignore'
    container "${params.container}"
    containerOptions "--network=host"

    input:
        tuple val(sample_id), val(srr)
    output:
        tuple val(sample_id), val(srr), path("${srr}/${srr}_1.fastq.gz"), path("${srr}/${srr}_2.fastq.gz"), path("${srr}/${srr}.fastq.gz"), emit: fastq
        tuple val(srr), path(metadata), emit: meta

    script:
    metadata = "${srr}_info.json"
    """
    if prefetch -L 1 ${srr}; then
        echo 'Prefetched'
        ffq -o ${metadata} ${srr} 2>&1 || echo 'No metadata downloaded.' > ${metadata}
        fasterq-dump -L 1 -f --threads ${task.cpus} -O ${srr} ${srr} 2>&1
        find ./${srr} -name "*.fastq" -exec pigz {} \\;
    else
        touch ${metadata}
    fi
    if ! [ -d ${srr} ]; then
        mkdir ${srr}
    fi
    touch ${srr}/${srr}_1.fastq.gz
    touch ${srr}/${srr}_2.fastq.gz
    touch ${srr}/${srr}.fastq.gz
    """

}

workflow downloadFiles {
    main:
        ids_channel = set_key_for_group_tuple(Channel.fromPath(params.samples_file)
            .splitCsv(header:true, sep:'\t')
            .map(row -> tuple(row.sample_id, row.align_id))).unique { it[1] }
        reads = symlink_or_download(ids_channel).fastq
        // Check if is_paired and convert to trimming pipeline format
        output = reads.filter { file(it[2]).size() + file(it[3]).size() + file(it[4]).size() != 0 }
        .map( 
            it -> (file(it[2]).size() == 0) ?
              tuple(it[0], it[1], it[4], "${projectDir}", "", "", false) :
              tuple(it[0], it[1], it[2], it[3], "", "", true)
        )
    emit:
        output
}

workflow alignTrimmed {
    take:
        data
    main:
        out = data | alignReads | call_hotspots
    emit:
        out.peaks
}

workflow alignFromSRA {
    downloadFiles() | trimReads | alignTrimmed
}

workflow {
    trimReadsFromFile() | alignTrimmed
}
