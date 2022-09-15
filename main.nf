#!/usr/bin/env nextflow
include { alignReads; get_container; set_key_for_group_tuple } from "./aligning"
include { callHotspots } from "./hotspots_calling"
include { trimReadsFromFile; trimReads } from "./trimming"


process symlink_or_download {
    publishDir "${params.outdir}/${sample_id}/stats", pattern: "${metadata}"
    publishDir "${params.outdir}/${sample_id}", pattern: "${srr}/*.fastq.gz"
    cpus params.threads
    tag "${srr}"
    container "${params.container}"

    input:
        tuple val(sample_id), val(srr)
    output:
        tuple val(sample_id), val(srr), path("${srr}/${srr}_1.fastq.gz"), path("${srr}/${srr}_2.fastq.gz"), path("${srr}/${srr}.fastq.gz"), emit: fastq
        tuple val(srr), path(metadata), emit: meta

    script:
    metadata = "${srr}_info.json"
    """
    prefetch -L 1 ${srr}
    ffq -o ${metadata} ${srr} 2>&1 || echo 'No metadata downloaded.' > ${metadata}
    fasterq-dump -L 1 -f --threads ${task.cpus} -O ${srr} ${srr} 2>&1
    find ./${srr} -name "*.fastq" -exec pigz {} \\;
    """

}

workflow downloadFiles {
    main:
        ids_channel = set_key_for_group_tuple(Channel.fromPath(params.samples_file)
            .splitCsv(header:true, sep:'\t')
            .map(row -> tuple(row.sample_id, row.align_id)))
        reads = symlink_or_download(ids_channel).fastq
        // Check if is_paired and convert to trimming pipeline format
        output = reads.map( 
            it -> (!file(it[2]).exists() || file(it[2]).size() == 0) ?
              tuple(it[0], it[1], it[4], path("./"), "", "", false) :
              tuple(it[0], it[1], it[2], it[3], "", "", true)
        )
    emit:
        output
}

workflow alignTrimmed {
    take:
        data
    main:
        alignReads(data) | callHotspots
    emit:
        callHotspots.out
}

workflow alignFromSRA {
    downloadFiles() | trimReads | alignTrimmed
}

workflow alignOnly {
    trimReadsFromFile() | alignReads
}

workflow {
    trimReadsFromFile() | alignTrimmed
}
