#!/usr/bin/env nextflow
include { alignReads } from "./aligning"
include { callHotspots } from "./hotspots_calling"
include { trimReadsFromFile; trimReads } from "./trimming"



process symlink_or_download {
    publishDir "${outdir}/fasta", pattern: "${metadata}"
    cpus params.threads
    input:
        tuple val(sample_id), val(srr)
    output:
        tuple val(sample_id), val(srr), path("${srr}_1.fastq.gz"), path("${srr}_2.fastq.gz"), path("${srr}.fastq.gz"), emit: fastq
        tuple val(srr), path(metadata), emit: meta

    script:
    metadata = "${srr}_info.json"
    """
    if [ ! -d ${params.readdirectory}/${align_id}/ ] || test -n "\$(find . -maxdepth 1 -wholename "${params.readdirectory}/${align_id}/*.fastq.gz" -print -quit)"
    then
        echo "${params.readdirectory} does not contain expected FastQ files. Downloading"
        prefetch ${srr}
        cd ${srr}
        ffq -o ${metadata} ${srr} || echo 'No metadata downloaded.' > no_metadata.json
        fasterq-dump --threads ${task.cpus} ${srr}
        find . -name "*.fastq" -exec pigz {} \\;
    else
        mkdir ${srr}
        ln -s ${params.readdirectory}/${srr}/* ${srr}
    fi
    """

}

workflow downloadFiles {
    main:
        ids_channel = Channel.fromPath(params.samples_file)
            .splitCsv(header:true, sep:'\t')
            .map(row -> tuple(row.sample_id, row.align_id))
        reads = symlink_or_download(ids_channel).fastq
        // Check if is_paired and convert to trimming pipeline format
        output = reads.map( 
            it -> (!it[2].exists() || it[2].length() == 0) ?
              tuple(it[0], it[1], it[4], path("./"), "", "", false) :
              tuple(it[0], it[1], it[2], it[3], "", "", true)
        )
    emit:
        output
}

workflow alignTrimmed {
    alignReads() | callHotspots
}

workflow alignFromSRA {
    downloadFiles() | trimReads | alignTrimmed
}

workflow {
    trimReadsFromFile() | alignTrimmed
}
