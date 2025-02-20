#!/usr/bin/env nextflow
include { alignReads } from "./aligning"
include { call_hotspots } from "./hotspots_calling"
include { trimReads } from "./trimming"
include { get_container; set_key_for_group_tuple } from "./helpers"


def remove_ambiguous_bases(adapter) {
    def x = adapter ? adapter.takeWhile { it in ['A', 'C', 'T', 'G'] } : ""
    if (x && (x != adapter)) {
        println("WARN: Adapter '${adapter}' contains ambiguous bases, using '${x}' instead")
    }
    return x
}


process symlink_or_download {
    publishDir "${params.outdir}/${ag_id}/stats", pattern: "${metadata}"
    publishDir "${params.outdir}/${ag_id}", pattern: "${srr}/*.fastq.gz"
    cpus params.threads
    tag "${srr}"
    maxForks 4
    scratch true
    errorStrategy 'ignore'
    container "${params.container}"
    containerOptions "--network=host"

    input:
        tuple val(ag_id), val(srr)
    output:
        tuple val(ag_id), val(srr), path("${srr}/${srr}_1.fastq.gz"), path("${srr}/${srr}_2.fastq.gz"), path("${srr}/${srr}.fastq.gz"), emit: fastq
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


workflow alignFromSRA {
    ids_channel = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.ag_id, row.align_id))
        | unique { it[1] }
        | set_key_for_group_tuple


    output = symlink_or_download(ids_channel).fastq
        | filter { file(it[2]).size() + file(it[3]).size() + file(it[4]).size() != 0 }
        | map( 
            it -> (file(it[2]).size() == 0) ?
                tuple(it[0], it[1], it[4], "${projectDir}", "", "", false) :
                    tuple(it[0], it[1], it[2], it[3], "", "", true)
            ) // Check if is_paired and convert to trimming pipeline format
        | trimReads
        | alignReads 
        | call_hotspots
}


workflow {
    Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(
            row -> tuple(
                row.ag_id,
                row.align_id,
                row.reads1,
                row.type == 'paired' ? row.reads2 : "${row.reads1}.mock",
                remove_ambiguous_bases(row?.adapterP7),
                row.type == 'paired' ? remove_ambiguous_bases(row?.adapterP5) : "",
                row.type == 'paired'
            )
        )
        | unique { it[1] }
        | set_key_for_group_tuple
        | trimReads
        | alignReads 
        //| call_hotspots
}
