#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { set_key_for_group_tuple } from "./aligning"

def remove_ambiguous_bases(adapter) {
    x = adapter.takeWhile { it in ['A', 'C', 'T', 'G'] }
    if (x != adapter) {
        println("WARN: Adapter '${adapter}' contains ambiguous bases, using '${x}' instead")
    }
    return x
}

process fastp_adapter_trim {
    cpus params.threads
    scratch true
    container "${params.container}"
    tag "${align_id}"
    publishDir "${params.outdir}/${sample_id}/stats/${align_id}", pattern: "fastp*"

    input:
        tuple val(sample_id), val(align_id), path(r1), path(r2), val(adapterP7), val(adapterP5), val(is_paired)

    output:
        tuple val(sample_id), path(name1), path(name2), val(is_paired), emit: fastq
        tuple val(sample_id), path('fastp.json'), emit: json
        tuple val(sample_id), path('fastp.html'), emit: html

    script:
    name1 = "${align_id}.1.trimmed.fastq.gz"
    if (adapterP7 == "") {
        command = is_paired ? "--detect_adapter_for_pe" : ""
    } else {
        command = "--adapter_sequence ${adapterP7}" + (is_paired ? " --adapter_sequence_r2 ${adapterP5}" : "")
    }
    
    if (is_paired) {
        name2 = "${align_id}.2.trimmed.fastq.gz"
        """
        fastp --in1 "${r1}" \
            --in2 "${r2}" \
            ${command} \
            --out1 "${name1}" \
            --out2 "${name2}" \
            --disable_quality_filtering \
            --disable_length_filtering \
            --thread ${task.cpus}
        """
    } else {
        name2 = './'
        """
        fastp -i "${r1}" \
            ${command} \
            -o ${name1} \
            --disable_quality_filtering \
            --disable_length_filtering \
            --thread ${task.cpus}
        """
    }
}

workflow trimReads {
    take: // [sample_id, align_id, r1, r2, adapter7, adapter5, is_paired]
        data
    main:
        trimmed = fastp_adapter_trim(data).fastq
    emit:
        trimmed
}

workflow trimReadsFromFile {
    main:
        Channel.fromPath(params.samples_file)
            | splitCsv(header:true, sep:'\t')
            | map(
                row -> tuple(
                    row.sample_id,
                    row.align_id,
                    row.reads1,
                    row.type == 'paired' ? row.reads2 : file('./'),
                    remove_ambiguous_bases(row.adapterP7),
                    row.type == 'paired' ? remove_ambiguous_bases(row.adapterP5) : "",
                    row.type == 'paired'
                )
            )
            | set_key_for_group_tuple
            | trimReads
    emit:
        trimReads.out
}

workflow {
    trimReadsFromFile()
}