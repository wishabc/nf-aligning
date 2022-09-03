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

// process split_fasta_paired {
//     container "${params.container}"
//     scratch true
//     input:
//         tuple val(sample_id), val(align_id), path(fastq1), path(fastq2)
//     output:
//         tuple val(sample_id), path("out/${name_prefix1}*.fastq.gz"), path("out/${name_prefix2}*.fastq.gz")
//     script:
//     name_prefix1 = "${align_id}.1."
//     name_prefix2 = "${align_id}.2."
//     """
//     mkdir out
//     zcat "${fastq1}" \
//     | split -l "${params.chunk_size}" \
//       --filter='gzip -1 > out/\$FILE.fastq.gz' \
//       - "${name_prefix1}"
//     zcat "${fastq2}" \
//     | split -l "${params.chunk_size}" \
//       --filter='gzip -1 > out/\$FILE.fastq.gz' \
//       - "${name_prefix2}"
//     """

// }

// process split_fasta_single {
//     container "${params.container}"
//     scratch true
//     input:
//         tuple val(sample_id), val(align_id), path(fastq)
//     output:
//         tuple val(sample_id), path("out/${name_prefix}*")
//     script:
//     name_prefix = "${align_id}."
//     """
//     mkdir out
//     zcat "${fastq}" \
//     | split -l "${params.chunk_size}" \
//       --filter='gzip -1 > out/\$FILE.fastq.gz' \
//       - "${name_prefix}"
//     """

// }

process fastp_adapter_trim {
    cpus params.threads
    scratch true
    container "${params.container}"
    publishDir "${params.outdir}/${sample_id}/stats/"

    input:
        tuple val(sample_id), val(align_id), path(r1), path(r2), val(adapterP7), val(adapterP5), val(is_paired)

    output:
        tuple val(sample_id), path(name1), path(name2), val(is_paired), emit: fastq
        // tuple val(sample_id), path('fastp.json'), emit: json
        // tuple val(sample_id), path('fastp.html'), emit: html

    script:
    name1 = "${align_id}.1.trimmed.fastq.gz"
    if (is_paired) {
        name2 = "${align_id}.2.trimmed.fastq.gz"
        """
        fastp --in1 "${r1}" \
            --in2 "${r2}" \
            --adapter_sequence    "${adapterP7}" \
            --adapter_sequence_r2 "${adapterP5}" \
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
            --adapter_sequence "${adapterP7}" \
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
        // fasta_chunks = data.branch{ 
        //     paired: it[6]
        //     single: true 
        // }
        // split_single = split_fasta_single(
        //     fasta_chunks.single.map(it -> tuple(it[0], it[1], it[2]))
        // ).join(
        //     fasta_chunks.single.map(it -> tuple(it[0], file('./'),
        //         remove_ambiguous_bases(it[4]),
        //         it[5],
        //         it[6])
        //     )
        // ).transpose()
        // split_paired = split_fasta_paired(
        //     fasta_chunks.paired.map(it -> tuple(it[0], it[1], it[2], it[3]))
        // ).map(
        //     it -> tuple(it[0], it[1].sort(), it[2].sort())
        // ).join(
        //         fasta_chunks.paired.map(it -> tuple(it[0],
        //         remove_ambiguous_bases(it[4]), 
        //         remove_ambiguous_bases(it[5]), 
        //         it[6])
        //     )
        // ).transpose()

        // trimmed = fastp_adapter_trim(split_single.mix(split_paired)).fastq
        trimmed = fastp_adapter_trim(data).fastq
    emit:
        trimmed
}

workflow trimReadsFromFile {
    main:
        reads_ch = Channel.fromPath(params.samples_file)
            .splitCsv(header:true, sep:'\t')
            .map(row -> tuple(row.sample_id, row.align_id, row.reads1,
                                row.reads2 != '' ? row.reads2 : file('./'), row.adapterP5, 
                                row.adapterP7, row.type == 'paired'))
        trimReads(set_key_for_group_tuple(reads_ch))
    emit:
        trimReads.out
}

workflow {
    trimReadsFromFile()
}