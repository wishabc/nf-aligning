#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { set_key_for_group_tuple } from "./helpers"

def split_fasta_file(file_path):
    file(file_path).splitFasta(by: params.chunk_size, file: true)


def remove_ambiguous_bases(adapter) {
    def x = adapter.takeWhile { it in ['A', 'C', 'T', 'G'] }
    if (x != adapter) {
        println("WARN: Adapter '${adapter}' contains ambiguous bases, using '${x}' instead")
    }
    return x
}


def parse_legacy_adapter_file(adapter_file) {
  // returns two values, p7 and p5
  adapters = [:]
  adapter_file.readLines().each {
    columns = it.split()
    adapters[columns[0]] = remove_ambiguous_bases(columns[1])
  }
  return [adapters.P7, adapters.P5]
}

process fastp_adapter_trim {

    cpus 3

    input:
        tuple val(sample_id), path(r1), path(r2), val(is_paired)
        tuple val(adapterP7), val(adapterP5)

    output:
        tuple val(sample_id), path(name1), path(name2), val(is_paired), emit: fastq
        tuple val(meta), path('fastp.json'), emit: metrics_json
        tuple val(meta), path('fastp.html'), emit: metrics_html

    script:
    name1 = 'out.r1.fastq.gz'
    if is_paired
        name2 = 'out.r2.fastq.gz'
        """
        fastp \
            --in1 "${r1}" \
            --in2 "${r2}" \
            --adapter_sequence    "${adapterP7}" \
            --adapter_sequence_r2 "${adapterP5}" \
            --out1 "${name1}" \
            --out2 "${name2}" \
            --disable_quality_filtering \
            --disable_length_filtering \
            --thread ${task.cpus}
        """
    else
        name2 = './'
        """
        fastp \
        -i "${r1}" \
        --adapter_sequence "${adapterP7}" \
        -o ${name} \
        --disable_quality_filtering \
        --disable_length_filtering \
        --thread ${task.cpus}
        """
}

workflow trimReads {
    take: // [sample_id, r1, r2, is_paired]
        data
    main:
        fasta_chunks = data.map(it ->
            tuple(it[0], 
                  split_fasta_file(it[1]),
                  it[3] ? split_fasta_file(it[2]) : './',
                  it[3]
                )
            ).transpose()
        legacy_adapters = parse_legacy_adapter_file(params.adapter_file)
        fasta_files = fastp_adapter_trim(fasta_chunks, legacy_adapters)
    emit:
        fasta_files
}

workflow trimReadsFromFile {
    reads_ch = Channel.fromPath(params.samples_file)
        .splitCsv(header:true, sep:'\t')
		.map(row -> tuple( row.sample_id, row.reads1, row.reads2, row.is_paired))
    trimReads(reads_ch)
}

workflow {
    trimReadsFromFile()
}