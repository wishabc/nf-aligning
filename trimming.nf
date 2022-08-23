#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { set_key_for_group_tuple } from "./helpers"

def split_fasta_file(file_path):
    file(file_path).splitFasta(by: params.chunk_size, file: true)


def parse_legacy_adapter_file(adapter_file) {
  // returns two values, p7 and p5
  adapters = [:]
  adapter_file.readLines().each {
    columns = it.split()
    adapters[columns[0]] = columns[1]
  }
  return [adapters.P7, adapters.P5]
}


process fastp_adapter_trim {

  // TODO: Get a newer version in module system!
  // TODO: Set this module in the RNA pipeline
  //module 'fastp/0.21.0'
  cpus 3

  input:
    tuple val(meta), path(r1), path(r2), val(is_paired)

  output:
    tuple val(sample_id), path('out.r1.fastq.gz'), path('out.r2.fastq.gz'), val(is_paired), emit: fastq
    tuple val(meta), path('fastp.json'),      emit: metrics_json
    tuple val(meta), path('fastp.html'),      emit: metrics_html

  script:
    // TODO: Double-check adapter ordering
    a7 = remove_ambiguous_bases(adapterP7)
    a5 = remove_ambiguous_bases(adapterP5)

    """
    fastp \
      --in1 "${r1}" \
      --in2 "${r2}" \
      --adapter_sequence    "${a7}" \
      --adapter_sequence_r2 "${a5}" \
      --out1 "out.r1.fastq.gz" \
      --out2 "out.r2.fastq.gz" \
      --disable_quality_filtering \
      --disable_length_filtering \
      --thread ${task.cpus}
    """
}

workflow splitFasta {
    take: // [sample_id, r1, r2, is_paired]
        data
    main:
        fasta_chunks = data.map(it ->
            tuple(it[0], 
                  split_fasta_file(it[1]),
                  it[3] ? split_fasta_file(it[2]) : it[2],
                  it[3]
                )
            ).transpose()
    emit:
        fasta_chunks
}

workflow trimReads {
    take:
        reads
    main:
        // FIXME
        fasta_files = splitFasta(reads).map { it ->
            (p7, p5) = parse_legacy_adapter_file(params.adapter_file);
            return [*it, p7, p5]
        }
        //TODO split into two channels for single and paired end reads
        fastp_adapter_trim(fasta_files)

        fastp_adapter_trim.out.fastq
        .map { meta, reads -> [ meta, meta.trim_to_length, reads[0], reads[1] ] }
        | trim_to_length
        | map { meta, r1, r2 -> [meta, meta.umi, r1, r2] }
        | add_umi_info
    emit:
        add_umi_info.out
}