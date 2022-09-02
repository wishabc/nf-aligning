#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


def remove_ambiguous_bases(adapter) {
    x = adapter.takeWhile { it in ['A', 'C', 'T', 'G'] }
    if (x != adapter) {
        println("WARN: Adapter '${adapter}' contains ambiguous bases, using '${x}' instead")
    }
    return x
}

process split_fasta_paired {
    container "${params.container}"
    scratch true
    input:
        tuple val(sample_id), path(fastq1), path(fastq2)
    output:
        tuple val(sample_id), path("out/{${name_prefix1}}*.fastq.gz"), path("out/{${name_prefix2}}*.fastq.gz")
    script:
    name_prefix1 = "${fastq1.baseName}."
    name_prefix2 = "${fastq2.baseName}."
    """
    mkdir out
    zcat "${fastq1}" \
    | split -l "${params.chunk_size}" \
      --filter='gzip -1 > out/\$FILE.fastq.gz' \
      - "${name_prefix1}"
    zcat "${fastq2}" \
    | split -l "${params.chunk_size}" \
      --filter='gzip -1 > out/\$FILE.fastq.gz' \
      - "${name_prefix2}"
    """

}

process split_fasta_single {
    container "${params.container}"
    scratch true
    input:
        tuple val(sample_id), path(fastq)
    output:
        tuple val(sample_id), path("out/${name_prefix}*")
    script:
    name_prefix = "${fastq.baseName}."
    """
    mkdir out
    zcat "${fastq}" \
    | split -l "${params.chunk_size}" \
      --filter='gzip -1 > out/\$FILE.fastq.gz' \
      - "${name_prefix}"
    """

}

process fastp_adapter_trim {
    cpus params.threads
    scratch true
    container "${params.container}"
    publishDir "${params.outdir}/${sample_id}/stats/${simple_name}", pattern: "fastp*"

    input:
        tuple val(sample_id), path(r1), path(r2), val(adapterP7), val(adapterP5), val(is_paired)

    output:
        tuple val(sample_id), path(name1), path(name2), val(is_paired)
        //tuple val(sample_id), path('fastp.json')
        //tuple val(sample_id), path('fastp.html')

    script:
    simple_name = r1.simpleName
    name1 = "${r1.simpleName}.r1.trimmed.fastq.gz"
    if (is_paired) {
        name2 = "${r1.simpleName}.r2.trimmed.fastq.gz"
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
    take: // [sample_id, r1, r2, adapter7, adapter5, is_paired]
        data
    main:
        fasta_chunks = data.branch{ 
            paired: it[5]
            single: true 
        }
        split_single = split_fasta_single(
            fasta_chunks.single.map(it -> tuple(it[0], it[1]))
        ).join(
            fasta_chunks.single.map(it -> tuple(it[0], file('./'),
            remove_ambiguous_bases(it[3]),
            it[4], 
            it[5]))
        ).transpose()
        par = split_fasta_paired(
            fasta_chunks.paired.map(it -> tuple(it[0], it[1], it[2]))
        ).map(it -> tuple(it[0], it[1].sort(), it[2].sort())).join(
            fasta_chunks.paired.map(it -> tuple(it[0],
            remove_ambiguous_bases(it[3]), 
            remove_ambiguous_bases(it[4]), 
            it[5]))
        ).transpose()
        fastp_adapter_trim(split_single.mix(split_paired))
    emit:
        fastp_adapter_trim.out
}

workflow trimReadsFromFile {
    main:
        reads_ch = Channel.fromPath(params.samples_file)
            .splitCsv(header:true, sep:'\t')
            .map(row -> tuple(row.sample_id, row.reads1,
                                row.reads2, row.adapterP5, 
                                row.adapterP7, row.type == 'paired'))
        trimReads(reads_ch)
    emit:
        trimReads.out
}

workflow {
    trimReadsFromFile()
}