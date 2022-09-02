#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

def split_fasta_file(file_path) {
    chunks = Channel.fromPath(file_path).splitFasta(by: params.chunk_size, file: true)
    chunks.view()
    return chunks
}

def remove_ambiguous_bases(adapter) {
    x = adapter.takeWhile { it in ['A', 'C', 'T', 'G'] }
    if (x != adapter) {
        println("WARN: Adapter '${adapter}' contains ambiguous bases, using '${x}' instead")
    }
    return x
}

process fastp_adapter_trim {
    cpus params.threads
    container "${params.container}"
    publishDir "${params.outdir}/${sample_id}/stats/${simple_name}", pattern: "fastp*"
    //scratch true

    input:
        tuple val(sample_id), val(r1), val(r2), val(adapterP7), val(adapterP5), val(is_paired)

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
        fasta_chunks = data.map{ it ->
            tuple(it[0], 
                  split_fasta_file(it[1]),
                  it[5] ? split_fasta_file(it[2]) : './',
                   remove_ambiguous_bases(it[3]),
                   it[5] ? remove_ambiguous_bases(it[4]) : it[4],
                  it[5]
                )
        }.transpose()
        fastp_adapter_trim(fasta_chunks)
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