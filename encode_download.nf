#!/usr/bin/env nextflow

include { convert_to_cram } from "./aligning"
nextflow.enable.dsl = 2

params.conda = "$moduleDir/environment.yml"


process download_encode {
    conda params.conda
    maxForks 2
    tag "${encode_id}"
    scratch true

    input:
        tuple val(encode_id), val(download_path), val(md5)

    output:
        tuple val(encode_id), path(name), path("${name}.bai")

    script:
    name = "${encode_id}.bam"
    """
    wget ${download_path} -O ${name}
    if [[ `md5sum ${name} | awk '{print \$1}'` != "${md5}" ]]; then
        exit 1
    fi
    samtools index ${name}
    """
}


workflow downloadEncode {
    take:
        metadata
    main:
        download_encode(metadata) | convert_to_cram
    emit:
        convert_to_cram.out
}


workflow {
    params.encode_meta = "/home/sabramov/encode_chipseq.tsv"
    metadata = Channel
		.fromPath(params.encode_meta)
		.splitCsv(header:true, sep:'\t')
		.map(row -> tuple( row.id, row.link, row.md5 ))
    
    downloadEncode(metadata)
}