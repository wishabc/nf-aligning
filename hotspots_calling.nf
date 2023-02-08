#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { get_container } from "./aligning"

process call_hotspots {
	tag "${id}"
	publishDir "${params.outdir}/${id}", pattern: "${name}"
    //container "${params.container}"
	//containerOptions "${get_container(params.nuclear_chroms)} ${get_container(params.chrom_sizes_bed)} ${get_container(params.mappable)} ${get_container(params.centers)}"
	//scratch true
	//errorStrategy 'ignore'
	module "hotspot2/2.1.1:kentutil/302:bedops/2.4.35-typical:bedtools/2.25.0:modwt/1.0"

	input:
	    tuple val(id), path(bam_file), path(bam_file_index)

	output:
	    tuple val(id), path(name)

	script:
	name = "${id}.peaks.fdr0.001.starch"
	"""
	samtools view -H ${bam_file} > header.txt
	cat ${params.nuclear_chroms} \
		| xargs samtools view -b ${bam_file} \
		| samtools reheader header.txt - \
		> nuclear.bam

	hotspot2.sh -F 0.001 -f 0.001 \
		-p "varWidth_20_${params.id}" \
		-M "${params.mappable}" \
		-c "${params.chrom_sizes_bed}" \
		-C "${params.centers}" \
		nuclear.bam \
		'.'

	mv ${params.id}.peaks.starch ${name}
	"""
}

workflow callHotspots {
	take:
		bam_files
	main:
		call_hotspots(bam_files)
	emit:
		call_hotspots.out
}

workflow hotspots {
	params.basepath = "/net/seq/data2/projects/sabramov/ENCODE4/atac_aligning/output"
	metadata = Channel.fromPath(params.samples_file)
		| splitCsv(header:true, sep:'\t')
		| map(row -> tuple( row.sample_id, 
			file("${params.basepath}/${row.sample_id}/${row.sample_id}.filtered.cram"),
			file("${params.basepath}/${row.sample_id}/${row.sample_id}.filtered.cram.crai"),
			))
		| filter { it[1].exists() }
	callHotspots(metadata)
}
workflow {
	metadata = Channel
      .fromPath(params.samples_file)
      .splitCsv(header:true, sep:'\t')
		  .map(row -> tuple( row.ag_id, row.bam_file, "${row.bam_file}.{crai,bai}"))
	callHotspots(metadata)
}