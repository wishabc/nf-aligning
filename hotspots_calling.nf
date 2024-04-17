#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { get_container } from "./aligning"

// Need to contenirize at some point
process call_hotspots {
	tag "${id}"
	label 'high_mem'
	publishDir "${params.outdir}/${id}", pattern: "${name}"
	publishDir "${params.outdir}/${id}", pattern: "nuclear.SPOT.txt", saveAs: { "${id}.SPOT.txt" } 
	publishDir "${params.outdir}/${id}", pattern: "nuclear.cleavage.total", saveAs: { "${id}.cleavage.total" }
	publishDir "${params.outdir}/${id}", pattern: "nuclear.density.bw", saveAs: { "${id}.density.bw" }
	publishDir "${params.outdir}/${id}", pattern: "nuclear.hotspot2.info", saveAs: { "${id}.hotspot2.info" }
	publishDir "${params.outdir}/${id}", pattern: "nuclear.cutcounts.starch", saveAs: { "${id}.cutcounts.starch" }
    // publishDir "${params.outdir}/${id}", pattern: "nuclear.hotspots.fdr0.05.starch", saveAs: { "${id}.cutcounts.starch" }

	//container "${params.container}"
	//containerOptions "${get_container(params.nuclear_chroms)} ${get_container(params.chrom_sizes_bed)} ${get_container(params.mappable)} ${get_container(params.centers)}"
	module "modwt/1.0:kentutil/302:bedops/2.4.35-typical:bedtools/2.25.0:samtools/1.3:hotspot2/2.1.1"

	input:
	    tuple val(id), path(bam_file), path(bam_file_index)

	output:
	    tuple val(id), path(name), path(spot), path("nuclear.cleavage.total"), path("nuclear.density.bw"), path("nuclear.hotspot2.info"), path("nuclear.cutcounts.starch")

	script:
	name = "${id}.peaks.fdr0.001.starch"
	spot = "nuclear.SPOT.txt"
	renamed_input = "nuclear.${bam_file.extension}"
	"""
	export TMPDIR=\$(mktemp -d)

	# workaround for hotspots2 naming scheme
	ln -s ${bam_file} ${renamed_input}
	ln -s ${bam_file_index} ${renamed_input}.${bam_file_index.extension}

	hotspot2.sh -F 0.001 -f 0.001 \
		-p "varWidth_20_${id}" \
		-M "${params.mappable}" \
		-c "${params.chrom_sizes_bed}" \
		-C "${params.centers}" \
		${renamed_input} \
		'.'

	mv nuclear.peaks.starch ${name}
	echo -e "hotspot2-num-bases\t\$(unstarch --bases    "${name}")" >> nuclear.hotspot2.info
  	echo -e "hotspot2-num-spots\t\$(unstarch --elements "${name}")" >> nuclear.hotspot2.info
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
		| map(row -> tuple(
			row.sample_id, 
			file("${params.basepath}/${row.sample_id}/${row.sample_id}.filtered.cram"),
			file("${params.basepath}/${row.sample_id}/${row.sample_id}.filtered.cram.crai"),
			))
		| filter { it[1].exists() }
	callHotspots(metadata)
}
workflow {
	Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
		| map(row -> tuple( row.ag_id, row.bam_file, "${row.bam_file}.crai"))
        | callHotspots
}
