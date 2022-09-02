#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { nuclearChromsContainer } from "./aligning"

process call_hotspots {
	tag "${id}"
	publishDir "${params.outdir}/${id}"
    container "${params.container}"
	scratch true

	input:
	    tuple val(id), path(bam_file)

	output:
	    tuple val(id), path(name), path("${id}.hotspots.fdr005.starch"), path("${id}.hotspots.fdr001.starch")

	script:
    name = "${id}.varw_peaks.fdr0.001.starch"
	"""
	samtools view -H ${bam_file} > header.txt
	cat ${params.nuclear_chroms} \
		| xargs samtools view -b ${bam_file} \
		| samtools reheader header.txt - \
		> nuclear.bam

	hotspot2.sh -F 0.05 -f 0.05 -p varWidth_20_${id} \
		-M ${params.mappable} \
    	-c ${params.chrom_sizes} \
    	-C ${params.centers} \
    	nuclear.bam \
    	peaks

	cd peaks
	hsmerge.sh -f 0.001 nuclear.allcalls.starch nuclear.hotspots.fdr0.001.starch
	rm -f nuclear.varw_peaks.*
	
    density-peaks.bash \
		./ \
		varWidth_20_${id} \
		nuclear.cutcounts.starch \
		nuclear.hotspots.fdr0.001.starch \
		${params.chrom_sizes} \
		nuclear.varw_density.fdr0.001.starch \
		nuclear.varw_peaks.fdr0.001.starch \
		\$(cat nuclear.cleavage.total)
		
        mv nuclear.varw_peaks.fdr0.001.starch ../${name}
    	mv nuclear.hotspots.fdr0.05.starch ../${id}.hotspots.fdr0.05.starch
    	mv nuclear.hotspots.fdr0.001.starch ../${id}.hotspots.fdr0.001.starch
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

workflow {
	metadata = Channel
      .fromPath(params.samples_file)
      .splitCsv(header:true, sep:'\t')
		  .map(row -> tuple( row.ag_id, row.bam_file))
	callHotspots(metadata)
}