#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


nuclear_chroms = "${params.genome}.nuclear.txt"
chrom_sizes = "${params.genome}.chrom_sizes.bed"
mappable = "${params.genome}.K76.mappable_only.bed"
centers = "${params.genome}.K76.center_sites.n100.nuclear.starch"

process call_hotspots {
	tag "${id}"

	// only publish varw_peaks and hotspots
	publishDir "${params.outdir}/hotspots"

    container "fwip/hotspot2:latest" // FIXME, add info to nextflow.config

	scratch true

	input:
	    tuple val(id), val(bam_file)

	output:
	    tuple val(id), val(bam_file), file(name)
	    tuple val(id), file("${id}.hotspots.fdr005.starch"), file("${id}.hotspots.fdr001.starch")

	script:
    name = "${id}.varw_peaks.fdr0.001.starch"
	"""
	samtools view -H ${bam_file} > header.txt
	cat ${nuclear_chroms} \
	| xargs samtools view -b ${bam_file} \
	| samtools reheader header.txt - \
	> nuclear.bam

	hotspot2.sh -F 0.05 -f 0.05 -p varWidth_20_${id} \
		-M ${mappable} \
    	-c ${chrom_sizes} \
    	-C ${centers} \
    	nuclear.bam \
    	peaks

	cd peaks
	hsmerge.sh -f 0.001 nuclear.allcalls.starch nuclear.hotspots.fdr0.001.starch
	rm -f nuclear.varw_peaks.*
	
    density-peaks.bash \
		./ \
		varWidth_20_${indiv_id}_${cell_type} \
		nuclear.cutcounts.starch \
		nuclear.hotspots.fdr0.001.starch \
		${chrom_sizes} \
		nuclear.varw_density.fdr0.001.starch \
		nuclear.varw_peaks.fdr0.001.starch \
		\$(cat nuclear.cleavage.total)
		
        cp nuclear.varw_peaks.fdr0.001.starch ../${name}
    	cp nuclear.hotspots.fdr0.05.starch ../${id}.hotspots.fdr0.05.starch
    	cp nuclear.hotspots.fdr0.001.starch ../${id}.hotspots.fdr0.001.starch
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