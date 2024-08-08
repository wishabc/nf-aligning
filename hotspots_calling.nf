#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { get_container } from "./aligning"

// Need to contenirize at some point
process call_hotspots {
	tag "${id}"
	label 'high_mem'
	publishDir "${params.outdir}/${id}", pattern: "${id}*"
	publishDir "${params.outdir}/${id}", pattern: "nuclear.SPOT.txt", saveAs: { "${id}.SPOT.txt" } 
	publishDir "${params.outdir}/${id}", pattern: "nuclear.cleavage.total", saveAs: { "${id}.cleavage.total" }
	publishDir "${params.outdir}/${id}", pattern: "nuclear.density.bw", saveAs: { "${id}.density.bw" }
	publishDir "${params.outdir}/${id}", pattern: "nuclear.hotspot2.info", saveAs: { "${id}.hotspot2.info" }
	publishDir "${params.outdir}/${id}", pattern: "nuclear.cutcounts.starch", saveAs: { "${id}.cutcounts.starch" }
    publishDir "${params.outdir}/${id}", pattern: "nuclear.allcalls.starch", saveAs: { "${id}.allcalls.starch" }

	//container "${params.container}"
	//containerOptions "${get_container(params.nuclear_chroms)} ${get_container(params.chrom_sizes_bed)} ${get_container(params.mappable)} ${get_container(params.centers)}"
	module "modwt/1.0:kentutil/302:bedops/2.4.35-typical:bedtools/2.25.0:samtools/1.3:hotspot2/2.1.1"

	input:
	    tuple val(id), path(bam_file), path(bam_file_index)

	output:
	    tuple val(id), path(spot), path("nuclear.cleavage.total"), path("nuclear.density.bw"), path("nuclear.hotspot2.info"), path("nuclear.cutcounts.starch"), path("nuclear.allcalls.starch"), emit: out
        tuple val(id), path("${id}.hotspots.fdr0.001.starch"), path(name), emit: point1pr
        tuple val(id), path(hotspot1pr), path(peaks1pr), emit: 1pr
        tuple val(id), path(hotspot5pr), path(peaks5pr), emit: 5pr

	script:
    hotspot5pr = "${id}.hotspots.fdr0.05.starch"
    peaks5pr = "${id}.peaks.fdr0.05.starch"
    hotstpot1pr = "${id}.hotspots.fdr0.01.starch"
    peaks1pr = "${id}.peaks.fdr0.01.starch"
 
	name = "${id}.peaks.fdr0.001.starch"
	spot = "nuclear.SPOT.txt"
	renamed_input = "nuclear.bam"
	"""
    mkdir tmp # hotspot2 needs a tmp directory, scratch doesn't work due to limited tmp space
	export TMPDIR="\$PWD/tmp"

	# workaround for hotspots2 naming scheme
	ln -sf ${bam_file} ${renamed_input}
	ln -sf ${bam_file_index} ${renamed_input}.bai

	hotspot2.sh -F 0.05 -f 0.05 \
		-p "varWidth_20_${id}" \
		-M "${params.mappable}" \
		-c "${params.chrom_sizes_bed}" \
		-C "${params.centers}" \
		${renamed_input} \
		'.'

    hsmerge.sh -f 0.001 -m 50 nuclear.allcalls.starch ${id}.hotspots.fdr0.001.starch
    hsmerge.sh -f 0.01 -m 50 nuclear.allcalls.starch ${peaks1pr}

    bash density-peaks.bash \
        \$TMPDIR \
        "varWidth_20_${id}" \
        nuclear.cutcounts.starch \
        ${peaks1pr} \
        ${params.chrom_sizes_bed} \
        \$TMPDIR/nuclear.density.0.05.starch \
        ${hotspot1pr} \
        `cat nuclear.cleavage.total`

    bash density-peaks.bash \
        \$TMPDIR \
        "varWidth_20_${id}" \
        nuclear.cutcounts.starch \
        ${id}.hotspots.fdr0.001.starch \
        ${params.chrom_sizes_bed} \
        \$TMPDIR/nuclear.density.0.05.starch \
        ${name} \
        `cat nuclear.cleavage.total`

    mv nuclear.hotspots.starch ${hotspot5pr}
    mv nuclear.peaks.starch ${peaks5pr}

	echo -e "hotspot2-num-bases\t\$(unstarch --bases ${name})" >> nuclear.hotspot2.info
  	echo -e "hotspot2-num-spots\t\$(unstarch --elements ${name})" >> nuclear.hotspot2.info

    rm -r tmp
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
	Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
		| map(row -> tuple( row.ag_id, row.cram_file, row.cram_index ?: "${row.cram_file}.crai"))
        | callHotspots
}

workflow tmp {
    Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
		| map(row -> tuple( row.ag_id, row.cram_file, row.cram_index ?: "${row.cram_file}.crai", file("${params.outdir}/${row.ag_id}/${row.ag_id}.peaks.fdr0.001.starch")))
        | filter { !it[3].exists() }
        | map(it -> tuple(it[0], it[1], it[2]))
        | callHotspots
}

// workflow hotspots {
// 	params.basepath = "/net/seq/data2/projects/sabramov/ENCODE4/atac_aligning/output"
// 	metadata = Channel.fromPath(params.samples_file)
// 		| splitCsv(header:true, sep:'\t')
// 		| map(row -> tuple(
// 			row.sample_id, 
// 			file("${params.basepath}/${row.sample_id}/${row.sample_id}.filtered.cram"),
// 			file("${params.basepath}/${row.sample_id}/${row.sample_id}.filtered.cram.crai"),
// 			))
// 		| filter { it[1].exists() }
// 	    | callHotspots
// }

// process calc_different_fdr {
// 	tag "${id}"
//     publishDir "${params.outdir}/${id}", pattern: "${name}"
// 	label 'high_mem'
//     module "modwt/1.0:kentutil/302:bedops/2.4.35-typical:bedtools/2.25.0:samtools/1.3:hotspot2/2.1.1"

//     input:
// 	    tuple val(id), path(allcalls_file)

//     output:
//         tuple val(id), path(name)
    
//     script:
//     name = "${id}.hotspots.fdr0.05.starch"
//     """
//     hsmerge.sh -f 0.05 -m 50 ${allcalls_file} ${name}
//     """
// }

// workflow tmp {
//     Channel.fromPath(params.samples_file)
//         | splitCsv(header:true, sep:'\t')
// 		| map(row -> row.ag_id)
//         | map(it -> tuple(it, file("/net/seq/data2/projects/sabramov/SuperIndex/dnase-index0415/matrices/downsampled_projected_atac/test_atac_peaks/output/${it}.allcalls.starch")))
// }