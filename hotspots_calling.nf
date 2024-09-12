#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { get_container } from "./aligning"


process spot_score {

    // Impossible to use anywhere except Altius cluster
    //conda params.conda
    module "bedops/2.4.35-typical:samtools/1.3:modwt/1.0:kentutil/302:hotspot2/2.1.1:jdk/1.8.0_92:gcc/4.7.2:R/3.2.5:picard/2.8.1:git/2.3.3:coreutils/8.25:bedtools/2.25.0:python/3.5.1:pysam/0.9.0:htslib/1.6.0:numpy/1.11.0:atlas-lapack/3.10.2:scipy/1.0.0:scikit-learn/0.18.1:preseq:/2.0.3:gsl/2.4"
    publishDir "${params.outdir}/${ag_id}"


    input:
        tuple val(ag_id), path(bam_file), path(bam_file_index)

    output:
        tuple val(ag_id), path("r1.*")

    script:
    renamed_input = "r1.${bam_file.extension}"
    genome_file = file(params.genome_fasta_file)
    genome_prefix = "${genome_file.parent}/${genome_file.simpleName}"

    """
    # workaround for hotspots1 naming scheme...
    # might need to check bam conversion
	ln -s ${bam_file} ${renamed_input}
	ln -s ${bam_file_index} ${renamed_input}.${bam_file_index.extension}
    bash $moduleDir/bin/runhotspot.bash \
      "${params.hotspots_dir}" \
      "\$PWD" \
      "${renamed_input}" \
      "${genome_prefix}" \
      ${params.chrominfo} \
      "${params.readlength}" \
      DNaseI

    starch --header \
        r1-both-passes/r1.hotspot.twopass.zscore.wig > r1.spots.starch

    bash $moduleDir/bin/info.sh \
      r1.spots.starch hotspot1 r1.spot.out \
      > r1.hotspot.info
    """
}


// Need to contenirize at some point
process call_hotspots {
	tag "${id}:${max_fdr}"
	label 'high_mem'
	publishDir "${params.outdir}/${id}", pattern: "${id}*"

	//container "${params.container}"
	//containerOptions "${get_container(params.nuclear_chroms)} ${get_container(params.chrom_sizes_bed)} ${get_container(params.mappable)} ${get_container(params.centers)}"
	module "modwt/1.0:kentutil/302:bedops/2.4.35-typical:bedtools/2.25.0:samtools/1.3:hotspot2/2.1.1"

	input:
	    tuple val(id), path(bam_file), path(bam_file_index), val(max_fdr)

	output:
        tuple val(id), val(max_fdr), path("${id}.allcalls.starch"), path("${id}.cutcounts.starch"), path("${id}.cleavage.total"), emit: peak_calling
	    tuple val(id), val(max_fdr), path("${id}.SPOT.txt"), path("${id}.density.bw"), emit: extra
        tuple val(id), val(max_fdr), path(hotspots), path(peaks), emit: hotspots

	script:
    hotspots = "${id}.hotspots.fdr${max_fdr}.starch"
    peaks = "${id}.peaks.fdr${max_fdr}.starch"
	renamed_input = "nuclear.bam"
	"""
    # hotspot2 needs a tmp directory, scratch doesn't work due to limited tmp space
    mkdir -p tmp
	export TMPDIR="\$PWD/tmp"

	# workaround for hotspots2 naming scheme
	ln -sf ${bam_file} ${renamed_input}
	ln -sf ${bam_file_index} ${renamed_input}.bai

	hotspot2.sh -F ${max_fdr} -f ${max_fdr} \
		-p "varWidth_20_${id}" \
		-M "${params.mappable}" \
		-c "${params.chrom_sizes_bed}" \
		-C "${params.centers}" \
		${renamed_input} \
		'.'

    mv nuclear.hotspots.fdr${max_fdr}.starch ${hotspots}
    mv nuclear.peaks.starch ${peaks}
    mv nuclear.SPOT.txt ${id}.SPOT.txt
    mv nuclear.cleavage.total ${id}.cleavage.total
    mv nuclear.density.bw ${id}.density.bw
    mv nuclear.cutcounts.starch ${id}.cutcounts.starch
    mv nuclear.allcalls.starch ${id}.allcalls.starch

    rm -r tmp
	"""
}


process hotspots_other_fdr {
	tag "${id}:${fdr}"
	label 'high_mem'
	publishDir "${params.outdir}/${id}"

	//container "${params.container}"
	//containerOptions "${get_container(params.nuclear_chroms)} ${get_container(params.chrom_sizes_bed)} ${get_container(params.mappable)} ${get_container(params.centers)}"
	module "modwt/1.0:kentutil/302:bedops/2.4.35-typical:bedtools/2.25.0:samtools/1.3:hotspot2/2.1.1"

    input:
	    tuple val(id), path(all_calls), path(cutcounts), path(total_counts), val(fdr)

    output:
        tuple val(id), val(fdr), path(hotspots), path(peaks)
    
    script:
    hotspots = "${id}.hotspots.fdr${fdr}.starch"
    peaks = "${id}.peaks.fdr${fdr}.starch"
    """
    # hotspot2 needs a tmp directory, scratch doesn't work due to limited tmp space
    mkdir -p tmp 
	export TMPDIR="\$PWD/tmp"

    hsmerge.sh -f ${fdr} -m 50 ${all_calls} ${hotspots}
    
    bash density-peaks.bash \
        \$TMPDIR \
        "varWidth_20_${id}" \
        ${cutcounts} \
        ${hotspots} \
        ${params.chrom_sizes_bed} \
        \$TMPDIR/nuclear.density.${fdr}.starch \
        ${peaks} \
        `cat ${total_counts}`
    
    rm -r tmp
    """
}

workflow callHotspots {
	take:
		bam_files
	main:
        fdrs = Channel.of(params.hotspot2_fdr)
            | flatMap(it -> it.tokenize(','))
            | map { it as Float }
        
        data = bam_files
            | combine(fdrs.max())
            | call_hotspots
        
        extra_fdrs = fdrs
            | combine(fdrs.max())
            | filter { it[0] != it[1] }
            | view()
            | map(it -> it[0])

        out = data.peak_calling
            | map(it -> tuple(it[0], it[4], it[5], it[2]))
            | combine(extra_fdrs)
            | hotspots_other_fdr
            | mix(data.hotspots)
	emit:
		out

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