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
	tag "${id}"
	//label 'high_mem'
	publishDir "${params.outdir}/${id}", pattern: "${id}*"
    memory { 80.GB + 20.GB * task.attempt }

	cpus 6
    conda "/home/sabramov/miniconda3/envs/jupyterlab"

	input:
	    tuple val(id), path(bam_file), path(bam_file_index)

	output:
        tuple val(id), path("${id}.*"), emit: all
        tuple val(id), path("${id}.*.parquet"), emit: tmp_files
        tuple val(id), path("${id}.hotspots*"), path("${id}.peaks*"), emit: peaks_hotspots

	script:
    fdrs = params.hotspot2_fdr.tokenize(',').join(' ')
	"""
    python3 ~/packages/hotspot2/hotspot2/track_memory.py \
        ${id} \
        --bam ${bam_file} \
        --fdrs ${fdrs} \
        --mappable_bases /net/seq/data2/projects/sabramov/SuperIndex/GRCh38_no_alts.K36.center_sites.n100.nuclear.merged.bed.gz \
        --chrom_sizes ${params.chrom_sizes}  \
        --cpus ${task.cpus} \
        --save_density

    # TODO rm pvals smoothed_signal parquets
	"""
}

workflow callHotspots {
	take:
		bam_files
	main:
        out = bam_files
            | call_hotspots
	emit:
		out.peaks_hotspots

}


workflow {
	Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
		| map(row -> tuple( row.ag_id, row.cram_file, row.cram_index ?: "${row.cram_file}.crai"))
        | callHotspots
}


// defunc

workflow hotspotLowerFdr {
    fdrs = Channel.of(params.hotspot2_fdr)
            | flatMap(it -> it.toString().tokenize(','))
            | map { it as Float }
    Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
		| map(
            row -> tuple(
                row.ag_id,
                file(row.allcalls),
                file(row.cutcounts),
                file(row.cleavage_total)
            )
        )
        | combine(fdrs)
        | map(it -> tuple(it[4], *it[0..3]))
        | hotspots_other_fdr
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