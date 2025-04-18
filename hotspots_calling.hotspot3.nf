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


// Need to containerize at some point
process call_hotspots {
	tag "${id}"
	//label 'high_mem'
	publishDir "${params.outdir}/${id}"
    memory { 60.GB + 20.GB * task.attempt }

	cpus 6
    conda "/home/sabramov/miniconda3/envs/jupyterlab"

	input:
	    tuple val(id), path(bam_file), path(bam_file_index)

	output:
        tuple val(id), path("${id}.*"), emit: all
        tuple val(id), path("fdr*/*"), emit: peaks
        tuple val(id), path("debug"), emit: debug

	script:
    fdrs = params.fdrs.tokenize(',').join(' ')
    save_debug = params.save_debug ? "--debug" : ""
	"""
    hotspot3 \
        ${id} \
        --bam ${bam_file} \
        --fdrs ${fdrs} \
        --mappable_bases ${params.mappable_bases} \
        --chrom_sizes ${params.nuclear_chrom_sizes}  \
        --cpus ${task.cpus} \
        --save_density \
        --reference ${params.genome_fasta_file} \
        --debug 2>&1 > ${id}.peak_calling.log

    if [ "${save_debug}" == "" ]; then
        rm -r debug/*
    fi
	"""
}

process call_hotspots_from_cutcounts {
	tag "${id}"
	//label 'high_mem'
	publishDir "${params.outdir}/${id}"
    memory { 60.GB + 20.GB * task.attempt }

	cpus 6
    conda "/home/sabramov/miniconda3/envs/jupyterlab"

	input:
	    tuple val(id), path(cutcounts), path(cutcounts_index), path(total_cutcounts)

	output:
        tuple val(id), path("${id}.*"), emit: all
        tuple val(id), path("fdr*/*"), emit: peaks
        tuple val(id), path("debug"), emit: debug

	script:
    fdrs = params.fdrs.tokenize(',').join(' ')
    save_debug = params.save_debug ? "--debug" : ""
	"""
    mkdir debug
    hotspot3 \
        ${id} \
        --cutcounts ${cutcounts} \
        --fdrs ${fdrs} \
        --mappable_bases ${params.mappable_bases} \
        --chrom_sizes ${params.nuclear_chrom_sizes}  \
        --tempdir debug \
        --cpus ${task.cpus} \
        --save_density \
        --debug 2>&1 > ${id}.peak_calling.log

    if [ "${save_debug}" == "" ]; then
        rm -r debug/*
    fi
	"""
}

process call_hotspots_from_pvals {
	tag "${id}"
	//label 'high_mem'
	publishDir "${params.outdir}/${id}"
    memory { 50.GB + 20.GB * task.attempt }

	cpus 8
    conda "/home/sabramov/miniconda3/envs/jupyterlab"

	input:
	    tuple val(id), path(cutcounts), path(cutcounts_index), path(total_cutcounts), path(pvals_parquet)

	output:
        tuple val(id), path("${id}.*"), emit: all
        tuple val(id), path("fdr*/*"), emit: peaks
        tuple val(id), path("debug"), emit: debugging

	script:
    fdrs = params.fdrs.tokenize(',').join(' ')
    save_debug = params.save_debug ? "--debug" : ""
	"""
    hotspot3 \
        ${id} \
        --cutcounts ${cutcounts} \
        --fdrs ${fdrs} \
        --mappable_bases ${params.mappable_bases} \
        --chrom_sizes ${params.nuclear_chrom_sizes}  \
        --pvals_parquet ${pvals_parquet} \
        --cpus ${task.cpus} \
        --debug 2>&1 > ${id}.peak_calling.log

    if [ "${save_debug}" == "" ]; then
        rm -r debug/*
    fi
	"""
}

workflow callHotspots {
	take:
		bam_files
	main:
        out = bam_files
            | call_hotspots
	emit:
		out.peaks

}


workflow {
	Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
		| map(row -> tuple(
                row.ag_id,
                row.cram_file,
                row.cram_index ?: "${row.cram_file}.crai",
            )
        )
        | map(it -> tuple(it[0], it[1], it[2]))
        | callHotspots
}


// defunc
workflow tmp_debug {
    Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
		| map(row -> tuple(
                row.ag_id,
                file("${params.outdir}/${row.ag_id}/${row.ag_id}.cutcounts.bed.gz"),
                file("${params.outdir}/${row.ag_id}/${row.ag_id}.cutcounts.bed.gz.tbi"),
                file("${params.outdir}/${row.ag_id}/${row.ag_id}.total_cutcounts"),
                file("${params.outdir}/${row.ag_id}/${row.ag_id}.pvals.parquet")
            )
        )
        | call_hotspots_from_pvals
}

// defunc
workflow fromCutcounts {
    prev_run_dir = "/net/seq/data2/projects/sabramov/SuperIndex/hotspot3/w_babachi_new.v5/output/"
    Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
		| map(row -> tuple(
                row.ag_id,
                file("${prev_run_dir}/${row.ag_id}/${row.ag_id}.cutcounts.bed.gz"),
                file("${prev_run_dir}/${row.ag_id}/${row.ag_id}.cutcounts.bed.gz.tbi"),
                file("${prev_run_dir}/${row.ag_id}/${row.ag_id}.total_cutcounts")
            )
        )
        | call_hotspots_from_cutcounts
}

process extract_pval {
    tag "${id}"
    //publishDir "${params.outdir}/${id}"
    conda "/home/sabramov/miniconda3/envs/jupyterlab"
    label "high_mem"

    input:
        tuple val(id), path(pvals_parquet), path(bed_file)

    output:
        tuple val(id), path(name)

    script:
    name = "${id}.max_pvals.npy"
    """
    echo 1
    hotspot3-pvals \
        ${pvals_parquet} \
        ${bed_file} \
        ${name} \
        --chrom_sizes ${params.nuclear_chrom_sizes} \
        --format npy
    """
}

process create_matrix {
    publishDir params.outdir
    conda "/home/sabramov/miniconda3/envs/jupyterlab"
    label "bigmem"

    input:
        path pvals
    
    output:
        path name
    
    script:
    name = "neglog10_pvals.npy"
    """
    echo "${pvals.join('\n')}" > filelist.txt
    python ${moduleDir}/bin/create_matrix.py \
        ${params.samples_file} \
        filelist.txt \
        ${name} \
    """
    
}

workflow extractMaxPvalue {
    prev_run_dir = "/net/seq/data2/projects/sabramov/SuperIndex/hotspot3/peak_calls.v23/"
    params.bed_file = ""
    Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
                row.ag_id,
                file("${prev_run_dir}/${row.ag_id}/${row.ag_id}.pvals.parquet"),
                file(params.bed_file)
            )
        )
        | extract_pval
        | map(it -> it[1])
        | collect(sort: true, flat: true)
        | create_matrix
}