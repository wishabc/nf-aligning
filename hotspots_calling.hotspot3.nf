#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { get_container } from "./aligning"


// Need to containerize at some point
process call_hotspots {
	tag "${id}"
	//label 'high_mem'
	publishDir "${params.outdir}/${id}"
    memory { 60.GB + 40.GB * task.attempt }

	cpus 15
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
