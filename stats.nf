include { get_container } from "./helpers"
include { filter_nuclear; total_bam_stats } from "./aligning"

fastaContainer = get_container(params.genome_fasta_file)


process run_preseq {
    conda "/home/sabramov/miniconda3/envs/super-index"
    tag "${sample_id}:${read_type}"
    publishDir "${params.outdir}"

    input:
        tuple val(sample_id), path(cram_file), path(cram_file_index), val(read_type)
    
    output:
        tuple val(sample_id), path(name)
    
    script:
    name = "${sample_id}.preseq_hist.txt"
    flag = (read_type == "paired") ? "-P" : ""
    """
    samtools view \
        -h ${cram_file} \
        \$(cat ${params.nuclear_chroms} | tr '\n' ' ') \
        | preseq c_curve -B ${flag} -v - 2>&1 > ${name}
    """
}



workflow nuclearStats {
    Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
            row.sample_id,
            file(row.cram_file),
            file(row.cram_index),
            )
        )
        | filter_nuclear
        | total_bam_stats
}

workflow preseq {
    Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
            row.sample_id,
            file(row.cram_file),
            file(row.cram_index),
            row.read_type
            )
        )
        | run_preseq
}


process percent_dup {
    scratch true
    container "${params.container}"
    containerOptions "${fastaContainer}"
    tag "${sample_id}"
    cpus 2
    publishDir "${params.outdir}/${sample_id}"

    input:
        tuple val(sample_id), path(bam_file), path(bam_file_index)
    
    output:
        tuple val(sample_id), path(name), path(perc_dup)

    script:
    name = "${sample_id}.spotdups.txt"
    perc_dup = "${sample_id}.percent_dup.txt"
    """
    samtools view \
        -u \
        --reference ${params.genome_fasta_file} \
        ${bam_file} \
        | samtools sort \
        -@${task.cpus} > sorted.bam
    samtools index sorted.bam

    picard RevertSam \
        INPUT=sorted.bam \
        OUTPUT=clear.bam \
        RESTORE_HARDCLIPS=false \
        VALIDATION_STRINGENCY=SILENT \
        REMOVE_DUPLICATE_INFORMATION=true\
        SORT_ORDER=coordinate \
        RESTORE_ORIGINAL_QUALITIES=false \
        REMOVE_ALIGNMENT_INFORMATION=false

    picard MarkDuplicatesWithMateCigar \
        INPUT=clear.bam \
        METRICS_FILE=${name} \
        OUTPUT=/dev/null \
        ASSUME_SORTED=true \
        MINIMUM_DISTANCE=300 \
        VALIDATION_STRINGENCY=SILENT \
        READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
    
    grep -A 1 "PERCENT_DUPLICATION" ${name} \
        | awk -F'\t' '{print \$(NF-1)}' \
        | tail -n +2 > ${perc_dup}
    """
}

process collect_basic_stats {
    container "${params.container}"
    containerOptions "${fastaContainer}"
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}"

    input:
        tuple val(sample_id), path(cram_file), path(cram_file_index)

    output:
        tuple val(sample_id), path(name)
    
    script:
    name = "${sample_id}.sequencing_stats.txt"
    """
    echo -e "sample_id\tname\tvalue" > ${name}
    echo -e "${sample_id}\tfiltered_aligned\t\$(samtools view --reference ${params.genome_fasta_file} -c ${cram_file})" >> ${name}
    echo -e "${sample_id}\tduplicates\t\$(samtools view -f 1024 --reference ${params.genome_fasta_file} -c ${cram_file})" >> ${name}
    """
}


workflow {
    Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
            row.sample_id,
            file(row.cram_file),
            file(row.cram_index),
            )
        )
        | (collect_basic_stats & percent_dup)
}