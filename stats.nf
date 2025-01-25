process percent_dup {
    scratch true
    container "${params.container}"
    tag "${ag_id}"
    cpus 2
    publishDir "${params.outdir}/${ag_id}"

    input:
        tuple val(ag_id), path(bam_file), path(bam_file_index)
    
    output:
        tuple val(ag_id), path(name), path(perc_dup)

    script:
    name = "${ag_id}.spotdups.txt"
    perc_dup = "${ag_id}.percent_dup.txt"
    """
    samtools sort -@${task.cpus} ${bam_file} > sorted.bam
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
    tag "${ag_id}"
    publishDir "${params.outdir}/${ag_id}"

    input:
        tuple val(ag_id), path(cram_file), path(cram_file_index)

    output:
        tuple val(ag_id), path(name)
    
    script:
    name = "${ag_id}.sequencing_stats.txt"
    """
    echo -e "ag_id\tname\tvalue" > ${name}
    echo -e "${ag_id}\tfiltered_aligned\t\$(samtools view -c ${cram_file})" >> ${name}
    echo -e "${ag_id}\tduplicates\t\$(samtools view -f 1024 -c ${cram_file})" >> ${name}
    """
}

process total_bam_stats {
    conda params.conda
    tag "${ag_id}"
    publishDir "${params.outdir}/${ag_id}"

     input:
        tuple val(ag_id), path(bam_file), path(bam_index)

    output:
        tuple val(ag_id), path(name)
    
    script:
    name = "${ag_id}.total_sequencing_stats.txt"
    """
    python3 $moduleDir/bin/bamcounts.py \
        ${bam_file} \
        ${name}
    """
}

process run_preseq {
    conda "/home/sabramov/miniconda3/envs/super-index"
    tag "${ag_id}:${read_type}"
    publishDir "${params.outdir}"

    input:
        tuple val(ag_id), path(cram_file), path(cram_file_index), val(read_type)
    
    output:
        tuple val(ag_id), path(name)
    
    script:
    name = "${ag_id}.preseq_hist.txt"
    flag = (read_type == "paired") ? "-P" : ""
    """
    samtools view \
        -h ${cram_file} \
        \$(cat ${params.nuclear_chroms} | tr '\n' ' ') \
        | preseq c_curve -B ${flag} -v - 2>&1 > ${name}
    """
}

workflow {
    Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
            row.ag_id,
            file(row.cram_file),
            file(row.cram_index),
            )
        )
        | (collect_basic_stats & percent_dup)
}

workflow preseq {
    Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
            row.ag_id,
            file(row.cram_file),
            file(row.cram_index),
            row.read_type
            )
        )
        | run_preseq
}