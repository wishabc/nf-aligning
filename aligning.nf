#!/usr/bin/env nextflow
include { get_container; set_key_for_group_tuple } from "./helpers"
nextflow.enable.dsl = 2


nuclearChromsContainer = get_container(params.nuclear_chroms)
chromSizesContainer = get_container(params.chrom_sizes)
fastaContainer = get_container(params.genome_fasta_file)

genome_fasta_file = file(params.genome_fasta_file)



process align_reads_single {
  cpus params.threads
  tag "${group_key}:${name}"
  scratch true
  container "${params.container}"
  containerOptions fastaContainer
  
  input:
    tuple val(group_key), path(trimmed_r1)

  output:
    tuple val(group_key), path(name)

  script:
  name = "${trimmed_r1.baseName}.bam"
  bwa_index = "${genome_fasta_file.parent}/${genome_fasta_file.baseName}"
  """
  bwa aln \
    -Y -l 32 -n 0.04 \
    -t "${task.cpus}" \
    "${bwa_index}" \
    "${trimmed_r1}" \
    > out.sai

  bwa samse \
    -n 10 \
    "${bwa_index}" \
    out.sai \
    "${trimmed_r1}" \
  | samtools view -b -T ${params.genome_fasta_file} - \
  > ${name}
  """
}

process align_reads_paired {

  cpus params.threads
  tag "${group_key}:${name}"
  scratch true
  container "${params.container}"
  containerOptions fastaContainer 

  input:
    tuple val(group_key), path(trimmed_r1), path(trimmed_r2)

  output:
    tuple val(group_key), path(name)

  script:
  name = "${trimmed_r1.baseName}.bam"
  bwa_index = "${genome_fasta_file.parent}/${genome_fasta_file.baseName}"
  """
  bwa aln \
    -Y -l 32 -n 0.04 \
    -t "${task.cpus}" \
    "${bwa_index}" \
    "${trimmed_r1}" \
    > out1.sai

  bwa aln \
    -Y -l 32 -n 0.04 \
    -t "${task.cpus}" \
    "${bwa_index}" \
    "${trimmed_r2}" \
    > out2.sai
  bwa sampe \
    -n 10 -a 750 \
    "${bwa_index}" \
    out1.sai out2.sai \
    "${trimmed_r1}" "${trimmed_r2}" \
  | samtools view -b -T ${params.genome_fasta_file} - \
  > ${name}
  """
}

process filter_and_sort {
    scratch true
    cpus params.threads
    tag "${group_key}"
    container "${params.container}"
    containerOptions nuclearChromsContainer

    input:
        tuple val(group_key), path(bam_file)

    output:
        tuple val(group_key), path(name)

    script:
    name = "${bam_file.baseName}.sorted.bam"
    """
    # filter
    python3 $projectDir/bin/filter_reads.py \
        ${bam_file} \
        filtered.bam \
        ${params.nuclear_chroms}
    # sort
    samtools sort -@"${task.cpus}" filtered.bam > ${name}
    """
}
/*
 * Step 3: Merge alignments into one big ol' file
 */
process merge_bam {
    tag "${group_key}"
    container "${params.container}"
    scratch true
    cpus 2
    errorStrategy "ignore"

    input:
        tuple val(group_key), path(bamfiles)

    output:
        tuple val(group_key), file(name), path("${name}.bai")

    script:
    name = "${group_key}.bam"
    if (group_key.size > 1) {
    """
    samtools merge merged.bam ${bamfiles}
    samtools sort -@"${task.cpus}" merged.bam > ${name}
    samtools index ${name}
    """
    } else {
    """
    ln -s ${bamfiles} ${name} 
    samtools index ${name}
    """
    }
}

/*
 * Step 4: Mark duplicates with Picard
 */
process mark_duplicates {
  tag "${ag_id}"
  scratch true
  label "high_mem"

  publishDir "${params.outdir}/${ag_id}/stats", pattern: "${metric_name}"
  container "${params.container}"

  input:
    tuple val(ag_id), path(merged_bam), path(merged_bam_index)

  output:
    tuple val(ag_id), path(name), path("${name}.bai"), path(metric_name)
  
  script:
  name = "${ag_id}.marked.bam"
  metric_name = "${ag_id}.MarkDuplicates.picard"
    """
    picard RevertOriginalBaseQualitiesAndAddMateCigar \
        INPUT="${merged_bam}" OUTPUT=cigar.bam \
        VALIDATION_STRINGENCY=SILENT RESTORE_ORIGINAL_QUALITIES=false \
        SORT_ORDER=coordinate MAX_RECORDS_TO_EXAMINE=0 \
        TMP_DIR=${workDir}

    picard MarkDuplicatesWithMateCigar \
        INPUT=cigar.bam OUTPUT=${name} \
        MINIMUM_DISTANCE=300 \
        TMP_DIR=${workDir} \
        METRICS_FILE=${metric_name} ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
        READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
    samtools index ${name}
  """
}

/**
Step 5: Filter down to nuclear reads passing filter
**/
process filter_nuclear {
    container "${params.container}"
    containerOptions nuclearChromsContainer
    tag "${ag_id}"
    cpus 3
    scratch true

    input:
        tuple val(ag_id), path(bam), path(bam_index)

    output:
        tuple val(ag_id), path("${name}"), path("${name}.bai")

    script:
    name = "${ag_id}.filtered.bam"
    """
    if [[ "${params.save_cram_mode}" == "nuclear" ]]; then
        samtools view -@ ${task.cpus} -b \
            -F 516 \
            ${bam} \
            \$(cat "${params.nuclear_chroms}") \
            -o ${name}
    else
        samtools view -@ ${task.cpus} -b -F 516 ${bam} -o ${name}
    fi

    samtools index -@ ${task.cpus} ${name}
    """
}


// Works only with paired end data
process insert_size {
  tag "${ag_id}"
  scratch true
  container "${params.container}"
  publishDir "${params.outdir}/${ag_id}/stats"

  input:
    tuple val(ag_id), path(bam), path(bai), val(is_paired)

  output:
    tuple val(ag_id), path(stats_name), path(pdf_name)

  when:
    is_paired

  script:
  stats_name = "${ag_id}.CollectInsertSizeMetrics.picard"
  pdf_name = "${stats_name}.pdf"
  """
  samtools idxstats "${bam}" \
  | cut -f 1 \
  | grep -v chrM \
  | grep -v chrC \
  | xargs samtools view -b "${bam}" -o nuclear.bam
  touch ${stats_name}
  touch ${pdf_name}
  picard CollectInsertSizeMetrics \
    INPUT=nuclear.bam \
    OUTPUT=${stats_name} \
    TMP_DIR=${workDir} \
    HISTOGRAM_FILE=${pdf_name} \
    VALIDATION_STRINGENCY=LENIENT \
    ASSUME_SORTED=true
  """
}

process total_bam_stats {
    container "${params.container}"
    containerOptions "${fastaContainer}"
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
        ${name} \
        --reference ${params.genome_fasta_file}
    """
}


process macs2 {
  tag "${ag_id}"
  publishDir "${params.outdir}/${ag_id}/stats"
  container "${params.container}"
  containerOptions "${chromSizesContainer}"
  scratch true

  when:
    params.do_macs

  input:
    tuple val(ag_id), path(bam), path(bai), val(is_paired)

  output:
    tuple val(ag_id), path("macs2.${ag_id}*")

  script:
  mode = 'BAM'
  """
  gen_size=\$(cat ${params.chrom_sizes} | awk '{print \$2}' | paste -sd+ | bc)
  macs2 callpeak \
    -t "${bam}" \
    -f ${mode} \
    -n macs2.${ag_id}\
    -g \$gen_size -p 0.01 \
    --shift 75 --extsize 150 \
    --nomodel -B --SPMR \
    --keep-dup all --call-summits > "macs2.${ag_id}.err"
  """
}


/**
Step 6: Convert Filtered Bam to cram file
**/
process convert_to_cram {
  tag "${ag_id}"
  publishDir "${params.outdir}/${ag_id}"
  cpus params.threads
  container "${params.container}"
  containerOptions fastaContainer

  input:
    tuple val(ag_id), path(bam), path(bam_index)

  output:
    tuple val(ag_id), path(cramfile), path("${cramfile}.crai")

  script:
  cramfile = bam.baseName + ".cram"
  """
  samtools view "${bam}" \
    -C -O cram,version=3.0,level=7,lossy_names=0 \
    -T "${params.genome_fasta_file}" \
    --threads "${task.cpus}" \
    --write-index \
    -o "${cramfile}"
  """
}


workflow alignBwa {
    take:
        trimmed_reads
    main:
        reads_divided = trimmed_reads.branch{
            paired: it[3]
            single: true
        }
        paired_bam = reads_divided.paired
            | map(it -> tuple(it[0], it[1], it[2]))
            | align_reads_paired

        all_bam = reads_divided.single 
            | map(it -> tuple(it[0], it[1]))
            | align_reads_single
            | mix(paired_bam)
    emit:
        all_bam
}


workflow alignReads {
    take:
        trimmed_reads
    main:
        marked_bam_files = trimmed_reads
            | alignBwa
            | filter_and_sort
            | groupTuple()
            | merge_bam
            | mark_duplicates
            | map(it -> tuple(it[0], it[1], it[2]))
        

        marked_bam_files
            | total_bam_stats

        is_paired_dict = trimmed_reads
            | map(it -> tuple(it[0], it[3]))
            | distinct()
        
        marked_bam_files
            | filter_nuclear
            | join(is_paired_dict)
            | (insert_size & macs2)
    
        filter_nuclear.out
            | convert_to_cram
    emit:
        convert_to_cram.out
}

workflow {
    fastq_trimmed_paired = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
	    | map(row -> tuple(row.ag_id, row.reads1, row.reads2, row.type == 'paired'))
        | set_key_for_group_tuple
        | alignReads
}
