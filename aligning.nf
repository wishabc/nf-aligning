#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Workaround, so when we groupTuple later, 
// it knows how many objects in the group are going to be
def set_key_for_group_tuple(ch) {
  s = ch.groupTuple()
  s.view()
  a = s.map{ key, files -> tuple(groupKey(key, files.size()), files) }
  a.view()
  a.transpose()
}


process align_reads_single {
  cpus params.threads
  tag "${group_key}:${name}"
  scratch true
  container "${params.container}"

  input:
    tuple val(group_key), path(trimmed_r1)

  output:
    tuple val(group_key), path(name)

  script:
  name = "${trimmed_r1.baseName}.bam"
  """
  bwa aln \
    -Y -l 32 -n 0.04 \
    -t "${task.cpus}" \
    "${params.genome_fasta_file}" \
    "${trimmed_r1}" \
    > out.sai

  bwa samse \
    -n 10 \
    "${params.genome_fasta_file}" \
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

  input:
    tuple val(group_key), path(trimmed_r1), path(trimmed_r2)

  output:
    tuple val(group_key), path(name)

  script:
  name = "${trimmed_r1.baseName}.bam"
  """
  bwa aln \
    -Y -l 32 -n 0.04 \
    -t "${task.cpus}" \
    "${params.genome_fasta_file}" \
    "${trimmed_r1}" \
    > out1.sai

  bwa aln \
    -Y -l 32 -n 0.04 \
    -t "${task.cpus}" \
    "${params.genome_fasta_file}" \
    "${trimmed_r2}" \
    > out2.sai
  bwa sampe \
    -n 10 -a 750 \
    "${params.genome_fasta_file}" \
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

  input:
    tuple val(group_key), path(bam_file)

  output:
    tuple val(group_key), path(name)
  
  script:
  name = "${align_id}.sorted.bam"
  """
  # filter
  python3 $projectDir/bin/filter_reads.py \
    ${bam_file} \
    filtered.bam \
    ${params.nuclear_chroms}
  # sort
  samtools sort \
    -l 0 -m 1G -@"${task.cpus}" filtered.bam \
    > ${name}
  """
}
/*
 * Step 3: Merge alignments into one big ol' file
 */
process merge_bam {
  tag "${group_key}"
  container "${params.container}"
  scratch true

  input:
    tuple val(group_key), path(bamfiles)

  output:
    tuple val(group_key), file(name), path("${name}.bai")

  script:
  name = "${group_key}.bam"
  """
  samtools merge ${name} ${bamfiles}
  samtools index ${name}
  """
}

/*
 * Step 4: Mark duplicates with Picard
 */
process mark_duplicates {
  tag "${sample_id}"
  scratch true
  publishDir "${params.outdir}/${sample_id}/stats", pattern: "${metric_name}"
  container "${params.container}"

  input:
    tuple val(sample_id), path(merged_bam), path(merged_bam_index)

  output:
    tuple val(sample_id), path(name), path("${name}.bai"), path(metric_name)
  
  script:
  name = "${sample_id}.marked.bam"
  metric_name = "${sample_id}.MarkDuplicates.picard"
  """
  picard RevertOriginalBaseQualitiesAndAddMateCigar \
    INPUT="${merged_bam}" OUTPUT=cigar.bam \
    VALIDATION_STRINGENCY=SILENT RESTORE_ORIGINAL_QUALITIES=false SORT_ORDER=coordinate MAX_RECORDS_TO_EXAMINE=0
  picard MarkDuplicatesWithMateCigar \
      INPUT=cigar.bam OUTPUT=${name} \
      MINIMUM_DISTANCE=300 \
      METRICS_FILE=${metric_name} ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
      READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
  samtools index ${name}
  """
}

/**
Step 5: Filter down to nuclear reads passing filter
**/
process filter {
  container "${params.container}"
  tag "${sample_id}"

  input:
    tuple val(sample_id), path(bam), path(bam_index), path(picard_dup_file) 

  output:
    tuple val(sample_id), path("${name}"), path("${name}.bai")

  script:
  name = "${sample_id}.filtered.bam"
  """
  samtools view -b -F 512 ${bam} > filtered.bam
  samtools index filtered.bam
  cat "${params.nuclear_chroms}" \
  | xargs samtools view -b filtered.bam > ${name}

  samtools index ${name}
  """
}
// Works only with paired end data
process insert_size {
  tag "${sample_id}"
  publishDir "${params.outdir}"
  scratch true
  container "${params.container}"
  publishDir "${params.outdir}/${sample_id}/stats"

  input:
    tuple val(sample_id), path(bam), path(bai), val(is_paired)

  output:
    tuple val(sample_id), path(stats_name), path(pdf_name)

  when:
    is_paired

  script:
  stats_name = "${sample_id}.CollectInsertSizeMetrics.picard"
  pdf_name = "${stats_name}.pdf"
  """
  samtools idxstats "${bam}" \
  | cut -f 1 \
  | grep -v chrM \
  | grep -v chrC \
  | xargs samtools view -b "${bam}" -o nuclear.bam
  picard CollectInsertSizeMetrics \
    INPUT=nuclear.bam \
    OUTPUT=${stats_name} \
    HISTOGRAM_FILE=${pdf_name} \
    VALIDATION_STRINGENCY=LENIENT \
    ASSUME_SORTED=true
  """
}

process density_files {
  publishDir "${params.outdir}/${sample_id}"
  tag "${sample_id}"
  container "${params.container}"

  input:
    tuple val(sample_id), path(bam), path(bai)

  output:
    tuple val(sample_id), path("${sample_id}.density.bed.starch"), emit: starch
    tuple val(sample_id), path("${sample_id}.density.bw"), emit: bigwig
    tuple val(sample_id), path("${sample_id}.density.bed.bgz"), emit: bgzip

  script:
    """
    bam2bed -d \
      < "${bam}" \
      | cut -f1-6 \
      | awk '{ if( \$6=="+" ){ s=\$2; e=\$2+1 } else { s=\$3-1; e=\$3 } print \$1 "\t" s "\t" e "\tid\t" 1 }' \
      | sort-bed - \
      > sample.bed
    
    unstarch "${params.density_buckets}" \
      | bedmap --faster --echo --count --delim "\t" - sample.bed \
      | awk -v binI=${params.density_step_size} -v win="${params.density_window_width}" \
          'BEGIN{ halfBin=binI/2; shiftFactor=win-halfBin } { print \$1 "\t" \$2 + shiftFactor "\t" \$3-shiftFactor "\tid\t" i \$4}' \
      | starch - \
      > density.bed.starch
      
    unstarch density.bed.starch | awk -v binI="${params.density_step_size}" -f "${moduleDir}/awk/bedToWig.awk" > density.wig
    wigToBigWig -clip density.wig "${params.genome_fasta_file}" density.bw
    
    unstarch density.bed.starch | bgzip > density.bed.bgz
    tabix -p bed density.bed.bgz
    """
}


/**
Step 6: Convert Filtered Bam to cram file
**/
process convert_to_cram {
  tag "${sample_id}"
  publishDir "${params.outdir}/${sample_id}"
  cpus params.threads
  container "${params.container}"

  input:
    tuple val(sample_id), path(bam), path(bam_index)

  output:
    tuple val(sample_id), path(cramfile), path("${cramfile}.crai")

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
      paired: it[4]
      single: true
    }
    paired_bam = align_reads_paired(reads_divided.paired.map(it -> tuple(it[0], it[1], it[2])))
    single_bam = align_reads_single(reads_divided.single.map(it -> tuple(it[0], it[1])))
    all_bam = paired_bam.mix(single_bam)
  emit:
    all_bam
}


workflow alignReads {
  take:
    trimmed_reads
  main:
    aligned_files = trimmed_reads | alignBwa | filter_and_sort
    filtered_bam_files = merge_bam(aligned_files.groupTuple()) 
    | mark_duplicates 
    | filter

    is_paired_dict = trimmed_reads.map(it -> tuple(it[0], it[4])).distinct()
    
    insert_size(filtered_bam_files.join(is_paired_dict))
    density_files(filtered_bam_files)
    convert_to_cram(filtered_bam_files)
  emit:
    convert_to_cram.out
}

workflow {
    fastq_trimmed_paired = Channel
      .fromPath(params.samples_file)
      .splitCsv(header:true, sep:'\t')
		  .map(row -> tuple( row.sample_id, row.reads1, row.reads2, row.is_paired))
    alignReads(set_key_for_group_tuple(fastq_trimmed_paired))
}