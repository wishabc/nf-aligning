#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { set_key_for_group_tuple } from "./helpers"

// TODO check publishDirs

process align_reads_single {
  cpus params.threads
  tag "${group_key}:${name}"
  scratch true

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

process filter {
  scratch true
  cpus params.threads
  tag "${group_key}"

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
// TODO: single end
process mark_duplicates {

  tag "${sample_id}"
  scratch true

  conda params.conda

  input:
    tuple val(sample_id), path(merged_bam), path(merged_bam_index)

  output:
    tuple val(sample_id), path(name), path("${name}.bai"), path("MarkDuplicates.picard")
  
  script:
  name = "${sample_id}.marked.bam"
  cmd = "MarkDuplicatesWithMateCigar"
  extra = "MINIMUM_DISTANCE=300"
  """
  picard RevertOriginalBaseQualitiesAndAddMateCigar \
    "INPUT=${merged_bam}" OUTPUT=cigar.bam \
    VALIDATION_STRINGENCY=SILENT RESTORE_ORIGINAL_QUALITIES=false SORT_ORDER=coordinate MAX_RECORDS_TO_EXAMINE=0
  picard "${cmd}" \
      INPUT=cigar.bam OUTPUT=${name}.bam \
      $extra \
      METRICS_FILE=MarkDuplicates.picard ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
      READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
  samtools index ${name}
  """
}

/**
Step 5: Filter down to nuclear reads passing filter
**/
process filter {
  conda params.conda
  tag "$sample_id"

  input:
    tuple val(sample_id), path(bam), path(bam_index), path(picard_dup_file) 

  output:
    file val(sample_id), path("${name}"), path("${name}.bai")

  script:
  flag = 512
  name = "${sample_id}.filtered.bam"
  """
  samtools view -b -F "${flag}" ${bam} > filtered.bam
  samtools index filtered.bam
  cat "${params.nuclear_chroms}" \
  | xargs samtools view -b filtered.bam \
  > ${name}
  samtools index ${name}
  """
}
// Works only with paired end data
process insert_size {
  // TODO: Fix chomosome removal!!
  tag "${sample_id}"
  publishDir params.outdir
  scratch true

  input:
    tuple val(sample_id), path(bam), path(bai), val(is_paired)

  output:
    tuple val(sample_id), path('CollectInsertSizeMetrics.picard'), emit: plaintext
    tuple val(sample_id), path('CollectInsertSizeMetrics.picard.pdf'), emit: pdf

  when:
    is_paired

  script:
  """
  samtools idxstats "${bam}" \
  | cut -f 1 \
  | grep -v chrM \
  | grep -v chrC \
  | xargs samtools view -b "${bam}" -o nuclear.bam
  picard CollectInsertSizeMetrics \
    INPUT=nuclear.bam \
    OUTPUT=CollectInsertSizeMetrics.picard \
    HISTOGRAM_FILE=CollectInsertSizeMetrics.picard.pdf \
    VALIDATION_STRINGENCY=LENIENT \
    ASSUME_SORTED=true
  """
}

process density_files {
  publishDir params.outdir
  tag "${sample_id}"

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
      
    unstarch density.bed.starch | awk -v binI="${params.density_step_size}" -f "\$STAMPIPES/awk/bedToWig.awk" > density.wig
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
  conda params.conda
  publishDir params.outdir
  cpus params.threads

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
    aligned_files = set_key_for_group_tuple(trimmed_reads) | alignBwa
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
    alignReads(fastq_trimmed_paired)
}