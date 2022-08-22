#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { set_key_for_group_tuple } from "./helpers"

params.conda = "$moduleDir/environment.yml"

process align_reads {

  cpus params.threads
  conda params.conda
  scratch true

  input:
    tuple val(group_key), val(align_id), file(trimmed_r1), file(trimmed_r2), val(is_paired)

  output:
    tuple val(group_key), val(align_id), path(name)

  script:
  name = "${align_id}.bam"
  if is_paired
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
  else
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

process filter {
  scratch true
  conda params.conda

  cpus params.threads

  input:
    tuple val(group_key), val(align_id), path(bam_file)

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
    -l 0 -m 1G -@ "${task.cpus}" filtered.bam \
    > ${name}
  """
}
/*
 * Step 3: Merge alignments into one big ol' file
 */
process merge_bam {
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
      INPUT=cigar.bam OUTPUT=marked.bam \
      $extra \
      METRICS_FILE=MarkDuplicates.picard ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
      READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
  samtools index marked.bam
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
    file val(sample_id), path("${name}")

  script:
  flag = 512
  name = "${sample_id}.filtered.bam"
  """
  samtools view -b -F "${flag}" ${bam} > filtered.bam
  samtools index filtered.bam
  cat "${params.nuclear_chroms}" \
  | xargs samtools view -b filtered.bam \
  > ${name}
  """
}

/**
Step 6: Convert Filtered Bam to cram file
**/
process convert_to_cram {
  conda params.conda
  publishDir params.outdir
  cpus params.threads

  input:
  tuple val(sample_id), path(bam), path(bam_index)

  output:
  tuple val(sample_id), path(cramfile), path("${cramfile}.crai")

  script:
  cramfile = bam.name.replace("bam", "cram")
  """
  samtools view "${bam}" \
    -C -O cram,version=3.0,level=7,lossy_names=0 \
    -T "${params.genome_fasta_file}" \
    --threads "${task.cpus}" \
    --write-index \
    -o "${cramfile}"
  """
}

workflow alignReads {
  take:
    trimmed_reads
  main:
    aligned_files = set_key_for_group_tuple(bam_files) | align_reads
    marked_dups_files = merge_bam(aligned_files.groupTuple()) | mark_duplicates | filter | cram
  emit:
    convert_to_cram.out
}

workflow {
    fastq_trimmed_paired = Channel
      .fromPath(params.samples_file)
      .splitCsv(header:true, sep:'\t')
		  .map(row -> tuple( row.sample_id, row.align_id, row.reads1, 
      row.reads2, row.is_paired))
    alignReads(fastq_trimmed_paired)
}