#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Workaround, so when we groupTuple later, 
// New key contains info on how many objects are in the group
def set_key_for_group_tuple(ch) {
  ch.groupTuple()
  .map{ it -> tuple(groupKey(it[0], it[1].size()), *it[1..(it.size()-1)]) }
  .transpose()
}

def get_container(file_name) {
  parent = file(file_name).parent
  old_parent = file(file_name).toRealPath().parent
  container = "--bind ${parent},${old_parent}"
}

fastaContainer = get_container(params.genome_fasta_file)
nuclearChromsContainer = get_container(params.nuclear_chroms)
chromSizesContainer = get_container(params.chrom_sizes)
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
  tag "${sample_id}"
  scratch true
  label "high_mem"

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
  tag "${sample_id}"
  scratch true

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
process macs2 {
  tag "${sample_id}"
  publishDir "${params.outdir}/${sample_id}/stats"
  container "${params.container}"
  containerOptions "${chromSizesContainer}"
  scratch true

  when:
    params.do_macs

  input:
    tuple val(sample_id), path(bam), path(bai), val(is_paired)

  output:
    tuple val(sample_id), path("macs2.${sample_id}*")

  script:
  mode = 'BAM'
  """
  gen_size=\$(cat ${params.chrom_sizes} | awk '{print \$2}' | paste -sd+ | bc)
  macs2 callpeak \
    -t "${bam}" \
    -f ${mode} \
    -n macs2.${sample_id}\
    -g \$gen_size -p 0.01 \
    --shift 75 --extsize 150 \
    --nomodel -B --SPMR \
    --keep-dup all --call-summits > "macs2.${sample_id}.err"
  """
}
process density_files {
  publishDir "${params.outdir}/${sample_id}/stats"
  tag "${sample_id}"
  errorStrategy 'ignore'
  container "${params.container}"
  label "high_mem"
  containerOptions "${chromSizesContainer} ${get_container(params.density_buckets)}"

  input:
    tuple val(sample_id), path(bam), path(bai)

  output:
    tuple val(sample_id), path(name), emit: starch
    tuple val(sample_id), path("${sample_id}.density.bw"), emit: bigwig
    tuple val(sample_id), path("${sample_id}.density.bed.bgz"), path("${sample_id}.density.bed.bgz.tbi"), emit: bgzip

  script:
  name = "${sample_id}.density.bed.starch"
  """
  bam2bed -d \
    < "${bam}" \
    | cut -f1-6 \
    | awk '{ if( \$6=="+" ){ s=\$2; e=\$2+1 } else { s=\$3-1; e=\$3 } print \$1 "\t" s "\t" e "\tid\t" 1 }' \
    | sort-bed - \
    > sample.bed
  
  unstarch "${params.density_buckets}" \
    | bedmap --faster --echo --count --delim "\t" - sample.bed \
    | awk -v OFS="\t" \
        -v binI=${params.density_step_size} \
        -v win="${params.density_window_width}" \
        'BEGIN{ halfBin=binI/2; shiftFactor=win-halfBin } { 
            print \$1,\$2+shiftFactor,\$3-shiftFactor,"id",\$4}' \
    | starch - \
    > "${name}"
    
  unstarch "${name}" \
    | awk -v binI="${params.density_step_size}" \
        -f "${moduleDir}/bin/bedToWig.awk" > density.wig
  wigToBigWig -clip density.wig "${params.chrom_sizes}" ${sample_id}.density.bw
  
  unstarch "${name}" | bgzip -c > ${sample_id}.density.bed.bgz
  tabix -p bed ${sample_id}.density.bed.bgz
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
  containerOptions fastaContainer

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
      paired: it[3]
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
    aligned_files = alignBwa(trimmed_reads) | filter_and_sort
    filtered_bam_files = merge_bam(aligned_files.groupTuple()) 
    | mark_duplicates 
    | filter_nuclear

    is_paired_dict = trimmed_reads.map(it -> tuple(it[0], it[3])).distinct()
    
    bam_files = filtered_bam_files.join(is_paired_dict)
    insert_size(bam_files)
    macs2(bam_files)
    density_files(filtered_bam_files)
    convert_to_cram(filtered_bam_files)
  emit:
    convert_to_cram.out
}

workflow mergeFiles {
  basepath = "/net/seq/data2/projects/sabramov/ENCODE4/atac_aligning/output"
  bam_files = Channel.fromPath(params.samples_file)
    | splitCsv(header:true, sep:'\t')
    | map(row -> tuple(row.sample_id,
      row.align_id,
      file("${basepath}/${row.align_id}/${row.align_id}.filtered.cram"),
      file("${basepath}/${row.align_id}/${row.align_id}.filtered.cram.crai"),
     ))
    | filter { it[2].exists() }
    | map(it -> tuple(it[0], it[2]))
    | set_key_for_group_tuple
    | groupTuple()
    | merge_bam
    | mark_duplicates 
    | filter_nuclear
    | convert_to_cram

    // is_paired_dict = trimmed_reads.map(it -> tuple(it[0], it[3])).distinct()
    
    // bam_files = filtered_bam_files.join(is_paired_dict)
    // insert_size(bam_files)
    // macs2(bam_files)
    // density_files(filtered_bam_files)
  emit:
    bam_files
}


workflow {
    fastq_trimmed_paired = Channel
      .fromPath(params.samples_file)
      .splitCsv(header:true, sep:'\t')
		  .map(row -> tuple(row.sample_id, row.reads1, row.reads2, row.type == 'paired'))
    alignReads(set_key_for_group_tuple(fastq_trimmed_paired))
}


process convert_to_bam {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/stats"
    container "${params.container}"
    scratch true

    when:
        params.do_macs

    input:
        tuple val(sample_id), path(cram), path(crai), val(is_paired)

    output:
        tuple val(sample_id), path(bam), path("${bam}.bai"), val(is_paired)

    script:
    bam = cram.baseName + '.bam'
    """
    samtools view ${cram} -b -h > ${bam}
    samtools index ${bam}
    """
}
workflow doMacs2 {
  Channel.fromPath(params.samples_file)
      | splitCsv(header:true, sep:'\t')
      | map(row -> tuple(row.ag_id, row.bam_file, "${row.bam_file}.crai", row.type == 'paired'))
      | convert_to_bam
      | macs2
}