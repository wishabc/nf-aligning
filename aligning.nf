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
  errorStrategy "ignore"

  input:
    tuple val(group_key), path(bamfiles)

  output:
    tuple val(group_key), file(name), path("${name}.bai")

  script:
  name = "${group_key}.bam"
  if (bamfiles.size() > 1) {
    """
    samtools merge ${name} ${bamfiles}
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
  maxForks 5
  time "2d"

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
  mode = is_paired ? 'BAMPE' : 'BAM'
  """
  gen_size=\$(cat ${params.chrom_sizes} | awk '{print \$2}' | paste -sd+ | bc)
  macs2 callpeak \
    -t "${bam}" \
    -f ${mode} \
    -n macs2.${sample_id}\
    -g \$gen_size -p 0.01 \
    --shift 75 --extsize 150 \
    --nomodel -B --SPMR \
    --keep-dup all --call-summits 2>"macs2.${sample_id}.err" || echo "MACS2 failed" >> "macs2.${sample_id}.err"
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
    | awk -v binI=${params.density_step_size} -v win="${params.density_window_width}" \
        'BEGIN{ halfBin=binI/2; shiftFactor=win-halfBin } { print \$1 "\t" \$2 + shiftFactor "\t" \$3-shiftFactor "\tid\t" i \$4}' \
    | starch - \
    > "${name}"
    
  unstarch "${name}" | awk -v binI="${params.density_step_size}" -f "${moduleDir}/bin/bedToWig.awk" > density.wig
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

process spot_score {
 scratch false
  input:
    tuple val(meta), path(bam), path(bai), path(mappable_only), path(chrom_info)

  output:
    tuple val(meta), file('subsample.r1.spot.out'), file('spotdups.txt')

  script:
  read_length = (mappable_only.name =~ /K([0-9]+)/)[0][1]
  // conditional to deal with single end data
  view_params = "-F 12 -f 3"
  if (params.paired)  
    """
    # random sample
    samtools view -h ${view_params} "${bam}" \
      | awk '{if( ! index(\$3, "chrM") && \$3 != "chrC" && \$3 != "random"){print}}' \
      | samtools view -1 - \
      -o sampled.bam

    bash $moduleDir/bin/random_sample.sh sampled.bam subsample.bam 5000000
    samtools view -1 -f 0x0040 subsample.bam -o subsample.r1.bam
    # hotspot
    bash runhotspot.bash \
      \$HOTSPOT_DIR \
      ./ \
      subsample.r1.bam \
      "${genome_name}" \
      "${read_length}" \
      DNaseI
    # Remove existing duplication marks
    picard RevertSam \
      INPUT=subsample.bam \
      OUTPUT=clear.bam \
      VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATE_INFORMATION=true SORT_ORDER=coordinate \
      RESTORE_ORIGINAL_QUALITIES=false REMOVE_ALIGNMENT_INFORMATION=false
    picard MarkDuplicatesWithMateCigar \
      INPUT=clear.bam \
      METRICS_FILE=spotdups.txt \
      OUTPUT=/dev/null \
      ASSUME_SORTED=true \
      MINIMUM_DISTANCE=300 \
      VALIDATION_STRINGENCY=SILENT \
      READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
    """
    else
    """
    # random sample
    samtools view -h -F 12 -s 0.6 "${bam}" \
      | awk '{if( ! index(\$3, "chrM") && \$3 != "chrC" && \$3 != "random"){print}}' \
      | samtools view -1 - \
      -o paired.bam
    bash \$STAMPIPES/scripts/bam/random_sample.sh paired.bam subsample.bam 5000000
          samtools view -1 subsample.bam -o subsample.r1.bam
    # hotspot
    bash \$STAMPIPES/scripts/SPOT/runhotspot.bash \
      \$HOTSPOT_DIR \
      \$PWD \
      \$PWD/subsample.r1.bam \
      "${genome_name}" \
      "${read_length}" \
      DNaseI
    
    # Remove existing duplication marks
    picard RevertSam \
      INPUT=subsample.bam \
      OUTPUT=clear.bam \
      VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATE_INFORMATION=true SORT_ORDER=coordinate \
      RESTORE_ORIGINAL_QUALITIES=false REMOVE_ALIGNMENT_INFORMATION=false
    picard MarkDuplicates \
      INPUT=clear.bam \
      METRICS_FILE=spotdups.txt \
      OUTPUT=/dev/null \
      ASSUME_SORTED=true \
      VALIDATION_STRINGENCY=SILENT \
      READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*' 
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

workflow {
    fastq_trimmed_paired = Channel
      .fromPath(params.samples_file)
      .splitCsv(header:true, sep:'\t')
		  .map(row -> tuple( row.sample_id, row.reads1, row.reads2, row.type == 'paired'))
    alignReads(set_key_for_group_tuple(fastq_trimmed_paired))
}