include { callHotspots } from "./hotspots_calling"

params.conda = "/home/sabramov/miniconda3/envs/babachi"


process take_r1_from_pair {
    tag "${uniq_id}"
    conda params.conda
    cpus 3

    input:
        tuple val(uniq_id), path(bam_file), path(bam_file_index)

    output:
        tuple val(uniq_id), path(name), path("${name}.bai")
    
    script:
    name = "${uniq_id}.r1.bam"
    """
    samtools view -F 140 -u -h \
        --reference ${params.genome_fasta_file} ${bam_file} \
        | samtools sort -@${task.cpus} -O bam - > ${name}
	samtools index ${name}
    """
    
}

process remove_duplicates {
    tag "${uniq_id}"
    conda params.conda

    input:
        tuple val(uniq_id), path(bam_file), path(bam_file_index)
    
    output:
        tuple val(new_id), path(name), path("${name}.bai")

    script:
    new_id = "${uniq_id}.dedupped"
    name = "${new_id}.bam"
    """
    samtools view -F 1024 -h -b ${bam_file} > ${name}
    samtools index ${name}
    """
}

process subsample {
    tag "${uniq_id}"
    conda params.conda
    publishDir "${params.outdir}/${uniq_id}"

    input:
        tuple val(uniq_id), path(bam_file), path(bam_file_index)
    
    output:
        tuple val(uniq_id), path(name), path("${name}.bai")

    script:
    name = "${uniq_id}.subsampled.bam"
    """
    python3 $moduleDir/bin/random_sample.py ${bam_file} ${name} \
        `samtools view -c ${bam_file}` ${params.subsampled_count}
    samtools index ${name}
    """
}

genome_file = file(params.genome_fasta_file)
genome_prefix = "${genome_file.parent}/${genome_file.simpleName}"

process spot_score {

    // Impossible to use anywhere except Altius cluster
    //conda params.conda
    module "bedops/2.4.35-typical:samtools/1.3:modwt/1.0:kentutil/302:hotspot2/2.1.1:jdk/1.8.0_92:gcc/4.7.2:R/3.2.5:picard/2.8.1:git/2.3.3:coreutils/8.25:bedtools/2.25.0:python/3.5.1:pysam/0.9.0:htslib/1.6.0:numpy/1.11.0:atlas-lapack/3.10.2:scipy/1.0.0:scikit-learn/0.18.1:preseq:/2.0.3:gsl/2.4"
    publishDir "${params.outdir}/${uniq_id}"


    input:
        tuple val(uniq_id), path(bam_file), path(bam_file_index)

    output:
        tuple val(uniq_id), path("r1.*")

    script:
    renamed_input = "r1.${bam_file.extension}"
    """
    # workaround for hotspots1 naming scheme...
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

    starch --header r1-both-passes/r1.hotspot.twopass.zscore.wig \
      > r1.spots.starch

    bash $moduleDir/bin/info.sh \
      r1.spots.starch hotspot1 r1.spot.out \
      > r1.hotspot.info
    """
}


process filter_nuclear {
  conda params.conda
  tag "${uniq_id}"
  scratch true

  input:
    tuple val(uniq_id), path(bam), path(bam_index)

  output:
    tuple val(uniq_id), path("${name}"), path("${name}.bai")

  script:
  name = "${uniq_id}.filtered.bam"
  """
  samtools view -b -F 512 ${bam} > filtered.bam
  samtools index filtered.bam
  cat "${params.nuclear_chroms}" \
  | xargs samtools view -b filtered.bam > ${name}

  samtools index ${name}
  """
}


workflow preprocessBams {
    take:
        data
    main:
        r1_data = take_r1_from_pair(data) | filter_nuclear
        out = remove_duplicates(r1_data) 
            | mix(r1_data)
            | subsample
        
    emit:
        out
}


workflow {
    bams = Channel.fromPath(params.samples_file)
			| splitCsv(header:true, sep:'\t')
			| map(row -> tuple(
                row.ag_id, 
                file(row.bam_file), 
                file("${row.bam_file}.crai")))
    data = preprocessBams(bams) 
    //callHotspots(data) // For hotspot2
    spot_score(data) // make it work for hotspot1

}

process percent_dup {
    scratch true
    conda params.conda
    tag "${ag_id}"
    cpus 2
    publishDir "${params.outdir}/${ag_id}"

    input:
        tuple val(ag_id), path(bam_file), path(bam_file_index)
    
    output:
        tuple val(ag_id), path(name)

    script:
    name = "${ag_id}.spotdups.txt"
    """
    samtools sort -@${task.cpus} ${bam_file} > sorted.bam
    samtools index sorted.bam
    picard RevertSam \
      INPUT=sorted.bam \
      OUTPUT=clear.bam \
      VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATE_INFORMATION=true SORT_ORDER=coordinate \
      RESTORE_ORIGINAL_QUALITIES=false REMOVE_ALIGNMENT_INFORMATION=false
    picard MarkDuplicatesWithMateCigar \
      INPUT=clear.bam \
      METRICS_FILE=${name} \
      OUTPUT=/dev/null \
      ASSUME_SORTED=true \
      MINIMUM_DISTANCE=300 \
      VALIDATION_STRINGENCY=SILENT \
      READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
    """
}

// nextflow /script.nf -entry percentDup -profile Altius --samples_file <>
workflow percentDup {
    Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
            row.ag_id, 
            file(row.bam_file), 
            file("${row.bam_file}.crai")))
        | percent_dup
    
}


// DEFUNC 

// Assume paired end
process subsample_with_pairs {
    conda params.conda
    tag "${ag_id}"
    cpus 2

    input:
        tuple val(ag_id), path(cram_file), path(cram_file_index), val(frac)
    
    output:
        tuple val(ag_id), path(name), path("${name}.bai")

    script:
    name = "${ag_id}.subsampled.bam"
    """
    samtools view ${cram_file} -h -b -s 42,${frac} \
        | samtools sort -@${task.cpus} > ${name}
    samtools index ${name}
    """
}

workflow subsampleTest {
    bams = Channel.fromPath("/net/seq/data2/projects/sabramov/SuperIndex/dnase_peak_density_analysis/alignments_for_unresolved.tsv")
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.ag_id, file(row.filtered_alignments_bam), file(row.bam_index), row.frac))
        | subsample_with_pairs
    callHotspots(bams)
    bams 
        | preprocessBams
        | spot_score
}
