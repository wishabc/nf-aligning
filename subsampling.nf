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
        `samtools view -c ${bam_file}` 5000000
    samtools index ${name}
    """
}

params.readlength = 36
genome_file = file(params.genome_fasta_file)
genome_prefix = "${genome_file.parent}/${genome_file.simpleName}"
params.hotspots_dir = "/home/solexa/hotspot/"
process spot_score {

    // Impossible to use anywhere except Altius cluster
    conda params.conda
    module "modwt/1.0:kentutil/302:jdk/1.8.0_92:gcc/4.7.2:R/3.2.5:coreutils/8.25:atlas-lapack/3.10.2:preseq:/2.0.3:gsl/2.4"
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
    bash $moduleDir/bin/runhotspot.bash" \
      "${params.hotspots_dir}" \
      "./" \
      "${bam_file}" \
      "${genome_prefix}" \
      "${params.readlength}" \
      DNaseI

    starch --header r1-both-passes/r1.hotspot.twopass.zscore.wig \
      > r1.spots.starch

    bash $moduleDir/bin/info.sh \
      r1.spots.starch hotspot1 r1.spot.out \
      > r1.hotspot.info
    """
}


workflow preprocessBams {
    take:
        data
    main:
        r1_data = take_r1_from_pair(data)
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
    callHotspots(data) // For hotspot2
    spot_score(data) // make it work for hotspot1

}