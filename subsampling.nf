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
    python3 $moduleDir/bin/rmdup.py ${bam_file} ${name}
    samtools index ${name}
    """
}

process subsample {
    tag "${uniq_id}"
    conda params.conda
    publishDir "${params.outdir}/subsampled"

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

workflow preprocessBams {
    take:
        data
    main:
        r1_data = take_r1_from_pair(data)
        out = remove_duplicates(r1_data ) 
            | mix(r1_data)
            | subsample
        
    emit:
        out
}

workflow {
    bams = Channel.fromPath(params.samples_file)
			| splitCsv(header:true, sep:'\t')
			| map(row -> tuple(
                row.uniq_id, 
                file(row.bam_file), 
                file("${row.bam_file}.crai")))
    preprocessBams(bams) // | callHotspots

}