include { callHotspots; spot_score } from "./hotspots_calling"
include { percent_dup } from "./stats"


params.conda = "/home/sabramov/miniconda3/envs/jupyterlab"


process take_r1_from_pair {
    tag "${ag_id}"
    conda params.conda
    cpus 3

    input:
        tuple val(ag_id), path(bam_file), path(bam_file_index)

    output:
        tuple val(ag_id), path(name), path("${name}.bai")
    
    script:
    name = "${ag_id}.r1.bam"
    """
    samtools view -F 140 -u -h \
        --reference ${params.genome_fasta_file} ${bam_file} \
        | samtools sort -@${task.cpus} -O bam - > ${name}
	samtools index ${name}
    """
}


process subsample {
    tag "${ag_id}"
    conda params.conda
    publishDir "${params.outdir}/${ag_id}"

    input:
        tuple val(ag_id), path(bam_file), path(bam_file_index)
    
    output:
        tuple val(ag_id), path(name), path("${name}.bai")

    script:
    name = "${ag_id}.subsampled.bam"
    """
    python3 $moduleDir/bin/random_sample.py ${bam_file} ${name} \
        `samtools view -c ${bam_file}` ${params.subsampling_spot1_depth}
    samtools index ${name}
    """
}


process filter_nuclear {
  conda params.conda
  tag "${ag_id}"
  scratch true

  input:
    tuple val(ag_id), path(bam), path(bam_index)

  output:
    tuple val(ag_id), path("${name}"), path("${name}.bai")

  script:
  name = "${ag_id}.filtered.bam"
  """
  cat "${params.nuclear_chroms}" \
    | xargs samtools view --reference ${params.genome_fasta_file} -F 4 -b ${bam} > ${name}

  samtools index ${name}
  """
}


workflow preprocessBams {
    take:
        data
    main:
        out = data
            | take_r1_from_pair
            | filter_nuclear
            | subsample
    emit:
        out
}


process subsample_with_pairs {
    conda params.conda
    tag "${ag_id}"
    cpus 2
    scratch true
    publishDir "${params.outdir}/${ag_id}"

    input:
        tuple val(ag_id), path(cram_file), path(cram_file_index)
    
    output:
        tuple val(ag_id), path(name), path("${name}.bai"), path(info)

    script:
    name = "${ag_id}.subsampled_pairs.bam"
    info = "${ag_id}.subsampling_info.txt"
    """
    total_reads=\$(samtools view -c "${cram_file}")
    frac=\$(echo "scale=4; if (\$total_reads <= ${params.subsample_hotspot2_depth}) 1 else ${params.subsample_hotspot2_depth}/\$total_reads" | bc)
    
    echo -e "Total_reads\tsampling_fraction" > ${info}
    echo -e "\$total_reads\t\$frac" >> ${info}
    samtools view ${cram_file} -h \
        --subsample-seed 42 \
        --reference ${params.genome_fasta_file} \
        --subsample \$frac \
        | samtools sort -@${task.cpus} > ${name}
    samtools index ${name}
    """
}

process subsample_with_pairs_frac {
    conda params.conda
    tag "${ag_id}"
    cpus 2
    scratch true
    publishDir "${params.outdir}/${ag_id}"

    input:
        tuple val(ag_id), path(cram_file), path(cram_file_index), val(frac)
    
    output:
        tuple val(ag_id), path(name), path("${name}.bai")

    script:
    name = "${ag_id}.subsampled_pairs.bam"
    """
    samtools view ${cram_file} -h \
        --subsample-seed 42 \
        --reference ${params.genome_fasta_file} \
        --subsample ${frac} \
        | samtools sort -@${task.cpus} > ${name}
    samtools index ${name}
    """
}

workflow spot1score {
    bams = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
            row.ag_id, 
            file(row.cram_file), 
            file("${row.cram_file}.crai")))
        | preprocessBams
        | spot_score
}


workflow subsampleTest {
    Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.ag_id, file(row.cram_file), file(row.cram_index)))
        | filter_nuclear
        | subsample_with_pairs
        | map(it -> tuple(it[0], it[1], it[2]))
        | callHotspots

    // bams 
    //     | preprocessBams
    //     | spot_score
}

workflow subsampleTest2 {
    input_data = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.ag_id, file(row.cram_file), file(row.cram_index), row.frac))
    
    input_data
        | map(it -> tuple(it[0], it[1], it[2]))
        | filter_nuclear
        | join(input_data.map(it -> tuple(it[0], it[3])))
        | subsample_with_pairs_frac
        | percent_dup

    // bams 
    //     | preprocessBams
    //     | spot_score
}

workflow filterAndCallHotspots {
    Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.ag_id, file(row.bam_file), file(row.bam_index)))
        | filter_nuclear
        | callHotspots
}


// DEFUNC 
def get_density_path(file_path, suffix=".subsampled_pairs.bam") {
    file_path.replace(suffix, ".density.bw")
}
workflow normalizeDensity {
    // FIXME: add density file path to samples_file
    Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
            row.ag_id, 
            file(get_density_path(row.cram_file, ".subsampled_pairs.bam")),
            file(row.cram_file),
            file(row.cram_index)
            )
        )
        | normalize_density
}

workflow tmp {
    metadata = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
            row.ag_id, 
            file(get_density_path(get_density_path(row.cram_file, ".filtered.cram"), ".cram")),
            file(row.cram_file),
            file(row.cram_index)
            )
        )
        | normalize_density
}

process normalize_density {
    tag "${ag_id}"
    publishDir "${params.outdir}/${ag_id}"
    conda params.conda

    input:
        tuple val(ag_id), path(density_starch), path(filtered_bam), path(bam_file_index)
    
    output:
        tuple val(ag_id), path(name)
    
    script:
    bin_size = 20
    scale = 1_000_000
    name = "${ag_id}.normalized.density.bw"
    """
    bigWigToBedGraph ${density_starch} tmp.bedGraph

    cat tmp.bedGraph \
        | awk \
            -v allcounts=\$(samtools view -c ${filtered_bam}) \
            -v extranuclear_counts=\$(samtools view -c "${filtered_bam}" chrM chrC) \
            -v scale=${scale} \
            -v OFS='\t' \
                'BEGIN{ tagcount=allcounts-extranuclear_counts } \
                { print \$1,\$2,\$3,"${ag_id}",(\$4/tagcount)*scale }' > tmp.bed
    cat tmp.bed | awk \
            -v "binI=${bin_size}" \
            -f "$moduleDir/bin/bedToWig.awk" \
        > tmp.wig

    wigToBigWig -clip tmp.wig ${params.chrom_sizes} ${name}
    """  
}
