include { callHotspots } from "./hotspots_calling"
include { mark_duplicates; filter_nuclear; total_bam_stats; convert_to_cram } from "./aligning"
include { get_container } from "./helpers"


fastaContainer = get_container(params.genome_fasta_file)



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


process spot_score {

    // Impossible to use anywhere except Altius cluster
    //conda params.conda
    module "bedops/2.4.35-typical:samtools/1.3:modwt/1.0:kentutil/302:hotspot2/2.1.1:jdk/1.8.0_92:gcc/4.7.2:R/3.2.5:picard/2.8.1:git/2.3.3:coreutils/8.25:bedtools/2.25.0:python/3.5.1:pysam/0.9.0:htslib/1.6.0:numpy/1.11.0:atlas-lapack/3.10.2:scipy/1.0.0:scikit-learn/0.18.1:preseq:/2.0.3:gsl/2.4"
    publishDir "${params.outdir}/${ag_id}"


    input:
        tuple val(ag_id), path(bam_file), path(bam_file_index)

    output:
        tuple val(ag_id), path("r1.spot.out")

    script:
    renamed_input = "r1.${bam_file.extension}"
    genome_file = file(params.genome_fasta_file)
    genome_prefix = "${genome_file.parent}/${genome_file.simpleName}"

    """
    # workaround for hotspots1 naming scheme...
    # might need to check bam conversion
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

    starch --header \
        r1-both-passes/r1.hotspot.twopass.zscore.wig > r1.spots.starch

    bash $moduleDir/bin/info.sh \
      r1.spots.starch hotspot1 r1.spot.out \
      > r1.hotspot.info
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


process mark_dups_subsample {
    container "${params.container}"
    containerOptions "${fastaContainer}"
    tag "${ag_id}"
    cpus 2

    input:
        tuple val(ag_id), path(bam_file), path(bam_file_index)
    
    output:
        tuple val(ag_id), path(name), path("${name}.bai")

    script:
    name = "${ag_id}.marked.bam"
    """
    picard RevertSam \
        INPUT=${bam_file} \
        OUTPUT=clear.bam \
        RESTORE_HARDCLIPS=false \
        VALIDATION_STRINGENCY=SILENT \
        REMOVE_DUPLICATE_INFORMATION=true\
        SORT_ORDER=coordinate \
        RESTORE_ORIGINAL_QUALITIES=false \
        REMOVE_ALIGNMENT_INFORMATION=false

    picard MarkDuplicatesWithMateCigar \
        INPUT=clear.bam \
        METRICS_FILE=perc_dups.txt \
        OUTPUT=${name} \
        ASSUME_SORTED=true \
        MINIMUM_DISTANCE=300 \
        VALIDATION_STRINGENCY=SILENT \
        READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
    samtools index ${name}
    """
}

workflow spot1score {
    bams = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
            row.sample_id, 
            file(row.cram_file), 
            file("${row.cram_file}.crai")))
        | preprocessBams
        | spot_score
}


workflow subsampleToFrac {
    input_data = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.sample_id, file(row.cram_file), file(row.cram_index), row.frac))
    
    input_data
        | map(it -> tuple(it[0], it[1], it[2]))
        | filter_nuclear
        | join(input_data.map(it -> tuple(it[0], it[3])))
        | subsample_with_pairs_frac
        | mark_dups_subsample
        | (callHotspots & total_bam_stats & convert_to_cram)
}

workflow filterAndCallHotspots {
    Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.sample_id, file(row.bam_file), file(row.bam_index)))
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
            row.sample_id, 
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
            row.sample_id, 
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
