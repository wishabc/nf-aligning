# nf-aligning

A Nextflow pipeline for aligning DNase-seq/ATAC-seq style short-read data to a reference genome with BWA. Reads can come from local FASTQ files already on the server, or be pulled directly from the Sequence Read Archive (SRA) by accession.

The pipeline handles adapter trimming (`fastp`), alignment (`bwa`), and downstream hotspot calling, and can optionally bootstrap its input by downloading raw reads from SRA before running the same trim/align/hotspot-calling steps.

## Pipeline Steps

1. **Download** *(SRA mode only, `-entry alignFromSRA`)* — fetch raw reads for each `align_id` from the SRA using `prefetch` and `fasterq-dump`.
2. **trimReads** — trim adapter sequences from the reads with `fastp`, quality and length filtering disabled.
3. **alignReads** — align each run with `bwa aln`/`bwa sampe` (paired-end) or `bwa samse` (single-end), sort and filter the alignments, then merge all runs belonging to the same `sample_id` into a single BAM and mark duplicates. From there: basic alignment stats (filtered-aligned read count and duplicate count, via `samtools view` on the CRAM, counting SAM flag `1024` for duplicates) are computed, nuclear reads are filtered out for downstream analysis, insert-size metrics and MACS2 peak calling are run, and the filtered alignments are converted to CRAM.
4. **call_hotspots** — call hotspots on the aligned, nuclear-filtered reads.

### Method summary

Reads were trimmed of adapter sequences (`fastp --disable_quality_filtering --disable_length_filtering --detect_adapter_for_pe`) and aligned to either the GRCh38/hg38 or mm10 genome assembly using `bwa aln`.

## Requirements

- [Nextflow](https://www.nextflow.io/)
- A container runtime (Docker/Singularity/Apptainer) able to run the pipeline's container image, referenced via `params.container`
- A compute cluster/executor configured in `nextflow.config` — see [Compute environment](#compute-environment) below

## Compute environment

This pipeline does not require access to any specific cluster. `-profile Altius` is just the example profile used at Altius (SLURM-based); adjust or add a profile in `nextflow.config` to match your own cluster/executor.

See the [Nextflow executor documentation](https://www.nextflow.io/docs/latest/executor.html) and [config profiles documentation](https://www.nextflow.io/docs/latest/config.html#config-profiles) for how to define a profile for your environment.

A minimal example for a SLURM cluster:

```groovy
profiles {
    slurm {
        process.executor = 'slurm'
        process.queue = 'your-queue-name'
        singularity.enabled = true
        singularity.autoMounts = true
    }
}
```

Run with `-profile slurm` (or whatever you name it) in place of `-profile Altius`.

## Container

The pipeline runs entirely inside a single Apptainer/Singularity container built from `aligning-pipeline.def` — no conda environment needs to be installed on the host. The image bundles HotSpot1, HotSpot2, hotspot3, kentUtils, and all of their Python/R dependencies, so alignment, subsampling, and hotspot calling (including hotspot3) all run from the same container.

### Building the `.sif` image

From the repository root:

```bash
singularity build --fakeroot aligning-pipeline.sif aligning-pipeline.def
```

This takes roughly **~15 minutes** to build (it compiles HotSpot1/HotSpot2/kentUtils from source and resolves two separate conda environments along the way).

`--fakeroot` lets you build without root access; your cluster account needs subuid/subgid ranges configured for this to work (ask your cluster admin if the build fails with a fakeroot-related permission error). If you have actual root/sudo instead, use `sudo singularity build aligning-pipeline.sif aligning-pipeline.def`.

Once built, point `params.container` at the resulting `aligning-pipeline.sif` path. Every process that references `params.container` runs from that single image — no separate conda setup or module loads needed for subsampling, hotspot calling, or anything else downstream.

Rebuild the image any time `aligning-pipeline.def` or one of the `environment*.yml` files it references changes.

## Setup

1. Clone the repository and build the container image referenced by `params.container` (see [Container](#container)).
2. Fill in the required file paths in `params.config` (reference genome, chrom sizes, hotspot annotation files, etc. — see [Parameters](#parameters) below).
3. Configure a profile for your compute environment in `nextflow.config` (see [Compute environment](#compute-environment) above).
4. Prepare a tab-separated samples file and point `params.samples_file` at it (see [Input files](#input-files) below).

## Usage

### Align local FASTQ files

```bash
nextflow run main.nf -profile Altius
```

This runs the default workflow:

1. Reads `params.samples_file` and groups records by `align_id`.
2. Trims adapters with `fastp` (`trimReads`).
3. Aligns trimmed reads to the reference genome with BWA (`alignReads`).

### Align reads downloaded from SRA

Use this mode to fetch raw reads for one or more SRA run accessions and run them through the same trimming/alignment/hotspot-calling pipeline, without needing to stage FASTQs on disk first.

```bash
nextflow run main.nf -profile Altius -entry alignFromSRA
```

This runs the `alignFromSRA` workflow:

1. Reads `params.samples_file` and, for each unique `align_id`, downloads the run with `prefetch` and converts it to FASTQ with `fasterq-dump` (paired- and single-end output are both attempted; whichever exists is used).
2. Also pulls run metadata via `ffq` and publishes it alongside the FASTQs under `<outdir>/<sample_id>/stats`.
3. Detects whether the run is paired- or single-end based on which FASTQ files are non-empty, and drops any run for which no reads were downloaded (e.g. failed/invalid accessions — download errors are ignored rather than failing the pipeline).
4. Trims, aligns, and calls hotspots exactly as in the default workflow.

Downloaded FASTQs are published to `<outdir>/<sample_id>/`, so re-running the pipeline can reuse previously downloaded data instead of re-fetching it from SRA.

## Input files

`params.samples_file` is a **tab-separated** file with a header row. The columns required depend on the entry point:

### Default workflow (local FASTQs)

| Column | Description |
|---|---|
| `sample_id` | Sample/group identifier used for output paths |
| `align_id` | Unique identifier for this alignment (used for grouping and file naming) |
| `reads1` | Path to read 1 FASTQ |
| `reads2` | Path to read 2 FASTQ (required only when `type` is `paired`) |
| `type` | `paired` or `single` |
| `adapterP7` | (Optional) Adapter sequence for read 1; auto-detected by `fastp` if omitted |
| `adapterP5` | (Optional) Adapter sequence for read 2, paired-end only; auto-detected if omitted |

Adapter sequences are sanitized automatically — any trailing ambiguous (non-ACGT) bases are stripped, with a warning printed to the console.

Example:

```tsv
sample_id align_id reads1 reads2 type adapterP7 adapterP5
LIB001 AG001 /data/reads/AG001_R1.fastq.gz /data/reads/AG001_R2.fastq.gz paired
LIB002 AG002 /data/reads/AG002.fastq.gz single
```

### `alignFromSRA` workflow

| Column | Description |
|---|---|
| `sample_id` | Sample/group identifier used for output paths |
| `align_id` | SRA run accession (e.g. `SRR12345678`) to download |

Example:

```tsv
sample_id align_id
SRX22512992 SRR26815981
SRX22512992 SRR26815982
```

## Parameters

Key parameters set in `params.config` (paths below are examples from a human GRCh38 configuration and should be adjusted for your genome/environment):

| Parameter | Description |
|---|---|
| `outdir` | Output directory (relative path recommended) |
| `samples_file` | Path to the tab-separated samples file |
| `container` | Path to the pipeline's container image |
| `threads` | CPU threads for read download/trimming steps |
| `genome_fasta_file` | Reference genome FASTA (indexed for BWA) |
| `nuclear_chroms` / `nuclear_chrom_sizes` | Nuclear chromosome name/size lists, used to restrict analysis to nuclear DNA |
| `chrom_sizes` / `chrom_sizes_bed` | Genome-wide chromosome sizes |
| `do_macs` | Whether to run MACS2 peak calling |
| `save_cram_mode` | Which reads to retain in output CRAMs (`all` or `nuclear`) |
| `fdrs` | Comma-separated FDR thresholds for hotspot calling |
| `save_debug` | Whether to keep intermediate debug output |
| `mappable_bases` / `mappable` | Mappability annotation files used in hotspot calling |
| `centers` | Center-site annotation file for hotspot calling |
| `density_buckets` | Chromatin density bucket file |
| `chrominfo` | Chromosome info file for the SPOT1 score step |
| `hotspots_dir` | Path to the HotSpot1 distribution |
| `readlength` | Read length used by HotSpot1 |
| `subsampling_spot1_depth` / `subsample_hotspot2_depth` | Read-depth targets used when subsampling for QC metrics |

## Output

Results are published under `<outdir>/<sample_id>/`. Raw FASTQs are nested under a subdirectory named by `align_id` (the SRA run accession in SRA mode); all downstream output files are prefixed with `sample_id` instead:

```
<outdir>/<sample_id>/
├── <align_id>/                                          # raw run FASTQs (SRA mode), named by SRA run accession, e.g. SRR13520417/
├── <sample_id>.background.bw                            # local background signal (bigWig)
├── <sample_id>.cutcounts.bed.gz(.tbi)                    # per-base cut counts, tabix-indexed
├── <sample_id>.fit_stats.tsv.gz(.tbi, .iter1, .iter1.tbi)  # background model fit statistics
├── <sample_id>.normalized_density.bw                     # normalized signal density (bigWig)
├── <sample_id>.peak_calling.log                          # peak calling log
├── <sample_id>.per_segment.background.bw                 # per-segment background signal
├── <sample_id>.pvals.parquet/                            # per-base p-values (partitioned parquet)
├── <sample_id>.thresholds.bw                             # significance thresholds (bigWig)
├── <sample_id>.total_cutcounts                           # total cut count
├── <sample_id>.total_sequencing_stats.txt                # sequencing depth/stats summary
├── debug/                                                # intermediate/debug output (only when `save_debug = true`)
├── fdr0.001/                                             # hotspot calls thresholded at FDR 0.001
├── fdr0.01/                                              # hotspot calls thresholded at FDR 0.01
├── fdr0.05/                                              # hotspot calls thresholded at FDR 0.05
└── stats/
    ├── <align_id>/                                       # per-run stats subdirectory
    ├── <align_id>_info.json                              # SRA run metadata (from `ffq`)
    ├── <sample_id>.CollectInsertSizeMetrics.picard(.pdf)  # Picard insert size metrics
    ├── <sample_id>.MarkDuplicates.picard                  # Picard duplicate-marking metrics
    └── macs2.<sample_id>_*                                # MACS2 peak calling output (narrowPeak, summits, pileup/lambda bedgraphs, xls)
```

The set of `fdr*` directories is controlled by `params.fdrs` (comma-separated FDR thresholds).

## Notes

- SRA downloads use `errorStrategy 'ignore'`: if `prefetch`/`fasterq-dump` fails for a given accession, the pipeline continues rather than aborting, and that run is silently dropped from downstream processing.
- `align_id` must be unique across the samples file — records are deduplicated on this column before processing.
- Hotspot calling can be re-run directly from precomputed cutcounts (skipping trimming and alignment entirely) via the `fromCutcounts` workflow in `hotspots_calling.nf` — useful for iterating on hotspot-calling parameters without re-aligning.

## TODO

- [ ] Document the stats and subsampling scripts (`spot1score`, `subsampleToFrac`, `filterAndCallHotspots`, `preprocessBams`) once cleaned up.
