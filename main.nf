#!/usr/bin/env nextflow
include { alignReads } from "./aligning"
include { callHotspots } from "./hotspots_calling"
include { trimReads } from "./trimming"


workflow {
    trimReads | alignReads | map(it -> tuple(it[0], it[1])) | callHotspots
}
