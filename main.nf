#!/usr/bin/env nextflow
include { alignReads } from "./aligning"
include { callHotspots } from "./hotspots_calling"
include { trimReadsFromFile } from "./trimming"

workflow fromFiles {
    trimReadsFromFile() | alignReads | callHotspots
}

// workflow fromSRA {
//     downloadFiles | processReadsFromFiles
// }

workflow {
    fromFiles()
}
