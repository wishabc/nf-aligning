#!/usr/bin/env nextflow
include { alignReads } from "./aligning"
include { callHotspots } from "./hotspots_calling"
include { trimReadsFromFile } from "./trimming"

workflow alignFromFiles {
    trimReadsFromFile() | alignReads | callHotspots
}

// workflow alignFromSRA {
//     downloadFiles | alignFromFiles
// }

workflow {
    alignFromFiles()
}
