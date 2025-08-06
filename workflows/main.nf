#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { mergeOtuTables }        from '../modules/mergeOtuTables.nf'
include { correlationMerged }     from '../modules/correlationMerged.nf'
include { correlationscsparxcc }   from '../modules/correlationscsparxcc.nf'

workflow {

    // Input channels
    Channel.fromPath(params.day_file)         .set { ch_day_file }
    Channel.fromPath(params.night_file)       .set { ch_night_file }
    Channel.fromPath(params.expression_file)  .set { ch_expression_file }


    // Run merge process
    ch_otu_table = mergeOtuTables(ch_day_file, ch_night_file)

    // Run correlationMerged
//    correlationMerged(ch_otu, ch_expression_file, ch_output_dir)

    // Run correlationscsparcc
    ch_scparxcc_mat = correlationscsparxcc(ch_otu_table, ch_expression_file)
}

