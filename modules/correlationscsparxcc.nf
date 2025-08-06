process correlationscsparxcc {

    tag 'cross_correlation_SparXCC'
    cpus 4

    input:
        path otu_table
        path expression_file

    output:
        path("SparXCC_output_day_common_samples.txt"),     emit: sparxcc_day_correlation_matrix
        path("SparXCC_output_night_common_samples.txt"),   emit: sparxcc_night_correlation_matrix

    script:
    """
    Rscript ${projectDir}/../bin/cross_cor_SparXCC.r $otu_table $expression_file
    """
}
