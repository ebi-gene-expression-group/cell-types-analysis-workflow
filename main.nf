#!/usr/bin/env nextflow 

PRED_LABELS_DIR = Channel.fromPath(params.input_dir).first()
REF_LABELS_FILE = Channel.fromPath(params.ref_labels_file).first()

process remove_tech_duplicates {
    publishDir "${params.results_dir}", mode: 'copy'
    conda "${baseDir}/envs/cell_types_analysis.yaml"
    input:
        file(ref_labels_file) from REF_LABELS_FILE
    output:
        file("metadata_filtered.tsv") into REF_LABELS_FILTERED

    """
    remove_tech_duplicates.R\
                --input-metadata ${ref_labels_file}\
                --barcode-col-ref ${params.barcode_col_ref}\
                --output-file metadata_filtered.tsv
    """
}


process get_tool_performance_table {
    publishDir "${params.results_dir}", mode: 'copy'
    conda "${baseDir}/envs/cell_types_analysis.yaml"
    input:
        file(pred_labels_dir) from PRED_LABELS_DIR
        file(ref_labels_file) from REF_LABELS_FILTERED

    output:
        file("${params.tool_perf_table}") into TOOL_PERF_TABLE

    """
    get_tool_performance_table.R\
                 --input-dir ${pred_labels_dir}\
                 --ref-file ${ref_labels_file}\
                 --ontology-graph ${params.ontology_graph}\
                 --cell-ontology-col ${params.cell_ontology_col}\
                 --barcode-col-ref ${params.barcode_col_ref}\
                 --label-column-ref ${params.label_column_ref}\
                 --semantic-sim-metric ${params.semantic_sim_metric}\
                 --output-path ${params.tool_perf_table}
    """
}

process generate_empirical_cdf {
    conda "${baseDir}/envs/cell_types_analysis.yaml"
    input:
        file(ref_labels_file) from REF_LABELS_FILTERED

    output: 
        file("${params.empirical_dist}") into EMP_DISTRIBUTION

    """
    get_empirical_dist.R\
            --input-ref-file ${ref_labels_file}\
            --num-iterations ${params.num_iter}\
            --label-column-ref ${params.label_column_ref}\
            --cell-ontology-col ${params.cell_ontology_col}\
            --semantic-sim-metric ${params.semantic_sim_metric}\
            --num-cores ${params.num_cores}\
            --ontology-graph ${params.ontology_graph}\
            --output-path ${params.empirical_dist}
    """
}

process get_pvals {
    publishDir "${params.results_dir}", mode: 'copy'
     conda "${baseDir}/envs/cell_types_analysis.yaml"
     input:
        file(tool_perf_table) from TOOL_PERF_TABLE
        file(emp_distr) from EMP_DISTRIBUTION

     output:
        file("${params.tool_table_pvals}") into TOOL_TABLE_PVALS

     """
     get_tool_pvals.R\
             --input-table ${tool_perf_table}\
             --emp-dist-list ${emp_distr}\
             --output-table ${params.tool_table_pvals}
     """
}

