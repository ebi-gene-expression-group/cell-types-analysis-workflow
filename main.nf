#!/usr/bin/env nextflow 

PRED_LABELS_DIR = Channel.fromPath(params.input_dir).first()
REF_LABELS_FILE = Channel.fromPath(params.ref_labels_file).first()

process get_tool_performance_table {
    conda 'envs/{my_package}'
    input:
        file(pred_labels_dir) from PRED_LABELS_DIR
        file(ref_labels_file) from REF_LABELS_FILE

    output:
        file("${params.tool_perf_table}") into TOOL_PERF_TABLE

    """
    get_tool_performance_table.R\
                 --input-dir ${pred_labels_dir}\
                 --ref-file ${ref_labels_file}\
                 --output-path ${params.tool_perf_table}
    """
}


process generate_empirical_cdf {
    conda 'envs/{my_package}'
    input:
        file(ref_labels_file) from REF_LABELS_FILE

    output: 
        file("${params.empirical_dist}") into EMP_DISTRIBUTION

    """
    get_empirical_dist.R\
            --input-ref-file ${ref_labels_file}\
            --num-iterations ${params.num_iter}\
            --num-cores ${params.num_cores}\
            --output-path ${params.empirical_dist}
    """
}

process get_pvals {
     conda 'envs/{my_package}'
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

process get_per_cell_stats {
    conda 'envs/{my_package}'
    input:
        file(input_dir) from PRED_LABELS_DIR
        file(ref_labels_file) from REF_LABELS_FILE
        file(tool_table) from TOOL_TABLE_PVALS

    output:
        file("${params.cell_anno_table}") into CELL_ANNO_TABLE

    """
    get_cell_annotations_table.R\
        --input-dir ${input_dir}\
        --ref-file ${ref_labels_file}\
        --tool-table ${tool_table}\
        --output-path ${params.cell_anno_table}
    """
}





