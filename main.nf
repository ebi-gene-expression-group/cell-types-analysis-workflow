#!/usr/bin/env nextflow 

PRED_LABELS_DIR = Channel.fromPath(params.input_dir).first()
REF_LABELS_FILE = Channel.fromPath(params.ref_labels_file).first()
ONTOLOGY_GRAPH = Channel.fromPath(params.ontology_graph).first()

process remove_tech_duplicates {
    publishDir "${params.results_dir}", mode: 'copy'
    conda "${baseDir}/envs/cell_types_analysis.yaml"
    input:
        file(ref_labels_file) from REF_LABELS_FILE
    output:
        file("metadata_filtered.tsv") into REF_LABELS_FILTERED
	file('metadata_filt') into REF_LABELS_FILTERED_DIR
    """
    remove_tech_duplicates.R\
                --input-metadata ${ref_labels_file}\
                --barcode-col-ref ${params.barcode_col_ref}\
                --output-file metadata_filtered.tsv
    
    mkdir -p metadata_filt/
    cp metadata_filtered.tsv metadata_filt 
    """
}

process build_cell_ontology_dict{
    publishDir "${params.results_dir}", mode: 'copy'
    conda "${baseDir}/envs/cell_types_analysis.yaml"
    input:
	file('metadata_filt') from REF_LABELS_FILTERED_DIR
    output: 
        file("${params.ontology_dict}") into ONTOLOGY_DICT

    """
    build_cell_ontology_dict.R\
    	--input-dir ${metadata_filt}\
	--condensed-sdrf ${params.condensed_sdrf}\
	--barcode-col-name ${params.barcode_col_ref}\
	--cell-label-col-name ${params.label_column_ref}\
	--cell-ontology-col-name ${params.cell_ontology_col}\
	--output-dict-path ${params.ontology_dict}\
	--output-text-path ${params.ontology_table}
    """
}

process get_tool_performance_table {
    publishDir "${params.results_dir}", mode: 'copy'
    conda "${baseDir}/envs/cell_types_analysis.yaml"
    input:
        file(pred_labels_dir) from PRED_LABELS_DIR
        file(ref_labels_file) from REF_LABELS_FILTERED
        file(cl_dictionary) from ONTOLOGY_DICT
	file(ontology_graph) from ONTOLOGY_GRAPH
    output:
        file("${params.tool_perf_table}") into TOOL_PERF_TABLE

    """
    get_tool_performance_table.R\
                 --input-dir ${pred_labels_dir}\
                 --ref-file ${ref_labels_file}\
                 --lab-cl-mapping ${cl_dictionary}\
                 --ontology-graph ${ontology_graph}\
                 --barcode-col-ref ${params.barcode_col_ref}\
                 --barcode-col-pred ${params.barcode_col_pred}\
                 --label-column-ref ${params.label_column_ref}\
                 --label-column-pred ${params.label_column_pred}\
                 --semantic-sim-metric ${params.semantic_sim_metric}\
                 --output-path ${params.tool_perf_table}
    """
}

process generate_empirical_cdf {
    conda "${baseDir}/envs/cell_types_analysis.yaml"
    input:
        file(ref_labels_file) from REF_LABELS_FILTERED
	file(ontology_dict) from ONTOLOGY_DICT

    output: 
        file("${params.empirical_dist}") into EMP_DISTRIBUTION

    """
    get_empirical_dist.R\
            --input-ref-file ${ref_labels_file}\
            --label-column-ref ${params.label_column_ref}\
            --lab-cl-mapping ${ontology_dict}\
            --num-iterations ${params.num_iter}\
            --num-cores ${params.num_cores}\
	    --ontology-graph ${params.ontology_graph}\
            --semantic-sim-metric ${params.semantic_sim_metric}\
            --output-path ${params.empirical_dist}
    """
}

process get_pvals {
    publishDir "${params.results_dir}", mode: 'copy'
     conda "${baseDir}/envs/cell_types_analysis.yaml"
     input:
        file(tool_perf_table) from TOOL_PERF_TABLE
        file(emp_dist) from EMP_DISTRIBUTION

     output:
        file("${params.tool_table_pvals}") into TOOL_TABLE_PVALS

     """
     get_tool_pvals.R\
             --input-table ${tool_perf_table}\
             --emp-dist-list ${emp_dist}\
             --output-table ${params.tool_table_pvals}
     """
}
