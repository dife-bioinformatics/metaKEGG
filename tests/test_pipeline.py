import os

from tests.test_config import (
    input_file_path_g, input_file_path_t, input_file_path_m,
    methylation_path, input_label_g, input_label_m,
    sheet_name_paths, sheet_name_genes, sheet_name_transcripts ,methylation_gene,
    methylation_pvalue, miRNA_pvalue, miRNA_gene, miRNA_path, methylation_pvalue_thresh, miRNA_pvalue_thresh,
    genes_column, log2fc_column, genes_sheet_name, pathways_sheet_name, input_file_path_bulk , analysis_type,
    save_to_eps, count_threshold, benjamini_threshold, output_folder_name, miRNA_path_quant,
    methylation_probe_column, miRNA_column, methylation_path_quant   
)
from src.metaKEGG.modules.pipeline import Pipeline

current_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_dir)


def test_single_input_genes():
    input_file_path = input_file_path_g
    input_label = input_label_g
    my_pipeline = Pipeline(input_file_path=input_file_path,
                           sheet_name_paths=sheet_name_paths,
                           sheet_name_genes=sheet_name_genes,
                           genes_column=genes_column,
                           log2fc_column=log2fc_column,
                           input_label=input_label,
                           save_to_eps=save_to_eps)

    my_pipeline.single_input_genes(benjamini_threshold=benjamini_threshold , count_threshold=count_threshold , pathway_pvalue=None)
    os.chdir(current_dir)

def test_single_input_transcripts():
    input_file_path = input_file_path_t
    input_label = input_label_g
    sheet_name_genes = sheet_name_transcripts
    my_pipeline = Pipeline(input_file_path=input_file_path,
                           sheet_name_paths=sheet_name_paths,
                           sheet_name_genes=sheet_name_genes,
                           genes_column=genes_column,
                           log2fc_column=log2fc_column,
                           input_label=input_label,
                           save_to_eps=save_to_eps)

    my_pipeline.single_input_transcripts(benjamini_threshold=benjamini_threshold , count_threshold=count_threshold)
    os.chdir(current_dir)

def test_multiple_inputs():
    os.chdir(current_dir)
    input_file_path = input_file_path_m
    input_label = input_label_m
    my_pipeline = Pipeline(input_file_path=input_file_path,
                           sheet_name_paths=sheet_name_paths,
                           sheet_name_genes=sheet_name_genes,
                           genes_column=genes_column,
                           log2fc_column=log2fc_column,
                           input_label=input_label,
                           save_to_eps=save_to_eps)

    my_pipeline.multiple_inputs(count_threshold=count_threshold, benjamini_threshold=benjamini_threshold)
    os.chdir(current_dir)

def test_single_input_with_methylation():
    os.chdir(current_dir)
    input_file_path = input_file_path_g
    input_label = input_label_g
    my_pipeline = Pipeline(
        input_file_path=input_file_path,
        sheet_name_paths=sheet_name_paths,
        sheet_name_genes=sheet_name_genes,
        genes_column=genes_column,
        log2fc_column=log2fc_column,
        input_label=input_label,
        save_to_eps=save_to_eps)

    my_pipeline.single_input_with_methylation(methylation_path=methylation_path, methylation_genes=methylation_gene, methylation_pvalue=methylation_pvalue,
                                              methylation_pvalue_thresh = methylation_pvalue_thresh,
                                              count_threshold=count_threshold, benjamini_threshold=benjamini_threshold)
    os.chdir(current_dir)

def test_single_input_with_miRNA():
    os.chdir(current_dir)
    input_file_path = input_file_path_g
    input_label = input_label_g
    my_pipeline = Pipeline(
        input_file_path=input_file_path,
        sheet_name_paths=sheet_name_paths,
        sheet_name_genes=sheet_name_genes,
        input_label=input_label,
        genes_column=genes_column,
        log2fc_column=log2fc_column,
        save_to_eps=save_to_eps)

    my_pipeline.single_input_with_miRNA(miRNA_path=miRNA_path, miRNA_genes=miRNA_gene,
                                        miRNA_pvalue=miRNA_pvalue, miRNA_pvalue_thresh=miRNA_pvalue_thresh,
                                        count_threshold=count_threshold, benjamini_threshold=benjamini_threshold)
    os.chdir(current_dir)

def test_single_input_with_methylation_and_miRNA():
    os.chdir(current_dir)
    input_file_path = input_file_path_g
    input_label = input_label_g

    my_pipeline = Pipeline(
        input_file_path=input_file_path,
        sheet_name_paths=sheet_name_paths,
        sheet_name_genes=sheet_name_genes,
        genes_column=genes_column,
        log2fc_column=log2fc_column,
        input_label=input_label,
        save_to_eps=save_to_eps)

    my_pipeline.single_input_with_methylation_and_miRNA(methylation_path=methylation_path, methylation_genes=methylation_gene,
                                                        methylation_pvalue=methylation_pvalue, methylation_pvalue_thresh=methylation_pvalue_thresh,
                                                         miRNA_path=miRNA_path, miRNA_genes=miRNA_gene,
                                                         miRNA_pvalue=miRNA_pvalue, miRNA_pvalue_thresh=miRNA_pvalue_thresh,
                                                         count_threshold=count_threshold, benjamini_threshold=benjamini_threshold)
    os.chdir(current_dir)

def test_single_input_bulk():
    os.chdir(current_dir)
    input_file_path = input_file_path_bulk
    input_label = None
    sheet_name_paths = pathways_sheet_name
    sheet_name_genes = genes_sheet_name
    my_pipeline = Pipeline(
        input_file_path=input_file_path,
        sheet_name_paths=sheet_name_paths,
        sheet_name_genes=sheet_name_genes,
        input_label=input_label,
        genes_column = genes_column,
        log2fc_column=log2fc_column,
        save_to_eps=save_to_eps)

    my_pipeline.single_input_genes_bulk_mapping()
    os.chdir(current_dir)


def test_output_folder_scheme():
    input_file_path = input_file_path_g
    input_label = input_label_g
    my_pipeline = Pipeline(input_file_path=input_file_path,
                           sheet_name_paths=sheet_name_paths,
                           sheet_name_genes=sheet_name_genes,
                           input_label=input_label,
                           save_to_eps=save_to_eps,
                           genes_column=genes_column,
                           log2fc_column=log2fc_column,
                           output_folder_name=output_folder_name,
                           folder_extension='with_extension')

    my_pipeline.single_input_genes(count_threshold=count_threshold, benjamini_threshold=benjamini_threshold,pathway_pvalue=None)
    os.chdir(current_dir)


    
def test_single_input_with_miRNA_quantification():
    os.chdir(current_dir)
    input_file_path = input_file_path_g
    input_label = input_label_g
    my_pipeline = Pipeline(
        input_file_path=input_file_path,
        sheet_name_paths=sheet_name_paths,
        sheet_name_genes=sheet_name_genes,
        genes_column=genes_column,
        log2fc_column=log2fc_column,
        input_label=input_label,
        save_to_eps=save_to_eps)

    my_pipeline.single_input_with_miRNA_quantification(miRNA_path=miRNA_path_quant, miRNA_genes=miRNA_gene,
                                                        miRNA_pvalue=miRNA_pvalue, miRNA_pvalue_thresh=miRNA_pvalue_thresh,
                                                        miRNA_ID_column=miRNA_column, 
                                                        count_threshold=count_threshold, benjamini_threshold=benjamini_threshold)
    os.chdir(current_dir)
    
def test_single_input_with_methylation_quantification():
    os.chdir(current_dir)
    input_file_path = input_file_path_g
    input_label = input_label_g
    my_pipeline = Pipeline(
        input_file_path=input_file_path,
        sheet_name_paths=sheet_name_paths,
        sheet_name_genes=sheet_name_genes,
        genes_column=genes_column,
        log2fc_column=log2fc_column,
        input_label=input_label,
        save_to_eps=save_to_eps)

    my_pipeline.single_input_with_methylation_quantification(methylation_path=methylation_path_quant, methylation_genes=methylation_gene,
                                                             methylation_pvalue=methylation_pvalue, methylation_pvalue_thresh=methylation_pvalue_thresh,
                                                            methylation_probe_column=methylation_probe_column,probes_to_cgs=False,
                                                             count_threshold=count_threshold, benjamini_threshold=benjamini_threshold)
    os.chdir(current_dir)


def test_single_input_with_methylation_quantification_correct_probes():
    os.chdir(current_dir)
    input_file_path = input_file_path_g
    input_label = input_label_g
    my_pipeline = Pipeline(
        input_file_path=input_file_path,
        sheet_name_paths=sheet_name_paths,
        sheet_name_genes=sheet_name_genes,
        input_label=input_label,
        genes_column=genes_column,
        log2fc_column=log2fc_column,
        save_to_eps=save_to_eps)

    my_pipeline.single_input_with_methylation_quantification(methylation_path=methylation_path_quant, methylation_genes=methylation_gene,
                                                             methylation_pvalue=methylation_pvalue, methylation_pvalue_thresh=methylation_pvalue_thresh,
                                                            methylation_probe_column=methylation_probe_column,probes_to_cgs=True,
                                                             count_threshold=count_threshold, benjamini_threshold=benjamini_threshold)
    os.chdir(current_dir)

def test_single_input_genes_with_compounds():
    input_label = input_label_g
    my_pipeline = Pipeline(input_file_path=input_file_path_g,
                           sheet_name_paths=sheet_name_paths,
                           sheet_name_genes=sheet_name_genes,
                           input_label=input_label,
                           log2fc_column=log2fc_column,
                           genes_column=genes_column,
                           save_to_eps=True,
                           compounds_list=['C00031' , 'C00162'] , folder_extension='compounds')

    my_pipeline.single_input_genes(count_threshold=1 , benjamini_threshold=benjamini_threshold)
    os.chdir(current_dir)