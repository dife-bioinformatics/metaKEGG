import os

from tests.test_config import (
    input_file_path_g, input_file_path_t, input_file_path_m,
    methylation_path, input_label_g, input_label_m,
    sheet_name_paths, sheet_name_genes, sheet_name_transcripts ,methylation_gene,
    methylation_pvalue, miRNA_pvalue, miRNA_gene, miRNA_path, methylation_pvalue_thresh, miRNA_pvalue_thresh,
    genes_column, log2fc_column, genes_sheet_name, pathways_sheet_name, input_file_path_bulk , analysis_type, save_to_eps, count_threshold, benjamini_threshold, output_folder_name, metadata_path_quant
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
                           input_label=input_label,
                           analysis_type=analysis_type,
                           save_to_eps=save_to_eps,
                           count_threshold=count_threshold,
                           benjamini_threshold=benjamini_threshold)

    my_pipeline.single_input_genes()
    os.chdir(current_dir)

def test_single_input_transcripts():
    input_file_path = input_file_path_t
    input_label = input_label_g
    sheet_name_genes = sheet_name_transcripts
    my_pipeline = Pipeline(input_file_path=input_file_path,
                           sheet_name_paths=sheet_name_paths,
                           sheet_name_genes=sheet_name_genes,
                           input_label=input_label,
                           analysis_type=analysis_type,
                           save_to_eps=save_to_eps,
                           count_threshold=count_threshold,
                           benjamini_threshold=benjamini_threshold)

    my_pipeline.single_input_transcripts()
    os.chdir(current_dir)

def test_multiple_inputs():
    os.chdir(current_dir)
    input_file_path = input_file_path_m
    input_label = input_label_m
    my_pipeline = Pipeline(input_file_path=input_file_path,
                           sheet_name_paths=sheet_name_paths,
                           sheet_name_genes=sheet_name_genes,
                           input_label=input_label,
                           analysis_type=analysis_type,
                           save_to_eps=save_to_eps,
                           count_threshold=count_threshold,
                           benjamini_threshold=benjamini_threshold)

    my_pipeline.multiple_inputs()
    os.chdir(current_dir)

def test_single_input_with_methylation():
    os.chdir(current_dir)
    input_file_path = input_file_path_g
    input_label = input_label_g
    my_pipeline = Pipeline(
        input_file_path=input_file_path,
        sheet_name_paths=sheet_name_paths,
        sheet_name_genes=sheet_name_genes,
        input_label=input_label,
        methylation_path=methylation_path,
        methylation_genes=methylation_gene,
        methylation_pvalue=methylation_pvalue,
        analysis_type=analysis_type ,
        save_to_eps=save_to_eps,
        count_threshold=count_threshold,
        benjamini_threshold=benjamini_threshold,
        methylation_pvalue_thresh = methylation_pvalue_thresh
    )

    my_pipeline.single_input_with_methylation()
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
        miRNA_path=miRNA_path,
        miRNA_genes=miRNA_gene,
        miRNA_pvalue=miRNA_pvalue,
        analysis_type=analysis_type,
        save_to_eps=save_to_eps,
        count_threshold=count_threshold,
        benjamini_threshold=benjamini_threshold,
        miRNA_pvalue_thresh=miRNA_pvalue_thresh
    )

    my_pipeline.single_input_with_miRNA()
    os.chdir(current_dir)

def test_single_input_with_methylation_and_miRNA():
    os.chdir(current_dir)
    input_file_path = input_file_path_g
    input_label = input_label_g

    my_pipeline = Pipeline(
        input_file_path=input_file_path,
        sheet_name_paths=sheet_name_paths,
        sheet_name_genes=sheet_name_genes,
        input_label=input_label,
        methylation_path=methylation_path,
        methylation_genes=methylation_gene,
        methylation_pvalue=methylation_pvalue,
        miRNA_path=miRNA_path,
        miRNA_genes=miRNA_gene,
        miRNA_pvalue=miRNA_pvalue,
        analysis_type=analysis_type,
        save_to_eps=save_to_eps,
        count_threshold=count_threshold,
        benjamini_threshold=benjamini_threshold,
        methylation_pvalue_thresh=methylation_pvalue_thresh,
        miRNA_pvalue_thresh=miRNA_pvalue_thresh
    )

    my_pipeline.single_input_with_methylation_and_miRNA()
    os.chdir(current_dir)

def test_single_input_bulk():
    os.chdir(current_dir)
    input_file_path = input_file_path_bulk
    input_label = None
    sheet_name_paths = pathways_sheet_name
    sheet_name_genes = genes_sheet_name
    count_threshold = None
    benjamini_threshold = None
    my_pipeline = Pipeline(
        input_file_path=input_file_path,
        sheet_name_paths=sheet_name_paths,
        sheet_name_genes=sheet_name_genes,
        input_label=input_label,
        methylation_path=methylation_path,
        methylation_genes=methylation_gene,
        methylation_pvalue=methylation_pvalue,
        miRNA_path=miRNA_path,
        miRNA_genes=miRNA_gene,
        miRNA_pvalue=miRNA_pvalue,
        genes_column = genes_column,
        log2fc_column=log2fc_column,
        analysis_type=analysis_type,
        save_to_eps=save_to_eps,
        count_threshold=count_threshold,
        benjamini_threshold=benjamini_threshold
    )

    my_pipeline.single_input_genes_bulk_mapping()
    os.chdir(current_dir)


def test_output_folder_scheme():
    input_file_path = input_file_path_g
    input_label = input_label_g
    my_pipeline = Pipeline(input_file_path=input_file_path,
                           sheet_name_paths=sheet_name_paths,
                           sheet_name_genes=sheet_name_genes,
                           input_label=input_label,
                           analysis_type=analysis_type,
                           save_to_eps=save_to_eps,
                           count_threshold=count_threshold,
                           benjamini_threshold=benjamini_threshold,
                           output_folder_name=output_folder_name,
                           folder_extension='with_extension')

    my_pipeline.single_input_genes()
    os.chdir(current_dir)


    
def test_single_input_with_miRNA_quantification():
    os.chdir(current_dir)
    input_file_path = input_file_path_g
    input_label = input_label_g
    my_pipeline = Pipeline(
        input_file_path=input_file_path,
        sheet_name_paths=sheet_name_paths,
        sheet_name_genes=sheet_name_genes,
        input_label=input_label,
        miRNA_path=metadata_path_quant,
        miRNA_genes=miRNA_gene,
        miRNA_pvalue=miRNA_pvalue,
        analysis_type=analysis_type,
        save_to_eps=save_to_eps,
        count_threshold=count_threshold,
        benjamini_threshold=benjamini_threshold,
        miRNA_pvalue_thresh=miRNA_pvalue_thresh
    )

    my_pipeline.single_input_with_miRNA_quantification()
    os.chdir(current_dir)
    
def test_single_input_with_methylation_quantification():
    os.chdir(current_dir)
    input_file_path = input_file_path_g
    input_label = input_label_g
    my_pipeline = Pipeline(
        input_file_path=input_file_path,
        sheet_name_paths=sheet_name_paths,
        sheet_name_genes=sheet_name_genes,
        input_label=input_label,
        methylation_path=metadata_path_quant,
        methylation_genes=methylation_gene,
        methylation_pvalue=methylation_pvalue,
        analysis_type=analysis_type,
        save_to_eps=save_to_eps,
        count_threshold=count_threshold,
        benjamini_threshold=benjamini_threshold,
        miRNA_pvalue_thresh=miRNA_pvalue_thresh
    )

    my_pipeline.single_input_with_methylation_quantification()
    os.chdir(current_dir)