import os
import pytest

from tests.test_config import (
    input_file_path_g, input_file_path_t, input_file_path_m,
    methylation_file_path, input_label_m,
    sheet_name_paths, sheet_name_genes, sheet_name_transcripts ,methylation_genes_column,
    methylation_pvalue_column, miRNA_pvalue_column, miRNA_genes_column, miRNA_file_path, methylation_pvalue_threshold, miRNA_pvalue_threshold,
    genes_column, log2fc_column, genes_sheet_name, pathways_sheet_name, input_file_path_bulk ,
    save_to_eps, count_threshold, benjamini_threshold, output_folder_name, miRNA_path_quant,
    methylation_probe_column, miRNA_column, methylation_path_quant   
)
from src.metaKEGG.modules.pipeline_async import PipelineAsync

current_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_dir)

@pytest.mark.asyncio
async def test_gene_expression():
    input_file_path = input_file_path_g
    my_pipeline = PipelineAsync(input_file_path=input_file_path,
                           sheet_name_paths=sheet_name_paths,
                           sheet_name_genes=sheet_name_genes,
                           genes_column=genes_column,
                           log2fc_column=log2fc_column,
                           save_to_eps=save_to_eps)

    await my_pipeline.gene_expression()
    os.chdir(current_dir)

@pytest.mark.asyncio
async def test_transcript_expression():
    input_file_path = input_file_path_t
    sheet_name_genes = sheet_name_transcripts
    my_pipeline = PipelineAsync(input_file_path=input_file_path,
                           sheet_name_paths=sheet_name_paths,
                           sheet_name_genes=sheet_name_genes,
                           genes_column=genes_column,
                           log2fc_column=log2fc_column,
                           save_to_eps=save_to_eps)

    await my_pipeline.transcript_expression()
    os.chdir(current_dir)

@pytest.mark.asyncio
async def test_multiple_inputs():
    os.chdir(current_dir)
    input_file_path = input_file_path_m
    input_label = input_label_m
    my_pipeline = PipelineAsync(input_file_path=input_file_path,
                           sheet_name_paths=sheet_name_paths,
                           sheet_name_genes=sheet_name_genes,
                           genes_column=genes_column,
                           log2fc_column=log2fc_column,
                           save_to_eps=save_to_eps)

    await my_pipeline.multiple_inputs(input_label=input_label, count_threshold=count_threshold, benjamini_threshold=benjamini_threshold)
    os.chdir(current_dir)

@pytest.mark.asyncio
async def test_methylated_genes():
    os.chdir(current_dir)
    input_file_path = input_file_path_g
    my_pipeline = PipelineAsync(
        input_file_path=input_file_path,
        sheet_name_paths=sheet_name_paths,
        sheet_name_genes=sheet_name_genes,
        genes_column=genes_column,
        log2fc_column=log2fc_column,
        save_to_eps=save_to_eps)

    await my_pipeline.methylated_genes(methylation_file_path=methylation_file_path, methylation_genes_column=methylation_genes_column, methylation_pvalue_column=methylation_pvalue_column,
                                              methylation_pvalue_threshold=methylation_pvalue_threshold,
                                              count_threshold=count_threshold, benjamini_threshold=benjamini_threshold)
    os.chdir(current_dir)

@pytest.mark.asyncio
async def test_demirs_per_gene():
    os.chdir(current_dir)
    input_file_path = input_file_path_g
    my_pipeline = PipelineAsync(
        input_file_path=input_file_path,
        sheet_name_paths=sheet_name_paths,
        sheet_name_genes=sheet_name_genes,
        genes_column=genes_column,
        log2fc_column=log2fc_column,
        save_to_eps=save_to_eps)

    await my_pipeline.demirs_per_gene(miRNA_file_path=miRNA_file_path, miRNA_genes_column=miRNA_genes_column,
                                        miRNA_pvalue_column=miRNA_pvalue_column, miRNA_pvalue_threshold=miRNA_pvalue_threshold,
                                        count_threshold=count_threshold, benjamini_threshold=benjamini_threshold)
    os.chdir(current_dir)

@pytest.mark.asyncio
async def test_methylated_and_mirna_target_genes():
    os.chdir(current_dir)
    input_file_path = input_file_path_g

    my_pipeline = PipelineAsync(
        input_file_path=input_file_path,
        sheet_name_paths=sheet_name_paths,
        sheet_name_genes=sheet_name_genes,
        genes_column=genes_column,
        log2fc_column=log2fc_column,
        save_to_eps=save_to_eps)

    await my_pipeline.methylated_and_mirna_target_genes(methylation_file_path=methylation_file_path, methylation_genes_column=methylation_genes_column,
                                                        methylation_pvalue_column=methylation_pvalue_column, methylation_pvalue_threshold=methylation_pvalue_threshold,
                                                        miRNA_file_path=miRNA_file_path, miRNA_genes_column=miRNA_genes_column, miRNA_pvalue_column=miRNA_pvalue_column, miRNA_pvalue_threshold=miRNA_pvalue_threshold,
                                                        count_threshold=count_threshold, benjamini_threshold=benjamini_threshold)
    os.chdir(current_dir)

@pytest.mark.asyncio
async def test_bulk_rnaseq_mapping():
    os.chdir(current_dir)
    input_file_path = input_file_path_bulk
    sheet_name_paths = pathways_sheet_name
    sheet_name_genes = genes_sheet_name
    my_pipeline = PipelineAsync(
        input_file_path=input_file_path,
        sheet_name_paths=sheet_name_paths,
        sheet_name_genes=sheet_name_genes,
        genes_column = genes_column,
        log2fc_column=log2fc_column,
        save_to_eps=save_to_eps
    )

    await my_pipeline.bulk_rnaseq_mapping()
    os.chdir(current_dir)

@pytest.mark.asyncio
async def test_output_folder_scheme():
    input_file_path = input_file_path_g
    my_pipeline = PipelineAsync(input_file_path=input_file_path,
                           sheet_name_paths=sheet_name_paths,
                           sheet_name_genes=sheet_name_genes,
                           save_to_eps=save_to_eps,
                           genes_column=genes_column,
                           log2fc_column=log2fc_column,
                           output_folder_name=output_folder_name,
                           folder_extension='with_extension')

    await my_pipeline.gene_expression(count_threshold=count_threshold, benjamini_threshold=benjamini_threshold,pathway_pvalue_threshold=None)
    os.chdir(current_dir)


@pytest.mark.asyncio    
async def test_demirs_per_gene():
    os.chdir(current_dir)
    input_file_path = input_file_path_g
    my_pipeline = PipelineAsync(
        input_file_path=input_file_path,
        sheet_name_paths=sheet_name_paths,
        sheet_name_genes=sheet_name_genes,
        genes_column=genes_column,
        log2fc_column=log2fc_column,
        save_to_eps=save_to_eps)

    await my_pipeline.demirs_per_gene(miRNA_file_path=miRNA_path_quant, miRNA_genes_column=miRNA_genes_column,
                                                        miRNA_pvalue_column=miRNA_pvalue_column, miRNA_pvalue_threshold=miRNA_pvalue_threshold,
                                                        miRNA_ID_column=miRNA_column, 
                                                        count_threshold=count_threshold, benjamini_threshold=benjamini_threshold)
    os.chdir(current_dir)

@pytest.mark.asyncio    
async def test_dmps_per_gene():
    os.chdir(current_dir)
    input_file_path = input_file_path_g
    my_pipeline = PipelineAsync(
        input_file_path=input_file_path,
        sheet_name_paths=sheet_name_paths,
        sheet_name_genes=sheet_name_genes,
        genes_column=genes_column,
        log2fc_column=log2fc_column,
        save_to_eps=save_to_eps)

    await my_pipeline.dmps_per_gene(methylation_file_path=methylation_path_quant, methylation_genes_column=methylation_genes_column,
                                                             methylation_pvalue_column=methylation_pvalue_column, methylation_pvalue_threshold=methylation_pvalue_threshold,
                                                            methylation_probe_column=methylation_probe_column,probes_to_cgs=False,
                                                             count_threshold=count_threshold, benjamini_threshold=benjamini_threshold)
    os.chdir(current_dir)

@pytest.mark.asyncio
async def test_dmps_per_gene_correct_probes():
    os.chdir(current_dir)
    input_file_path = input_file_path_g
    my_pipeline = PipelineAsync(
        input_file_path=input_file_path,
        sheet_name_paths=sheet_name_paths,
        sheet_name_genes=sheet_name_genes,
        genes_column=genes_column,
        log2fc_column=log2fc_column,
        save_to_eps=save_to_eps)

    my_pipeline.dmps_per_gene(methylation_file_path=methylation_path_quant, methylation_genes_column=methylation_genes_column,
                                                             methylation_pvalue_column=methylation_pvalue_column, methylation_pvalue_threshold=methylation_pvalue_threshold,
                                                            methylation_probe_column=methylation_probe_column,probes_to_cgs=True,
                                                             count_threshold=count_threshold, benjamini_threshold=benjamini_threshold)
    os.chdir(current_dir)

@pytest.mark.asyncio
async def test_gene_expression_with_compounds():
    my_pipeline = PipelineAsync(input_file_path=input_file_path_g,
                           sheet_name_paths=sheet_name_paths,
                           sheet_name_genes=sheet_name_genes,
                           log2fc_column=log2fc_column,
                           genes_column=genes_column,
                           save_to_eps=True,
                           compounds_list=['C00031' , 'C00162'] , folder_extension='compounds')

    await my_pipeline.gene_expression(count_threshold=1 , benjamini_threshold=benjamini_threshold)
    os.chdir(current_dir)