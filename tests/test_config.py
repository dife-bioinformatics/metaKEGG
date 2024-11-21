input_file_path_g = '../examples/single_input_genes.xlsx'
input_file_path_t = '../examples/single_input_transcripts.xlsx'

input_file_path_m = ['../examples/single_input_genes.xlsx',
                   '../examples/multiple_inputs_1.xlsx',
                   '../examples/multiple_inputs_2.xlsx']


input_label_m = ['input1' , 'input2' , 'input3']
sheet_name_paths = "pathways"
sheet_name_genes = "gene_metrics"
sheet_name_transcripts = "transcript_metrics"

methylation_file_path = '../examples/methylation.csv'
methylation_genes_column = 'methylation_gene_symbol'
methylation_pvalue_column = 'methylation_pval'
methylation_pvalue_column = None
methylation_probe_column = 'CG_ID'

miRNA_file_path = '../examples/miRNA.tsv'
miRNA_genes_column = 'miRNA_gene_symbol'
miRNA_pvalue_column = 'miRNA_pval'
miRNA_pvalue_column = None
miRNA_column = 'miRNA_ID'

input_file_path_bulk = '../examples/single_input_bulk.xlsx'
pathways_sheet_name = 'pathways'
genes_sheet_name = 'gene_metrics'
genes_column = 'gene_symbol'
log2fc_column = 'logFC'

analysis_type = None
save_to_eps = True
count_threshold=2
benjamini_threshold=None
methylation_pvalue_threshold=None
miRNA_pvalue_threshold=None
output_folder_name = 'my_folder_name'

miRNA_path_quant = '../examples/miRNA_for_quantification.tsv'
methylation_path_quant = '../examples/methylation_for_quantification.csv'