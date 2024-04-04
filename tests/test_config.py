input_file_path_g = '../examples/single_input_genes.xlsx'
input_file_path_t = '../examples/single_input_transcripts.xlsx'

input_file_path_m = ['../examples/single_input_genes.xlsx',
                   '../examples/multiple_inputs_1.xlsx',
                   '../examples/multiple_inputs_2.xlsx']


input_label_m = ['input1' , 'input2' , 'input3']
input_label_g = 'input'
sheet_name_paths = "pathways"
sheet_name_genes = "gene_metrics"
sheet_name_transcripts = "transcript_metrics"

methylation_path = '../examples/methylation.csv'
methylation_gene = 'methylation_gene_symbol'
methylation_pvalue = 'methylation_pval'
methylation_pvalue = None

miRNA_path = '../examples/miRNA.tsv'
miRNA_gene = 'miRNA_gene_symbol'
miRNA_pvalue = 'miRNA_pval'
miRNA_pvalue = None


input_file_path_bulk = '../examples/single_input_bulk.xlsx'
pathways_sheet_name = 'pathways'
genes_sheet_name = 'gene_metrics'
genes_column = 'gene_symbol'
log2fc_column = 'logFC'

analysis_type = None
save_to_eps = False
count_threshold=2
benjamini_threshold=None
methylation_pvalue_thresh=None
miRNA_pvalue_thresh=None
output_folder_name = 'my_folder_name'

metadata_path_quant = '../examples/metadata_for_quantification.tsv'
