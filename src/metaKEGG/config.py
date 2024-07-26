pale_yellow = '#ffffd1'
gray = '#d6d6d6'
dark_gray = '#8a8a8a'
white = '#ffffff'
lipid_yellow = '#e7ba52'

tsv_suffixes = ['.tsv']
csv_suffixes = ['.csv']
excel_suffixes = ['.xlsx', '.xls']

analysis_types = {1 : 'Gene expression' , 
                  2 : 'Transcript expression' ,
                  3 : 'Bulk RNAseq mapping', 
                  4 : 'Multiple inputs', 
                  5 : 'Methylated genes', 
                  6 : 'DMPs per gene',
                  7 : 'miRNA target genes',
                  8 : 'DEmiRs per gene',
                  9 : 'Methylated + miRNA target genes'}

analysis_types_to_execute = {1 : 'single_input_genes()' , 
                  2 : 'single_input_transcripts()' ,
                  3 : 'single_input_genes_bulk_mapping()', 
                  4 : 'multiple_inputs()', 
                  5 : 'single_input_with_methylation()', 
                  6 : 'single_input_with_methylation_quantification()',
                  7 : 'single_input_with_miRNA()',
                  8 : 'single_input_with_miRNA_quantification()',
                  9 : 'single_input_with_methylation_and_miRNA()'
                  }
