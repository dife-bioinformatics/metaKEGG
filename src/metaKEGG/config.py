pale_yellow = '#ffffd1'
gray = '#d6d6d6'
dark_gray = '#8a8a8a'
white = '#ffffff'
lipid_yellow = '#e7ba52'

tsv_suffixes = ['.tsv']
csv_suffixes = ['.csv']
excel_suffixes = ['.xlsx', '.xls']

analysis_types = {1 : 'Single input (Gene IDs)' , 
                  2 : 'Single input (Transcript IDs)' ,
                  3 : 'Multiple input (Gene IDs)', 
                  4 : 'Single input w Methylation (Gene IDs)', 
                  5 : 'Single input w miRNA (Gene IDs)', 
                  6 : 'Single input w Methylation & miRNA (Gene IDs)',
                  7 : 'Single input (Gene IDs) Bulk mapping ',
                  8 : 'Single input w Methylation & Quantification (Gene IDs)',
                  9 : 'Single input w miRNA & Quantification (Gene IDs)'}

analysis_types_to_execute = {1 : 'single_input_genes()' , 
                  2 : 'single_input_transcripts()' ,
                  3 : 'multiple_inputs()', 
                  4 : 'single_input_with_methylation()', 
                  5 : 'single_input_with_miRNA()', 
                  6 : 'single_input_with_methylation_and_miRNA()',
                  7 : 'single_input_genes_bulk_mapping()',
                  8 : 'single_input_with_methylation_quantification()',
                  9 : 'single_input_with_miRNA_quantification()'
                  }
