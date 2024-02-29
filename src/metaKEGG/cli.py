from .modules.pipeline import Pipeline
import argparse

def parse_args():
    """
    Parse command line arguments.

    Returns:
        argparse.Namespace: Parsed command line arguments.
    """
    parser = argparse.ArgumentParser(prog='metaKEGG',
                                     description='Pipeline for automatic enichement analysis mapping on KEGG pathways.')

    parser.add_argument('--input_file_path', required=True, help='Path to the input file (Excel format) or list of input files. Can be a David analysis output, or RNAseq')
    parser.add_argument('--sheet_name_paths', required=True, default='pathways', help='Sheet name containing the pathway information (see docs). Has to apply to all input files in case of multiple.')
    parser.add_argument('--sheet_name_genes', required=True, default='gene_metrics', help='Sheet name for gene information (see docs). Has to apply to all input files in case of multiple.')
    parser.add_argument('--genes_column', required=True, default='gene_symbol', help='Column name for gene symbols in the sheet_name_genes')
    parser.add_argument('--log2fc_column', required=True, default='logFC', help='Column name for log2fc values in the sheet_name_genes')
    parser.add_argument('--analysis_type', required=True, type=int, default=None, choices=[1, 2, 3, 4, 5, 6, 7] ,help='Analysis type (1-7)')
    parser.add_argument('--count_threshold', required=False, type=int, default=2, help='Minimum number of genes per pathway, for pathway to be drawn. Default value : 2')
    parser.add_argument('--input_label', required=False, default=None, help='Input label or list of labels for multiple inputs')
    parser.add_argument('--folder_extension', required=False, default=None, help='Folder extension to be appended to the default naming scheme. If None and default folder exists, will overwrite folder')
    parser.add_argument('--methylation_path', required=False, default=None, help='Path to methylation data (Excel , CSV or TSV format)')
    parser.add_argument('--methylation_pvalue', required=False, default=None, help='Column name for methylation p-value')
    parser.add_argument('--methylation_genes', required=False, default=None, help='Column name for methylation gene symbols')
    parser.add_argument('--methylation_pvalue_thresh', required=False, default=0.05, help='P-value threshold for the methylation values')
    parser.add_argument('--miRNA_path', required=False, default=None, help='Path to miRNA data (Excel , CSV or TSV format)')
    parser.add_argument('--miRNA_pvalue', required=False, default=None, help='Column name for miRNA p-value')
    parser.add_argument('--miRNA_genes', required=False, default=None, help='Column name for miRNA gene symbols')
    parser.add_argument('--miRNA_pvalue_thresh', required=False, default=0.05, help='P-value threshold for the miRNA values')
    parser.add_argument('--benjamini_threshold', required=False, default=None, help='Benjamini Hochberg p-value threshold for the pathway')
    parser.add_argument('--save_to_eps', required=False, default=False, help='True/False statement to save the maps and colorscales or legends as seperate .eps files in addition to the .pdf exports')

    return parser.parse_args()


def main():
    """
    Main function to execute the pipeline based on the parsed command line arguments.
    """
    args = parse_args()

    if args.analysis_type is not None:
        Pipeline(
            input_file_path=args.input_file_path,
            sheet_name_paths=args.sheet_name_paths,
            sheet_name_genes=args.sheet_name_genes,
            genes_column=args.genes_column,
            log2fc_column=args.log2fc_column,
            analysis_type=args.analysis_type,
            count_threshold=args.count_threshold,
            input_label=args.input_label,
            folder_extension=args.folder_extension,
            methylation_path=args.methylation_path,
            methylation_pvalue=args.methylation_pvalue,
            methylation_genes=args.methylation_genes,
            methylation_pvalue_thresh=args.methylation_pvalue_thresh,
            miRNA_path=args.miRNA_path,
            miRNA_pvalue=args.miRNA_pvalue,
            miRNA_genes=args.miRNA_genes,
            miRNA_pvalue_thresh=args.miRNA_pvalue_thresh,

            benjamini_threshold=args.benjamini_threshold,
            save_to_eps=args.save_to_eps,
        )
    else:
        print(f'Could not initiate the pipeline dues to missing value for the analysis_type. You provided {args.analysis_type}')

if __name__ == "__main__":
    main()