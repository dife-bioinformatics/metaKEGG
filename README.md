[![Stars](https://img.shields.io/github/stars/dife-bioinformatics/metaKEGG?style=flat&logo=GitHub&color=yellow)](https://github.com/dife-bioinformatics/metaKEGG)
[![PyPI](https://img.shields.io/pypi/v/metaKEGG?logo=PyPI)](https://pypi.org/project/metaKEG)


# `metaKEGG`

metaKEGG is a fully integrated solution with class-leading features to visualize the KEGG pathway enrichment analysis results from the DAVID Functional Annotation Tool, or RNAseq data.

## Table of Contents
- [Disclaimer](#disclaimer)
- [Installing metaKEGG](#installing-metakegg)
    - [Environment preparation](#environment-preparation)
    - [Install from PyPI](#install-from-pypi)
    - [Local installation with venv and requirement.txt in Windows](#local-installation-with-venv-and-requirementtxt-in-windows)
    - [Local installation with conda env and environment.yml](#local-installation-with-conda-env-and-environmentyml)
- [Getting started](#getting-started)
    - [CLI usage](#cli-usage)
    - [Programmatic/Library usage](#programmaticlibrary-usage)
- [Example usage](#example-usage)
    - [Gene expression](#gene-expression)
    - [Transcript expression](#transcript-expression)
    - [Bulk RNAseq mapping](#bulk-rnaseq-mapping)
    - [Multiple inputs](#multiple-inputs)
    - [Methylated genes](#methylated-genes)
    - [DMPs per gene](#dmps-per-genes)
    - [miRNA target genes](#mirna-target-genes)
    - [DEmiRs per gene](#demirs-per-gene)    
    - [Methylated + miRNA target genes](#methylated--mirna-target-genes)    
    - [Example of using multiple modules with one initialization](#example-of-using-multiple-modules-with-one-initialization)

## Disclaimer

metaKEGG uses the KEGG API, and is restricted to academic use by academic users belonging to academic institutions.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.
You should have received a copy of the GNU Affero General Public License along with this
program. If not, see [GNU Affero General Public License](https://www.gnu.org/licenses/#AGPL).

Author: Michail Lazaratos, Deutsches Institut für Ernährungsforschung Potsdam-Rehbrücke / German Institute of Human Nutrition Potsdam-Rehbrücke (DIfE)

## Installing metaKEGG

Either clone from this GitHub repository or install from 
[PyPI](https://pypi.org/) (recommended).


### Environment preparation

1. Install an [Anaconda Distribution](https://docs.anaconda.com/free/anaconda/install/) or [Miniconda](https://docs.anaconda.com/free/miniconda/miniconda-other-installer-links/)


2. Create a new conda environment to install metaKEGG in.

Example using python 3.9. The pakcage was tested with this version. Later version should work but there is no guarantee.

```
conda create -n metaKEGG_env python=3.9
conda activate metaKEGG_env
```

### Install from PyPI

Install the package directy from PyPI.

```
pip install metaKEGG
```

### Local installation with venv and requirement.txt in Windows

To create a copy of the dev environment using venv, after cloning this repository do the following:

```
cd /path/to/parent/folder
python -m venv venv
venv\Scripts\activate
python -m pip install -r requirements.txt
pip install -e .
```

### Local installation with conda env and environment.yml

To create a copy of the dev environment using conda, after cloning this repository do the following:

```
cd /path/to/parent/folder
conda env create -f .\conda_env\environment.yml
conda activate metaKEGG_env
pip install -e .
```

## Getting started

After successfully installing metaKEGG, you can use it two ways.

### CLI usage

1. Simply open a terminal and make sure you activate the environment metaKEGG is installed in.

2. Type `metaKEGG` in the console. This will prompt all the input arguments. Type `metaKEGG -h` for a detailed description of each argument.

CLI usage can be also wrapped in bash scripts and integrated in pipelines.

### Programmatic/Library usage

Import the class

Note: The Pipeline class requires specific arguments for initialization.
```
import metaKEGG
metaKEGG.Pipeline(**kwargs)
```

Alternatively
```
from metaKEGG import Pipeline
Pipeline(**kwargs)
```

## Example usage

In the `/examples` directory you can find example files to perform all analysis types provided by metaKEGG.
The paths to these files will be used for demonstration purposes below.

The analysis types available are:

```
1 : Gene expression 
2 : Transcript expression
3 : Bulk RNAseq mapping 
4 : Multiple inputs 
5 : Methylated genes 
6 : DMPs per gene
7 : miRNA target genes
8 : DEmiRs per gene
9 : Methylated + miRNA target genes
```

Default values are:

```
sheet_name_paths = "pathways"
sheet_name_genes = "gene_metrics"
genes_column = "gene_symbol"
log2fc_column = "logFC"
input_label = None
count_threshold = 2
pathway_pvalue = None
benjamini_threshold = None
save_to_eps = False
folder_extension = None
output_folder_name = None

methylation_path = None (required for 5, 6 & 9)
methylation_genes = None (required for 5, 6 & 9)
methylation_probe_column = None (required for 6)
methylation_pvalue = None 
methylation_pvalue_thresh=0.05
probes_to_cgs=False

miRNA_path = None (required for 7, 8 & 9)
miRNA_genes = None (required for  7, 8 & 9)
miRNA_ID_column = None (required for 8)
miRNA_pvalue = None
miRNA_pvalue_thresh=0.05

analysis_type = None (required to be set between 1 and 9 for CLI usage)                  
```

If usage is programmatic and analysis_type = None see section [Example of using multiple modules with one initialization](#one-initialization) below for executing the pipeline.

-----
### Gene expression

This function takes a single input file as an argument and maps the detected genes on the enriched KEGG reference pathway, and colors them according to their log2FC values.

1. Define input arguments

```
input_file_path = "examples/single_input_genes.xlsx"
input_label = "input1"
sheet_name_paths = "pathways"
sheet_name_genes = "gene_metrics"
genes_column = "gene_symbol"
log2fc_column = "logFC"
analysis_type = 1
count_threshold = 2
pathway_pvalue = None
benjamini_threshold = None
save_to_eps = False
folder_extension = None
```

2. Run analysis

```
import metaKEGG

metaKEGG.Pipeline(input_file_path=input_file_path, input_label=input_label, sheet_name_paths=sheet_name_paths, sheet_name_genes=sheet_name_genes, genes_column=genes_column, log2fc_column=log2fc_column, analysis_type=analysis_type, count_threshold=count_threshold, pathway_pvalue=pathway_pvalue ,benjamini_threshold=benjamini_threshold, save_to_eps=save_to_eps, folder_extension=folder_extension)
```

Alternatively using the CLI

```
metaKEGG --input_file_path="examples/single_input_genes.xlsx"  --input_label="input1"  --sheet_name_paths="pathways" --sheet_name_genes="gene_metrics" --genes_column="gene_symbol" --log2fc_column="logFC" --analysis_type=1 --count_threshold=2 --pathway_pvalue=None --benjamini_threshold=None --save_to_eps=False --folder_extension=None
```

-----
### Transcript expression

This function takes a single input file as an argument and maps the detected transcripts on the enriched KEGG reference pathway, and colors them according to their log2FC values.

**_NOTE:_**  Pathway enrichement analysis with the DAVID Functional Annotation Tool, should be performed using transcript IDs.

1. Define input arguments

```
input_file_path = "examples/single_input_transcripts.xlsx"
input_label = "input1"
sheet_name_paths = "pathways"
sheet_name_genes = "transcript_metrics"
genes_column = "gene_symbol"
log2fc_column = "logFC"
analysis_type = 2
count_threshold = 2
pathway_pvalue = None
benjamini_threshold = None
save_to_eps = False
folder_extension = None
```

2. Run analysis

```
import metaKEGG

metaKEGG.Pipeline(input_file_path=input_file_path, input_label=input_label, sheet_name_paths=sheet_name_paths, sheet_name_genes=sheet_name_genes, genes_column=genes_column, log2fc_column=log2fc_column, analysis_type=analysis_type, count_threshold=count_threshold, pathway_pvalue=pathway_pvalue, benjamini_threshold=benjamini_threshold, save_to_eps=save_to_eps, folder_extension=folder_extension)
```

Alternatively using the CLI

```
metaKEGG --input_file_path="examples/single_input_transcripts.xlsx" --input_label="input1" --sheet_name_paths="pathways" --sheet_name_genes="transcript_metrics" --genes_column="gene_symbol" --log2fc_column="logFC" --analysis_type=2 --count_threshold=2 --pathway_pvalue=None --benjamini_threshold=None --save_to_eps=False --folder_extension=None
```
-----
### Bulk RNAseq mapping

This function takes RANseq data, as single input file argument, maps the genes on a provided list of target pathways (assuming they are also found in the target pathways), and colors them according to their log2FC values.

1. Define input arguments

```
input_file_path = "examples/single_input_bulk.xlsx"
input_label = "input1"
genes_column = "gene_symbol"
log2fc_column = "logFC"
sheet_name_paths = "pathways"
sheet_name_genes = "gene_metrics"
analysis_type = 3
save_to_eps = False
folder_extension = None
```

2. Run analysis

```
import metaKEGG

metaKEGG.Pipeline(input_file_path=input_file_path, input_label=input_label, sheet_name_paths=sheet_name_paths, sheet_name_genes=sheet_name_genes, genes_column=genes_column, log2fc_column=log2fc_column, analysis_type=analysis_type, save_to_eps=save_to_eps, folder_extension=folder_extension)
```

Alternatively using the CLI

```
metaKEGG --input_file_path="examples/single_input_genes.xlsx" --input_label="input1" --sheet_name_paths="pathways" --sheet_name_genes="gene_metrics" --genes_column="gene_symbol" --log2fc_column="logFC" --analysis_type=3 --save_to_eps=False --folder_extension=None
```
-----
### Multiple inputs

This function takes a list of inputs file as an argument and only maps pathways that are found in all of the inputs files.
For a common pathway, it will generate all possible states for a gene, from each individual input, to all possible combinations and assignes a unique color code to each combination.
The detected genes are mapped enriched KEGG reference pathway, based on the state they're in.


1. Define input arguments

```
input_file_path = ["examples/single_input_genes.xlsx",
                   "examples/multiple_inputs_1.xlsx",
                   "examples/multiple_inputs_2.xlsx"]
                   
input_label = ["input1" , "input2" , "input3"]
sheet_name_paths = "pathways"
sheet_name_genes = "gene_metrics"
genes_column = "gene_symbol"
log2fc_column = "logFC"
analysis_type = 4
count_threshold = 2
pathway_pvalue = None
benjamini_threshold = None
save_to_eps = False
folder_extension = None
```

2. Run analysis

```
import metaKEGG

metaKEGG.Pipeline(input_file_path=input_file_path, input_label=input_label, sheet_name_paths=sheet_name_paths, sheet_name_genes=sheet_name_genes, genes_column=genes_column, log2fc_column=log2fc_column, analysis_type=analysis_type, count_threshold=count_threshold, pathway_pvalue=pathway_pvalue, benjamini_threshold=benjamini_threshold, save_to_eps=save_to_eps, folder_extension=folder_extension)
```

Alternatively using the CLI

```
metaKEGG --input_file_path=["examples/single_input_genes.xlsx",
                   "examples/multiple_inputs_1.xlsx",
                   "examples/multiple_inputs_2.xlsx"] 
         --input_label=["input1" , "input2" , "input3"] 
         --sheet_name_paths="pathways" --sheet_name_genes="gene_metrics" --genes_column="gene_symbol" --log2fc_column="logFC" --analysis_type=4 --count_threshold=2 --pathway_pvalue=None --benjamini_threshold=None --save_to_eps=False --folder_extension=None
```
-----
### Methylated genes

This function takes a single input file and a methylation metadata file as arguments and maps the detected genes on the enriched KEGG reference pathway, and colors them according to their methylation state. The state is defined as a binary reprsentation, depending if DMPs corresponding to a given gene are detected, or not.


1. Define input arguments

```
input_file_path = "examples/single_input_genes.xlsx"
input_label = "input1"
sheet_name_paths = "pathways"
sheet_name_genes = "gene_metrics"
genes_column = "gene_symbol"
log2fc_column = "logFC"
analysis_type = 5
methylation_path = "examples/methylation.csv"
methylation_genes = "methylation_gene_symbol"
methylation_pvalue = "methylation_pval"
methylation_pvalue_thresh = 0.05
count_threshold = 2
pathway_pvalue = None
benjamini_threshold = None
save_to_eps = False
folder_extension = None
```

2. Run analysis

```
import metaKEGG

metaKEGG.Pipeline(input_file_path=input_file_path, input_label=input_label, sheet_name_paths=sheet_name_paths, sheet_name_genes=sheet_name_genes, genes_column=genes_column, log2fc_column=log2fc_column, analysis_type=analysis_type, methylation_path=methylation_path, methylation_genes=methylation_genes, methylation_pvalue=methylation_pvalue, methylation_pvalue_thresh=methylation_pvalue_thresh, count_threshold=count_threshold, pathway_pvalue=pathway_pvalue, benjamini_threshold=benjamini_threshold, save_to_eps=save_to_eps, folder_extension=folder_extension)
```

Alternatively using the CLI

```
metaKEGG --input_file_path="examples/single_input_genes.xlsx" --input_label="input1" --sheet_name_paths="pathways" --sheet_name_genes="gene_metrics" 
         --genes_column="gene_symbol" --log2fc_column="logFC" --analysis_type=5 
         --methylation_path="examples/methylation.csv" --methylation_genes="methylation_gene_symbol" --methylation_pvalue="methylation_pval" --methylation_pvalue_thresh=0.05 
         --count_threshold=2 --pathway_pvalue=None --benjamini_threshold=None --save_to_eps=False --folder_extension=None
```
-----
### DMPs per gene

This function takes a single input file and a methylation metadata file as arguments and maps the detected genes on the enriched KEGG reference pathway. It generates bins to quantify the number of DMPs that correspond to a given gene, and colors a gege according to its DMP bin. The function also returns the quantification histogram plots, both in a grouped and an absolute count representation.

1. Define input arguments

```
input_file_path = "examples/single_input_genes.xlsx"
input_label = "input1"
sheet_name_paths = "pathways"
sheet_name_genes = "gene_metrics"
genes_column = "gene_symbol"
log2fc_column = "logFC"
analysis_type = 6
methylation_path = "examples/methylation_for_quantification.csv"
methylation_genes = "methylation_gene_symbol"
methylation_pvalue = "methylation_pval"
methylation_pvalue_thresh = 0.05
methylation_probe_column = "CG_ID"
probes_to_cgs=False
count_threshold = 2
pathway_pvalue = None
benjamini_threshold = None
save_to_eps = False
folder_extension = None
```

2. Run analysis

```
import metaKEGG

metaKEGG.Pipeline(input_file_path=input_file_path, input_label=input_label, sheet_name_paths=sheet_name_paths, sheet_name_genes=sheet_name_genes, genes_column=genes_column, log2fc_column=log2fc_column, analysis_type=analysis_type, methylation_path=methylation_path, methylation_genes=methylation_genes, methylation_pvalue=methylation_pvalue, methylation_pvalue_thresh=methylation_pvalue_thresh, methylation_probe_column=methylation_probe_column, probes_to_cgs=probes_to_cgs, count_threshold=count_threshold, pathway_pvalue=pathway_pvalue, benjamini_threshold=benjamini_threshold, save_to_eps=save_to_eps, folder_extension=folder_extension)
```

Alternatively using the CLI

```
metaKEGG --input_file_path="examples/single_input_genes.xlsx" --input_label="input1" --sheet_name_paths="pathways" --sheet_name_genes="gene_metrics" 
         --genes_column="gene_symbol" --log2fc_column="logFC" --analysis_type=6 
         --methylation_path="examples/methylation.csv" --methylation_genes="methylation_gene_symbol" --methylation_pvalue="methylation_pval" --methylation_pvalue_thresh=0.05 
         --methylation_probe_column="CG_ID" --probes_to_cgs=False 
         --count_threshold=2 --pathway_pvalue=None --benjamini_threshold=None --save_to_eps=False --folder_extension=None
```

**_NOTE:_**  When using probes_to_cgs=True, the pipeline will split the CG probes by the underscore '_' character and keep the first part, essentially correcting for different probe chemistry that could occur in the same position. Example format is cg00000000_BC21 and cg00000000_TC21, which would be counted as two separate probes targeting the same gene. Using the argument probes_to_cgs with True, the probes become cg00000000 & cg00000000, and duplicated entries per gene are eliminated, essentially counting one probe for the target gene.

------
### miRNA target genes

This function takes a single input file and a miRNA metadata file as arguments and maps the detected genes on the enriched KEGG reference pathway, and colors them according to their miRNA state. The state is defined as a binary reprsentation, depending if DEmiRs are targeting a given gene, or not.

1. Define input arguments

```
input_file_path = "examples/single_input_genes.xlsx"
input_label = "input1"
sheet_name_paths = "pathways"
sheet_name_genes = "gene_metrics"
genes_column = "gene_symbol"
log2fc_column = "logFC"
analysis_type = 7
miRNA_path = "examples/miRNA.tsv"
miRNA_genes = "miRNA_gene_symbol"
miRNA_pvalue = "miRNA_pval"
miRNA_pvalue_thresh=0.05
pathway_pvalue = None
count_threshold = 2
benjamini_threshold = None
save_to_eps = False
folder_extension = None
```

2. Run analysis

```
import metaKEGG

metaKEGG.Pipeline(input_file_path=input_file_path, input_label=input_label, sheet_name_paths=sheet_name_paths, sheet_name_genes=sheet_name_genes, genes_column=genes_column, log2fc_column=log2fc_column, analysis_type=analysis_type, miRNA_path=miRNA_path, miRNA_genes=miRNA_genes, miRNA_pvalue=miRNA_pvalue, miRNA_pvalue_thresh=miRNA_pvalue_thresh,count_threshold=count_threshold, pathway_pvalue=pathway_pvalue, benjamini_threshold=benjamini_threshold, save_to_eps=save_to_eps, folder_extension=folder_extension)
```

Alternatively using the CLI

```
metaKEGG --input_file_path="examples/single_input_genes.xlsx" --input_label="input1" --sheet_name_paths="pathways" --sheet_name_genes="gene_metrics" 
         --genes_column="gene_symbol" --log2fc_column="logFC" --analysis_type=7 
         --miRNA_path="examples/miRNA.tsv" --miRNA_genes="miRNA_gene_symbol" --miRNA_pvalue="miRNA_pval" --miRNA_pvalue_thresh=0.05 
         --count_threshold=2 --pathway_pvalue=None --benjamini_threshold=None --save_to_eps=False --folder_extension=None
```
------
### DEmiRs per gene

This function takes a single input file and a miRNA metadata file as arguments and maps the detected genes on the enriched KEGG reference pathway. It generates bins to quantify the number of DEmiRs that correspond to a given gene, and colors a gege according to its DEmiR bin. The function also returns the quantification histogram plots, both in a grouped and an absolute count representation.

1. Define input arguments

```
input_file_path = "examples/single_input_genes.xlsx"
input_label = "input1"
sheet_name_paths = "pathways"
sheet_name_genes = "gene_metrics"
genes_column = "gene_symbol"
log2fc_column = "logFC"
analysis_type = 8
miRNA_path = "examples/miRNA_for_quantification.tsv"
miRNA_genes = "miRNA_gene_symbol"
miRNA_pvalue = "miRNA_pval"
miRNA_pvalue_thresh=0.05
miRNA_column = "miRNA_ID"
pathway_pvalue = None
count_threshold = 2
benjamini_threshold = None
save_to_eps = False
folder_extension = None
```

2. Run analysis

```
import metaKEGG

metaKEGG.Pipeline(input_file_path=input_file_path, input_label=input_label, sheet_name_paths=sheet_name_paths, sheet_name_genes=sheet_name_genes, genes_column=genes_column, log2fc_column=log2fc_column, analysis_type=analysis_type, miRNA_path=miRNA_path, miRNA_genes=miRNA_genes, miRNA_pvalue=miRNA_pvalue, miRNA_pvalue_thresh=miRNA_pvalue_thresh,miRNA_column=miRNA_column, count_threshold=count_threshold, pathway_pvalue=pathway_pvalue, benjamini_threshold=benjamini_threshold, save_to_eps=save_to_eps, folder_extension=folder_extension)
```

Alternatively using the CLI

```
metaKEGG --input_file_path="examples/single_input_genes.xlsx" --input_label="input1" --sheet_name_paths="pathways" --sheet_name_genes="gene_metrics" 
         --genes_column="gene_symbol" --log2fc_column="logFC" --analysis_type=8 
         --miRNA_path="examples/miRNA.tsv" --miRNA_genes="miRNA_gene_symbol" --miRNA_pvalue="miRNA_pval" --miRNA_pvalue_thresh=0.05 --miRNA_column="miRNA_ID" 
         --count_threshold=2 --pathway_pvalue=None --benjamini_threshold=None --save_to_eps=False --folder_extension=None
```

-------
### Methylated + miRNA target genes

This function takes a single input file, a methylation, and a miRNA metadata file as arguments and maps the detected genes on the enriched KEGG reference pathway. Genes are colored according to their methylation and miRNA states. The states is defined as a binary reprsentations of the methylation and miRNA combinations.

1. Define input arguments

```
input_file_path = "examples/single_input_genes.xlsx"
input_label = "input1"
sheet_name_paths = "pathways"
sheet_name_genes = "gene_metrics"
genes_column = "gene_symbol"
log2fc_column = "logFC"
analysis_type = 9
methylation_path = "examples/methylation.csv"
methylation_genes = "methylation_gene_symbol"
methylation_pvalue = "methylation_pval"
methylation_pvalue_thresh = 0.05
miRNA_path = "examples/miRNA.tsv"
miRNA_genes = "miRNA_gene_symbol"
miRNA_pvalue = "miRNA_pval"
miRNA_pvalue_thresh = 0.05
count_threshold = 2
pathway_pvalue = None
benjamini_threshold = None
save_to_eps = False
folder_extension = None
```

2. Run analysis

```
import metaKEGG

metaKEGG.Pipeline(input_file_path=input_file_path, input_label=input_label, sheet_name_paths=sheet_name_paths, sheet_name_genes=sheet_name_genes, genes_column=genes_column, log2fc_column=log2fc_column, analysis_type=analysis_type, methylation_path=methylation_path, methylation_genes=methylation_genes, methylation_pvalue=methylation_pvalue, methylation_pvalue_thresh=methylation_pvalue_thresh, miRNA_path=miRNA_path, miRNA_genes=miRNA_gene, miRNA_pvalue=miRNA_pvalue, miRNA_pvalue_thresh=miRNA_pvalue_thresh, count_threshold=count_threshold, pathway_pvalue=pathway_pvalue, benjamini_threshold=benjamini_threshold, save_to_eps=save_to_eps, folder_extension=folder_extension)
```

Alternatively using the CLI

```
metaKEGG --input_file_path="examples/single_input_genes.xlsx" --input_label="input1" --sheet_name_paths="pathways" --sheet_name_genes="gene_metrics" 
         --genes_column="gene_symbol" --log2fc_column="logFC" --analysis_type=9 
         --methylation_path="examples/methylation.csv" --methylation_genes="methylation_gene_symbol" --methylation_pvalue="methylation_pval" --methylation_pvalue_thresh=0.05 
         --miRNA_path="examples/miRNA.tsv" --miRNA_genes="miRNA_gene_symbol" --miRNA_pvalue="miRNA_pval" --miRNA_pvalue_thresh=0.05 
         --count_threshold=2 --pathway_pvalue=None --benjamini_threshold=None --save_to_eps=False --folder_extension=None
```
-----
### Example of using multiple modules with one initialization

Define as many input arguments that can be used for multiple modules and leave the `analysis_type` to `None`

1. Define input arguments

```
input_file_path = "examples/single_input_genes.xlsx"
input_label = "input1"
sheet_name_paths = "pathways"
sheet_name_genes = "gene_metrics"
genes_column = "gene_symbol"
log2fc_column = "logFC"
methylation_path = "examples/methylation.csv"
methylation_genes = "methylation_gene_symbol"
methylation_pvalue = "methylation_pval"
methylation_pvalue_thresh = 0.05
miRNA_path = "examples/miRNA.tsv"
miRNA_genes = "miRNA_gene_symbol"
miRNA_pvalue = "miRNA_pval"
miRNA_pvalue_thresh = 0.05
count_threshold = 2
benjamini_threshold = None
save_to_eps = False
folder_extension = None

analysis_type = None
```

2. Run analysis

```
import metaKEGG

m = metaKEGG.Pipeline(input_file_path=input_file_path, input_label=input_label, sheet_name_paths=sheet_name_paths, sheet_name_genes=sheet_name_genes, genes_column=genes_column, log2fc_column=log2fc_column, analysis_type=analysis_type, methylation_path=methylation_path, methylation_genes=methylation_genes, methylation_pvalue=methylation_pvalue, methylation_pvalue_thresh=methylation_pvalue_thresh, miRNA_path=miRNA_path, miRNA_gene=miRNA_gene, miRNA_pvalue=miRNA_pvalue, miRNA_pvalue_thresh=miRNA_pvalue_thresh, count_threshold=count_threshold, benjamini_threshold=benjamini_threshold, save_to_eps=save_to_eps, folder_extension=folder_extension)
```

```
m.single_input_genes()
m.single_input_with_methylation()
m.single_input_with_miRNA()
m.single_input_with_methylation_and_miRNA()
```

-----