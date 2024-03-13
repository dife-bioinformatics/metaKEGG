# `metaKEGG`

- [Disclaimer](#disclaimer)
- [Installing metaKEGG](#installing-metaKEGG)
    + [Environment preparation](#environment-preparation)
    + [Install from PyPI](#install-from-PyPI)
    + [Local installation with venv and requirement.txt in Windows](#local-venv)
    + [Local installation with conda env and environment.yml](#local-conda)
- [Getting started](#getting-started)
    + [CLI usage](#cli-usage)
    + [Programmatic/Library usage](#programmatic-library-usage)
- [Example usage](#example-usage)
    + [Single Input Analysis (Gene IDs)](#single-input-analysis-genes)
    + [Single Input Analysis (Transcript IDs)](#single-input-analysis-transcripts)
    + [Multiple Input Analysis (Gene IDs)](#multiple-input-analysis)
    + [Single Input Analysis with Methylation (Gene IDs)](#single-input-analysis-methylation)
    + [Single Input Analysis with miRNA (Gene IDs)](#single-input-analysis-mirna)
    + [Single Input Analysis with Methylation and miRNA (Gene IDs)](#single-input-analysis-methylation-mirna)
    + [Single input (Gene IDs) Bulk mapping](#single-input-analysis-bulk)
    + [Example of using multiple modules with one initialization](#one-initialization)

## Disclaimer

metaKEGG uses the KEGG API, and is restricted to academic use by academic users belonging to academic institutions.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this
program. If not, see https://www.gnu.org/licenses/.

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

Import

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
1 : Single Input Analysis (Gene IDs)
2 : Single Input Analysis (Transcript IDs)
3 : Multiple Input Analysis (Gene IDs)
4 : Single Input Analysis with Methylation (Gene IDs)
5 : Single Input Analysis with miRNA (Gene IDs) 
6 : Single Input Analysis with Methylation and miRNA (Gene IDs)
7 : Single input (Gene IDs) Bulk mapping
```

Default values are:

```
sheet_name_paths = "pathways"
sheet_name_genes = "gene_metrics"
genes_column = "gene_symbol"
log2fc_column = "logFC"
input_label = None
count_threshold = 2
benjamini_threshold = None
save_to_eps = False
folder_extension = None

methylation_path = None (required for 4 and 6)
methylation_genes = None (required for 4 and 6)
methylation_pvalue = None 
methylation_pvalue_thresh=0.05

miRNA_path = None (required for 5 and 6)
miRNA_genes = None (required for 5 and 6)
miRNA_pvalue = None
miRNA_pvalue_thresh=0.05

analysis_type = None (required to be set between 1 and 7 for CLI usage,
                      if usage is programmatic and analysis_type = None see below for usage)
```

### Single Input Analysis (Gene IDs)

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
benjamini_threshold = None
save_to_eps = False
folder_extension = None
```

2. Run analysis

```
import metaKEGG

metaKEGG.Pipeline(input_file_path=input_file_path, input_label=input_label, sheet_name_paths=sheet_name_paths, sheet_name_genes=sheet_name_genes, genes_column=genes_column, log2fc_column=log2fc_column, analysis_type=analysis_type, count_threshold=count_threshold, benjamini_threshold=benjamini_threshold, save_to_eps=save_to_eps, folder_extension=folder_extension)
```

Alternatively using the CLI

```
metaKEGG --input_file_path="examples/single_input_genes.xlsx"  --input_label="input1"  --sheet_name_paths="pathways" --sheet_name_genes="gene_metrics" --genes_column="gene_symbol" --log2fc_column="logFC" --analysis_type=1 --count_threshold=2 --benjamini_threshold=None --save_to_eps=False --folder_extension=None
```

### Single Input Analysis (Transcript IDs)

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
benjamini_threshold = None
save_to_eps = False
folder_extension = None
```

2. Run analysis

```
import metaKEGG

metaKEGG.Pipeline(input_file_path=input_file_path, input_label=input_label, sheet_name_paths=sheet_name_paths, sheet_name_genes=sheet_name_genes, genes_column=genes_column, log2fc_column=log2fc_column, analysis_type=analysis_type, count_threshold=count_threshold, benjamini_threshold=benjamini_threshold, save_to_eps=save_to_eps, folder_extension=folder_extension)
```

Alternatively using the CLI

```
metaKEGG --input_file_path="examples/single_input_transcripts.xlsx"  --input_label="input1"  --sheet_name_paths="pathways" --sheet_name_genes="transcript_metrics" --genes_column="gene_symbol" --log2fc_column="logFC" --analysis_type=2 --count_threshold=2 --benjamini_threshold=None --save_to_eps=False --folder_extension=None
```

### Multiple Input Analysis (Gene IDs)

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
analysis_type = 3
count_threshold = 2
benjamini_threshold = None
save_to_eps = False
folder_extension = None
```

2. Run analysis

```
import metaKEGG

metaKEGG.Pipeline(input_file_path=input_file_path, input_label=input_label, sheet_name_paths=sheet_name_paths, sheet_name_genes=sheet_name_genes, genes_column=genes_column, log2fc_column=log2fc_column, analysis_type=analysis_type, count_threshold=count_threshold, benjamini_threshold=benjamini_threshold, save_to_eps=save_to_eps, folder_extension=folder_extension)
```

Alternatively using the CLI

```
metaKEGG --input_file_path=["examples/single_input_genes.xlsx",
                   "examples/multiple_inputs_1.xlsx",
                   "examples/multiple_inputs_2.xlsx"] 
         --input_label=["input1" , "input2" , "input3"] 
         --sheet_name_paths="pathways" --sheet_name_genes="gene_metrics" --genes_column="gene_symbol" --log2fc_column="logFC" --analysis_type=3 --count_threshold=2 --benjamini_threshold=None --save_to_eps=False --folder_extension=None
```

### Single Input Analysis with Methylation (Gene IDs)

1. Define input arguments

```
input_file_path = "examples/single_input_genes.xlsx"
input_label = "input1"
sheet_name_paths = "pathways"
sheet_name_genes = "gene_metrics"
genes_column = "gene_symbol"
log2fc_column = "logFC"
analysis_type = 4
methylation_path = "examples/methylation.csv"
methylation_gene = "methylation_gene_symbol"
methylation_pvalue = "methylation_pval"
methylation_pvalue_thresh = 0.05,
count_threshold = 2
benjamini_threshold = None
save_to_eps = False
folder_extension = None
```

2. Run analysis

```
import metaKEGG

metaKEGG.Pipeline(input_file_path=input_file_path, input_label=input_label, sheet_name_paths=sheet_name_paths, sheet_name_genes=sheet_name_genes, genes_column=genes_column, log2fc_column=log2fc_column, analysis_type=analysis_type, methylation_path=methylation_path, methylation_gene=methylation_gene, methylation_pvalue=methylation_pvalue, methylation_pvalue_thresh=methylation_pvalue_thresh, count_threshold=count_threshold, benjamini_threshold=benjamini_threshold, save_to_eps=save_to_eps, folder_extension=folder_extension)
```

Alternatively using the CLI

```
metaKEGG --input_file_path="examples/single_input_genes.xlsx"  --input_label="input1"  --sheet_name_paths="pathways" --sheet_name_genes="gene_metrics" 
         --genes_column="gene_symbol" --log2fc_column="logFC" --analysis_type=4 
         --methylation_path="examples/methylation.csv" --methylation_gene="methylation_gene_symbol" --methylation_pvalue="methylation_pval" --methylation_pvalue_thresh=0.05 
         --count_threshold=2 --benjamini_threshold=None --save_to_eps=False --folder_extension=None
```

### Single Input Analysis with miRNA (Gene IDs)

1. Define input arguments

```
input_file_path = "examples/single_input_genes.xlsx"
input_label = "input1"
sheet_name_paths = "pathways"
sheet_name_genes = "gene_metrics"
genes_column = "gene_symbol"
log2fc_column = "logFC"
analysis_type = 5
miRNA_path = "examples/miRNA.tsv"
miRNA_gene = "miRNA_gene_symbol"
miRNA_pvalue = "miRNA_pval"
miRNA_pvalue_thresh=0.05,
count_threshold = 2
benjamini_threshold = None
save_to_eps = False
folder_extension = None
```

2. Run analysis

```
import metaKEGG

metaKEGG.Pipeline(input_file_path=input_file_path, input_label=input_label, sheet_name_paths=sheet_name_paths, sheet_name_genes=sheet_name_genes, genes_column=genes_column, log2fc_column=log2fc_column, analysis_type=analysis_type, miRNA_path=miRNA_path, miRNA_gene=miRNA_gene, miRNA_pvalue=miRNA_pvalue, miRNA_pvalue_thresh=miRNA_pvalue_thresh,count_threshold=count_threshold, benjamini_threshold=benjamini_threshold, save_to_eps=save_to_eps, folder_extension=folder_extension)
```

Alternatively using the CLI

```
metaKEGG --input_file_path="examples/single_input_genes.xlsx"  --input_label="input1"  --sheet_name_paths="pathways" --sheet_name_genes="gene_metrics" 
         --genes_column="gene_symbol" --log2fc_column="logFC" --analysis_type=5 
         --miRNA_path="examples/miRNA.tsv" --miRNA_gene="miRNA_gene_symbol" --miRNA_pvalue="miRNA_pval" --miRNA_pvalue_thresh=0.05 
         --count_threshold=2 --benjamini_threshold=None --save_to_eps=False --folder_extension=None
```

### Single Input Analysis with Methylation and miRNA (Gene IDs)

1. Define input arguments

```
input_file_path = "examples/single_input_genes.xlsx"
input_label = "input1"
sheet_name_paths = "pathways"
sheet_name_genes = "gene_metrics"
genes_column = "gene_symbol"
log2fc_column = "logFC"
analysis_type = 6
methylation_path = "examples/methylation.csv"
methylation_gene = "methylation_gene_symbol"
methylation_pvalue = "methylation_pval"
methylation_pvalue_thresh = 0.05
miRNA_path = "examples/miRNA.tsv"
miRNA_gene = "miRNA_gene_symbol"
miRNA_pvalue = "miRNA_pval"
miRNA_pvalue_thresh = 0.05
count_threshold = 2
benjamini_threshold = None
save_to_eps = False
folder_extension = None
```

2. Run analysis

```
import metaKEGG

metaKEGG.Pipeline(input_file_path=input_file_path, input_label=input_label, sheet_name_paths=sheet_name_paths, sheet_name_genes=sheet_name_genes, genes_column=genes_column, log2fc_column=log2fc_column, analysis_type=analysis_type, methylation_path=methylation_path, methylation_gene=methylation_gene, methylation_pvalue=methylation_pvalue, methylation_pvalue_thresh=methylation_pvalue_thresh, miRNA_path=miRNA_path, miRNA_gene=miRNA_gene, miRNA_pvalue=miRNA_pvalue, miRNA_pvalue_thresh=miRNA_pvalue_thresh, count_threshold=count_threshold, benjamini_threshold=benjamini_threshold, save_to_eps=save_to_eps, folder_extension=folder_extension)
```

Alternatively using the CLI

```
metaKEGG --input_file_path="examples/single_input_genes.xlsx"  --input_label="input1"  --sheet_name_paths="pathways" --sheet_name_genes="gene_metrics" 
         --genes_column="gene_symbol" --log2fc_column="logFC" --analysis_type=6 
         --methylation_path="examples/methylation.csv" --methylation_gene="methylation_gene_symbol" --methylation_pvalue="methylation_pval" --methylation_pvalue_thresh=0.05 
         --miRNA_path="examples/miRNA.tsv" --miRNA_gene="miRNA_gene_symbol" --miRNA_pvalue="miRNA_pval" --miRNA_pvalue_thresh=0.05 
         --count_threshold=2 --benjamini_threshold=None --save_to_eps=False --folder_extension=None
```

### Single input (Gene IDs) Bulk mapping

1. Define input arguments

```
input_file_path = "examples/single_input_bulk.xlsx"
input_label = "input1"
genes_column = "gene_symbol"
log2fc_column = "logFC"
sheet_name_paths = "pathways"
sheet_name_genes = "gene_metrics"
analysis_type = 7
count_threshold = None
benjamini_threshold = None
save_to_eps = False
folder_extension = None
```

2. Run analysis

```
import metaKEGG

metaKEGG.Pipeline(input_file_path=input_file_path, input_label=input_label, sheet_name_paths=sheet_name_paths, sheet_name_genes=sheet_name_genes, genes_column=genes_column, log2fc_column=log2fc_column, analysis_type=analysis_type, count_threshold=count_threshold, benjamini_threshold=benjamini_threshold, save_to_eps=save_to_eps, folder_extension=folder_extension)
```

Alternatively using the CLI

```
metaKEGG --input_file_path="examples/single_input_genes.xlsx"  --input_label="input1"  --sheet_name_paths="pathways" --sheet_name_genes="gene_metrics" --genes_column="gene_symbol" --log2fc_column="logFC" --analysis_type=7 --count_threshold=None --benjamini_threshold=None --save_to_eps=False --folder_extension=None
```

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
methylation_gene = "methylation_gene_symbol"
methylation_pvalue = "methylation_pval"
methylation_pvalue_thresh = 0.05
miRNA_path = "examples/miRNA.tsv"
miRNA_gene = "miRNA_gene_symbol"
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

m = metaKEGG.Pipeline(input_file_path=input_file_path, input_label=input_label, sheet_name_paths=sheet_name_paths, sheet_name_genes=sheet_name_genes, genes_column=genes_column, log2fc_column=log2fc_column, analysis_type=analysis_type, methylation_path=methylation_path, methylation_gene=methylation_gene, methylation_pvalue=methylation_pvalue, methylation_pvalue_thresh=methylation_pvalue_thresh, miRNA_path=miRNA_path, miRNA_gene=miRNA_gene, miRNA_pvalue=miRNA_pvalue, miRNA_pvalue_thresh=miRNA_pvalue_thresh, count_threshold=count_threshold, benjamini_threshold=benjamini_threshold, save_to_eps=save_to_eps, folder_extension=folder_extension)
```

```
m.single_input_genes()
m.single_input_with_methylation()
m.single_input_with_miRNA()
m.single_input_with_methylation_and_miRNA()
```

