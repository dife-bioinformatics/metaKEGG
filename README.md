[![Stars](https://img.shields.io/github/stars/dife-bioinformatics/metaKEGG?style=flat&logo=GitHub&color=yellow)](https://github.com/dife-bioinformatics/metaKEGG)
[![PyPI](https://img.shields.io/pypi/v/metaKEGG?logo=PyPI)](https://pypi.org/project/metaKEGG)


# `metaKEGG`

metaKEGG is a fully integrated solution with class-leading features to visualize the KEGG pathway enrichment analysis results from the DAVID Functional Annotation Tool, or RNAseq data.

## Table of Contents
- [Disclaimer](#disclaimer)
- [Web application](#web-application)
- [Installing metaKEGG](#installing-metakegg)
    - [Environment preparation](#environment-preparation)
    - [Install from PyPI](#install-from-pypi)
    - [Local installation with venv and requirement.txt in Windows](#local-installation-with-venv-and-requirementtxt-in-windows)
    - [Local installation with conda env and environment.yml](#local-installation-with-conda-env-and-environmentyml)
- [Getting started](#getting-started)
    - [Programmatic/Library usage](#programmaticlibrary-usage)
- [File structure](#file-structure)
    - [Main input file](#main-input-file)
    - [Metadata files](#metadata-files)
- [Example files](#example-files)
- [Parameters / Arguments](#parameters--arguments) 
- [Example usage](#example-usage)
    - [Gene expression](#gene-expression)
    - [Transcript expression](#transcript-expression)
    - [Bulk RNAseq mapping](#bulk-rnaseq-mapping)
    - [Multiple inputs](#multiple-inputs)
    - [Differentially methylated genes](#differentially-methylated-genes)
    - [DMPs per gene](#dmps-per-genes)
        - [Probe correction](#probe-correction)
    - [miRNA target genes](#mirna-target-genes)
    - [DEmiRs per gene](#demirs-per-gene)    
    - [Methylated + miRNA target genes](#methylated--mirna-target-genes)    

## Disclaimer

metaKEGG uses the KEGG API, and is restricted to academic use by academic users belonging to academic institutions.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.
You should have received a copy of the GNU Affero General Public License along with this
program. If not, see [GNU Affero General Public License](https://www.gnu.org/licenses/#AGPL).

Author: Michail Lazaratos, Deutsches Institut für Ernährungsforschung Potsdam-Rehbrücke / German Institute of Human Nutrition Potsdam-Rehbrücke (DIfE)

## Web application

You can perform all analyses pipelines graphically through our [metaKEGG web application](https://metakegg.apps.dzd-ev.org/).

You can check out the metaKEGG web-wrapper project [here](https://github.com/DZD-eV-Diabetes-Research/meta-kegg-web-wrapper).



## Installing metaKEGG

Either clone from this GitHub repository or install from 
[PyPI](https://pypi.org/)


### Environment preparation

1. Install an [Anaconda Distribution](https://docs.anaconda.com/free/anaconda/install/) or [Miniconda](https://docs.anaconda.com/free/miniconda/miniconda-other-installer-links/)


2. Create a new conda environment to install metaKEGG in

Example using python 3.11. The pakcage was tested with versions 3.9 through 3.12. Later versions should work but there is no guarantee.

```
conda create -n metaKEGG_env python=3.11
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

After successfully installing metaKEGG, you can use it by importing the library.

### Programmatic/Library usage

The Pipeline class requires global arguments for instantiation.

```
import metaKEGG
metaKEGG.Pipeline(**kwargs)
```

Alternatively
```
from metaKEGG import Pipeline
Pipeline(**kwargs)
```

---
## File structure

metaKEGG requires the input files to be structure in a certain format before submitted for analysis.

### Main input file

The pathway enrichment results (DAVID format) should be in a Microsoft Excel file.
The sheeting containging the pathway enrichment results must contain the columns :

- **`Category`**: Enirichment analysis identifier. Must contain keyword *KEGG_PATHWAY* to be processed 
- **`Term`**: KEGG pathway ID and full name format  `organism code number:name`. Example  `hsa04932:Non-alcoholic fatty liver disease`
- **`Count`**: Number of genes in the pathway
- **`PValue`**: Raw pathway p-value 
- **`Genes`**: Genes in pathway seperated by `, ` (a comma followed by a whitespace)
- **`Benjamini`**: Adjusted p-values using the Benjamini-Hochberg method
 
 <br>

| Category     | Term              | Count | PValue  | Genes                | Benjamini |
|--------------|-------------------|-------|---------|----------------------|-----------|
| KEGG_PATHWAY | mmu00000:Pathway1 | 3     | 0.003   | Gene1, Gene2, Gene3  | 0.04      |
| KEGG_PATHWAY | mmu99999:Pathway2 | 3     | 0.02    | Gene1, Gene4, Gene5  | 0.06      |
<br>

In the same file, a different sheet contains the gene/transcript information and has a minimum of two columns, the official gene symbol and the log<sub>2</sub>FC values. Optionally, a column with gene IDs can be also present.

| ID                 | gene_symbol | log<sub>2</sub>FC |
|--------------------|-------------|-------------------|
| ENSMUSG00000000001 | Gene1       |-2.45              |
| ENSMUSG00000000002 | Gene2       | 1.33              |
| ENSMUSG00000000003 | Gene3       | 0.44              |
| ENSMUSG00000000004 | Gene4       |-1.16              |
| ENSMUSG00000000005 | Gene5       | 2.01              |
<br>

**_NOTE:_**  In the case of `Bulk RNAseq mapping` analysis, the pathway file should only contain the `Term` column.

| Term              |
|-------------------|
| mmu00000:Pathway1 |
| mmu99999:Pathway2 |
<br>

---

### Metadata files

The accepted file types are:

- **Excel (.xls or .xlsx)**: The file must contain exactly one sheet.
- **TSV (.tsv)**: Tab-separated value files.
- **CSV (.csv)**: Comma-separated value files.

Metadata files contain a minimum of three columns:

- **ID**
- **Gene Symbol**
- **P-value**

<br>

Example for methylation metadata
<br>

| CG_ID  | methylation_gene_symbol  | methylation_pval |
|--------|--------------------------|------------------|
| cg_1   | Gene1                    | 0.003            | 
| cg_2   | Gene2                    | 0.7              | 
| cg_3   | Gene3                    | 0.04             | 
<br>

Example for miRNA metadata
<br>

| miRNA_ID      | miRNA_gene_symbol  | miRNA_pval  |
|---------------|--------------------|-------------|
| mmu-miR-1-3p  | Gene1              | 0.01        | 
| mmu-miR-2-3p  | Gene2              | 0.03        | 
| mmu-miR-3-3p  | Gene3              | 0.8         | 
<br>

---

### Example Files

In the `/examples` directory, you can find example files to perform all analysis types provided by **metaKEGG**. The paths to these files will be used for demonstration purposes below.

---

### Analysis Types Available

1. **Gene expression**
2. **Transcript expression**  
3. **Bulk RNAseq mapping**  
4. **Multiple inputs**  
5. **Differentially methylated genes**  
6. **DMPs per gene**  
7. **miRNA target genes**  
8. **DEmiRs per gene**  
9. **Methylated + miRNA target genes**  

---

## Parameters / Arguments

### General Input Parameters
These parameters are applicable for all analysis types:

- **`input_file_path`**: Path to enrichment analysis and gene information, or bulk RNAseq & selected pathways (Excel format).  
  *(No default value, required)*  
- **`sheet_name_paths`**: Sheet name containing the pathway information.  
  *(Default: `"pathways"`)*
- **`sheet_name_genes`**: Sheet name for gene information.  
  *(Default: `"gene_metrics"`)*
- **`genes_column`**: Column name for gene symbols in the `sheet_name_genes`.  
  *(Default: `"gene_symbol"`)*
- **`log2fc_column`**: Column name for log2fc values in the `sheet_name_genes`.  
  *(Default: `"logFC"`)*
- **`save_to_eps`**: A boolean (`True`/`False`) indicating whether to save the maps and color scales or legends as separate `.eps` files in addition to the `.pdf` exports.  
  *(Default: `False`)*
- **`folder_extension`**: Folder extension to be appended to the default naming scheme. If `None` and the default folder exists, the folder will be overwritten.  
  *(Default: `None`)*
- **`output_folder_name`**: Folder name that will overwrite the default naming scheme. If `None`, the default naming scheme will be used.  
  *(Default: `None`)*
- **`input_label`**: Label of input file (Required for multiple inputs to create the legend).  
  *(No default value, required for multiple inputs analysis)*

---

### Analysis-Specific Parameters

#### Methylation Analysis Parameters
These parameters are specific to methylation-based analysis:

- **`methylation_file_path`**: Path to methylation data (Excel, CSV, or TSV format).  
  *(No default value, required for methylation-related analyses)*  
- **`methylation_genes_column`**: Column name for methylation gene symbols.  
  *(No default value, required for methylation-related analyses)*  
- **`methylation_probe_column`**: Column name for methylation probes.  
  *(No default value, required for 6: DMPs per Gene)*  
- **`methylation_pvalue_column`**: Column name for methylation p-values.  
  *(Default: `None`)*
- **`methylation_pvalue_threshold`**: P-value threshold for the probes.  
  *(Default: `0.05`)*
- **`probes_to_cgs`**: If `True`, will perform probe correction. See [Probe correction](#probe-correction).  
  *(Default: `False`)*

#### miRNA Analysis Parameters
These parameters are specific to miRNA-based analysis:

- **`miRNA_file_path`**: Path to miRNA data (Excel, CSV, or TSV format).  
  *(No default value, required for miRNA-related analyses)*  
- **`miRNA_genes_column`**: Column name for miRNA gene symbols.  
  *(No default value, required for miRNA-related analyses)*  
- **`miRNA_ID_column`**: Column name for miRNA IDs.  
  *(No default value, required for 8: DEmiRs per Gene)*  
- **`miRNA_pvalue_column`**: Column name for miRNA p-value.  
  *(Default: `None`)*
- **`miRNA_pvalue_threshold`**: P-value threshold for miRNAs.  
  *(Default: `0.05`)*

-----

### Filter Parameters
These parameters are used to filter the data for analysis:

- **`count_threshold`**: Minimum number of genes per pathway for pathway to be drawn.  
  *(Default: `2`)*
- **`pathway_pvalue_threshold`**: Raw p-value threshold for pathways to be drawn.  
  *(Default: `None`)*
- **`benjamini_threshold`**: Benjamini Hochberg p-value threshold for the pathways.  
  *(Default: `None`)*

---

## Example Usage

### Gene expression

This function takes a single input file as an argument and maps the detected genes on the enriched KEGG reference pathway, and colors them according to their log2FC values.

1. Define input arguments

```
input_file_path = "examples/single_input_genes.xlsx"
sheet_name_paths = "pathways"
sheet_name_genes = "gene_metrics"
genes_column = "gene_symbol"
log2fc_column = "logFC"
count_threshold = 2
pathway_pvalue_threshold = None
benjamini_threshold = None
save_to_eps = False
folder_extension = None
```

2. Run analysis

```
import metaKEGG

my_pipeline = metaKEGG.Pipeline(input_file_path=input_file_path,
                                sheet_name_paths=sheet_name_paths,
                                sheet_name_genes=sheet_name_genes,
                                genes_column=genes_column,
                                log2fc_column=log2fc_column,
                                save_to_eps=save_to_eps,
                                folder_extension=folder_extension)

my_pipeline.gene_expression(benjamini_threshold=benjamini_threshold,
                            count_threshold=count_threshold,
                            pathway_pvalue_threshold=pathway_pvalue_threshold)

```

-----
### Transcript expression

This function takes a single input file as an argument and maps the detected transcripts on the enriched KEGG reference pathway, and colors them according to their log2FC values.

**_NOTE:_**  Pathway enrichement analysis with the DAVID Functional Annotation Tool, should be performed using transcript IDs.

1. Define input arguments

```
input_file_path = "examples/single_input_transcripts.xlsx"
sheet_name_paths = "pathways"
sheet_name_genes = "transcript_metrics"
genes_column = "gene_symbol"
log2fc_column = "logFC"
count_threshold = 2
pathway_pvalue_threshold = None
benjamini_threshold = None
save_to_eps = False
folder_extension = None
```

2. Run analysis

```
import metaKEGG

my_pipeline = metaKEGG.Pipeline(input_file_path=input_file_path,
                                sheet_name_paths=sheet_name_paths,
                                sheet_name_genes=sheet_name_genes,
                                genes_column=genes_column,
                                log2fc_column=log2fc_column,
                                save_to_eps=save_to_eps,
                                folder_extension=folder_extension)

my_pipeline.transcript_expression(benjamini_threshold=benjamini_threshold,
                                  count_threshold=count_threshold,
                                  pathway_pvalue_threshold=pathway_pvalue_threshold)
```
-----
### Bulk RNAseq mapping

This function takes RANseq data, as single input file argument, maps the genes on a provided list of target pathways (assuming they are also found in the target pathways), and colors them according to their log2FC values.

1. Define input arguments

```
input_file_path = "examples/single_input_bulk.xlsx"
genes_column = "gene_symbol"
log2fc_column = "logFC"
sheet_name_paths = "pathways"
sheet_name_genes = "gene_metrics"
save_to_eps = False
folder_extension = None
```

2. Run analysis

```
import metaKEGG

my_pipeline = metaKEGG.Pipeline(input_file_path=input_file_path,
                                sheet_name_paths=sheet_name_paths,
                                sheet_name_genes=sheet_name_genes,
                                genes_column = genes_column,
                                log2fc_column=log2fc_column,
                                save_to_eps=save_to_eps)

my_pipeline.bulk_rnaseq_mapping()
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
count_threshold = 2
pathway_pvalue_threshold = None
benjamini_threshold = None
save_to_eps = False
folder_extension = None
```

2. Run analysis

```
import metaKEGG

my_pipeline = metaKEGG.Pipeline(input_file_path=input_file_path,
                                sheet_name_paths=sheet_name_paths,
                                sheet_name_genes=sheet_name_genes,
                                genes_column=genes_column,
                                log2fc_column=log2fc_column,
                                save_to_eps=save_to_eps)
                        
my_pipeline.multiple_inputs(input_label=input_label,
                            count_threshold=count_threshold,
                            benjamini_threshold=benjamini_threshold,
                            pathway_pvalue_threshold=pathway_pvalue_threshold)

```
-----
### Differentially methylated genes

This function takes a single input file and a methylation metadata file as arguments and maps the detected genes on the enriched KEGG reference pathway, and colors them according to their methylation state. The state is defined as a binary reprsentation, depending if DMPs corresponding to a given gene are detected, or not.


1. Define input arguments

```
input_file_path = "examples/single_input_genes.xlsx"
sheet_name_paths = "pathways"
sheet_name_genes = "gene_metrics"
genes_column = "gene_symbol"
log2fc_column = "logFC"
methylation_file_path = "examples/methylation.csv"
methylation_genes_column = "methylation_gene_symbol"
methylation_pvalue_column = "methylation_pval"
methylation_pvalue_threshold = 0.05
count_threshold = 2
pathway_pvalue_threshold = None
benjamini_threshold = None
save_to_eps = False
folder_extension = None
```

2. Run analysis

```
import metaKEGG

my_pipeline = metaKEGG.Pipeline(input_file_path=input_file_path,
                                sheet_name_paths=sheet_name_paths,
                                sheet_name_genes=sheet_name_genes,
                                genes_column=genes_column,
                                log2fc_column=log2fc_column,
                                save_to_eps=save_to_eps)

my_pipeline.methylated_genes(methylation_file_path=methylation_file_path,
                             methylation_genes_column=methylation_genes_column,
                             methylation_pvalue_column=methylation_pvalue_column,
                             methylation_pvalue_threshold=methylation_pvalue_threshold,
                             count_threshold=count_threshold,
                             benjamini_threshold=benjamini_threshold,
                             pathway_pvalue_threshold=pathway_pvalue_threshold)
```
-----
### DMPs per gene

This function takes a single input file and a methylation metadata file as arguments and maps the detected genes on the enriched KEGG reference pathway. It generates bins to quantify the number of DMPs that correspond to a given gene, and colors a gege according to its DMP bin. The function also returns the quantification histogram plots, both in a grouped and an absolute count representation.

1. Define input arguments

```
input_file_path = "examples/single_input_genes.xlsx"
sheet_name_paths = "pathways"
sheet_name_genes = "gene_metrics"
genes_column = "gene_symbol"
log2fc_column = "logFC"
methylation_file_path = "examples/methylation_for_quantification.csv"
methylation_genes_column = "methylation_gene_symbol"
methylation_pvalue_column = "methylation_pval"
methylation_pvalue_threshold = 0.05
methylation_probe_column = "CG_ID"
probes_to_cgs=False
count_threshold = 2
pathway_pvalue_threshold = None
benjamini_threshold = None
save_to_eps = False
folder_extension = None
```

2. Run analysis

```
import metaKEGG

my_pipeline = metaKEGG.Pipeline(input_file_path=input_file_path,
                                sheet_name_paths=sheet_name_paths,
                                sheet_name_genes=sheet_name_genes,
                                genes_column=genes_column,
                                log2fc_column=log2fc_column,
                                save_to_eps=save_to_eps)

my_pipeline.dmps_per_gene(methylation_file_path=methylation_file_path,
                          methylation_genes_column=methylation_genes_column,
                          methylation_pvalue_column=methylation_pvalue_column,
                          methylation_pvalue_threshold=methylation_pvalue_threshold,
                          methylation_probe_column=methylation_probe_column,
                          probes_to_cgs=probes_to_cgs,
                          count_threshold=count_threshold,
                          benjamini_threshold=benjamini_threshold,
                          pathway_pvalue_threshold=pathway_pvalue_threshold)
```

### Probe correction

**_NOTE:_**  When using probes_to_cgs=True, the pipeline will split the CG probes by the underscore '_' character and keep the first part, essentially correcting for different probe chemistry that could occur in the same position. Example format is cg00000000_BC21 and cg00000000_TC21, which would be counted as two separate probes targeting the same gene. Using the argument probes_to_cgs with True, the probes become cg00000000 & cg00000000, and duplicated entries per gene are eliminated, essentially counting one probe for the target gene.

------
### miRNA target genes

This function takes a single input file and a miRNA metadata file as arguments and maps the detected genes on the enriched KEGG reference pathway, and colors them according to their miRNA state. The state is defined as a binary reprsentation, depending if DEmiRs are targeting a given gene, or not.

1. Define input arguments

```
input_file_path = "examples/single_input_genes.xlsx"
sheet_name_paths = "pathways"
sheet_name_genes = "gene_metrics"
genes_column = "gene_symbol"
log2fc_column = "logFC"
miRNA_file_path = "examples/miRNA.tsv"
miRNA_genes_column = "miRNA_gene_symbol"
miRNA_pvalue_column = "miRNA_pval"
miRNA_pvalue_threshold=0.05
pathway_pvalue_threshold = None
count_threshold = 2
benjamini_threshold = None
save_to_eps = False
folder_extension = None
```

2. Run analysis

```
import metaKEGG

my_pipeline = metaKEGG.Pipeline(input_file_path=input_file_path,
                                sheet_name_paths=sheet_name_paths,
                                sheet_name_genes=sheet_name_genes,
                                genes_column=genes_column,
                                log2fc_column=log2fc_column,
                                save_to_eps=save_to_eps)

my_pipeline.mirna_target_genes(miRNA_file_path=miRNA_file_path,
                               miRNA_genes_column=miRNA_genes_column,
                               miRNA_pvalue_column=miRNA_pvalue_column,
                               miRNA_pvalue_threshold=miRNA_pvalue_threshold,
                               count_threshold=count_threshold,
                               benjamini_threshold=benjamini_threshold,
                               pathway_pvalue_threshold=pathway_pvalue_threshold)
```
------
### DEmiRs per gene

This function takes a single input file and a miRNA metadata file as arguments and maps the detected genes on the enriched KEGG reference pathway. It generates bins to quantify the number of DEmiRs that correspond to a given gene, and colors a gege according to its DEmiR bin. The function also returns the quantification histogram plots, both in a grouped and an absolute count representation.

1. Define input arguments

```
input_file_path = "examples/single_input_genes.xlsx"
sheet_name_paths = "pathways"
sheet_name_genes = "gene_metrics"
genes_column = "gene_symbol"
log2fc_column = "logFC"
miRNA_file_path = "examples/miRNA_for_quantification.tsv"
miRNA_genes_column = "miRNA_gene_symbol"
miRNA_pvalue_column = "miRNA_pval"
miRNA_pvalue_threshold = 0.05
miRNA_column = "miRNA_ID"
pathway_pvalue_threshold = None
count_threshold = 2
benjamini_threshold = None
save_to_eps = False
folder_extension = None
```

2. Run analysis

```
import metaKEGG

my_pipeline = metaKEGG.Pipeline(input_file_path=input_file_path,
                                sheet_name_paths=sheet_name_paths,
                                sheet_name_genes=sheet_name_genes,
                                genes_column=genes_column,
                                log2fc_column=log2fc_column,
                                save_to_eps=save_to_eps)

my_pipeline.demirs_per_gene(miRNA_file_path=miRNA_file_path,
                            miRNA_genes_column=miRNA_genes_column,
                            miRNA_pvalue_column=miRNA_pvalue_column,
                            miRNA_pvalue_threshold=miRNA_pvalue_threshold,
                            miRNA_ID_column=miRNA_column, 
                            count_threshold=count_threshold,
                            benjamini_threshold=benjamini_threshold,
                            pathway_pvalue_threshold=pathway_pvalue_threshold)
```
-------
### Methylated + miRNA target genes

This function takes a single input file, a methylation, and a miRNA metadata file as arguments and maps the detected genes on the enriched KEGG reference pathway. Genes are colored according to their methylation and miRNA states. The states is defined as a binary reprsentations of the methylation and miRNA combinations.

1. Define input arguments

```
input_file_path = "examples/single_input_genes.xlsx"
sheet_name_paths = "pathways"
sheet_name_genes = "gene_metrics"
genes_column = "gene_symbol"
log2fc_column = "logFC"
methylation_file_path = "examples/methylation.csv"
methylation_genes_column = "methylation_gene_symbol"
methylation_pvalue_column = "methylation_pval"
methylation_pvalue_threshold = 0.05
miRNA_file_path = "examples/miRNA.tsv"
miRNA_genes_column = "miRNA_gene_symbol"
miRNA_pvalue_column = "miRNA_pval"
miRNA_pvalue_threshold = 0.05
count_threshold = 2
pathway_pvalue_threshold = None
benjamini_threshold = None
save_to_eps = False
folder_extension = None
```

2. Run analysis

```
import metaKEGG

my_pipeline = metaKEGG.Pipeline(input_file_path=input_file_path,
                        sheet_name_paths=sheet_name_paths,
                        sheet_name_genes=sheet_name_genes,
                        genes_column=genes_column,
                        log2fc_column=log2fc_column,
                        save_to_eps=save_to_eps)

my_pipeline.methylated_and_mirna_target_genes(methylation_file_path=methylation_file_path,
                                              methylation_genes_column=methylation_genes_column,
                                              methylation_pvalue_column=methylation_pvalue_column,
                                              methylation_pvalue_threshold=methylation_pvalue_threshold,
                                              miRNA_file_path=miRNA_file_path,
                                              miRNA_genes_column=miRNA_genes_column,
                                              miRNA_pvalue_column=miRNA_pvalue_column,
                                              miRNA_pvalue_threshold=miRNA_pvalue_threshold,
                                              count_threshold=count_threshold,
                                              benjamini_threshold=benjamini_threshold,
                                              pathway_pvalue_threshold=pathway_pvalue_threshold)
```
-----