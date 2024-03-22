import numpy as np
from Bio.KEGG import REST
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import pandas as pd
import os
import time
from Bio.KEGG.REST import *
from matplotlib.backends.backend_pdf import PdfPages
import PyPDF2
import requests
import PIL.Image
from io import BytesIO
from reportlab.pdfgen import canvas
import os
import datetime
import shutil
from matplotlib.patches import Patch
import sys
from pylab import *
import pathlib

import ssl
ssl._create_default_https_context = ssl._create_unverified_context

from ..config import csv_suffixes, tsv_suffixes, excel_suffixes

def get_colors_from_colorscale(colorscale_list, how_many=None):
    """
    Generate a list of hex color codes from a given colorscale list.

    Parameters:
    - colorscale_list (list): A list of colorscale names or colormaps.
    - how_many (int, optional): The number of colors to generate for each colorscale. If not specified,
      the default behavior is to generate the full range of colors in each colorscale.

    Returns:
    list: A list of hex color codes, representing the colors in the specified colorscales.

    Example:
    >>> colorscales = ['viridis', 'plasma']
    >>> get_colors_from_colorscale(colorscales, how_many=5)
    ['#440154', '#3b528b', '#21918c', '#5ec962', '#fde725', '#0d0887', '#7e03a8', '#cc4778', '#f89540', '#f0f921']
    """
    color_list_out = []
    cmap = None

    if how_many:
        for cscale in colorscale_list:
            cmap = cm.get_cmap(cscale, how_many)
            color_list = [matplotlib.colors.rgb2hex(cmap(i)[:3]) for i in range(cmap.N)]
            color_list_out.extend(color_list)
    else:
        for cscale in colorscale_list:
            cmap = cm.get_cmap(cscale)
            color_list = [matplotlib.colors.rgb2hex(cmap(i)[:3]) for i in range(cmap.N)]
            color_list_out.extend(color_list)

    return color_list_out

def file_naming_scheme(input_data, id=None):
    """
    Generate a truncated and standardized file name based on the provided ID and parsed output.

    Parameters:
    - id (str or int): The identifier used to access information in the parsed output.
    - parsed_output (dict): A dictionary containing parsed information, typically obtained from some data processing.

    Returns:
    str: A truncated and standardized file name based on the specified ID and parsed output.

    Example:
    >>> parsed_data = {mmu04141: {'name': 'Protein processing in endoplasmic reticulum'} , mmu00010: {'name': 'Glycolysis / Gluconeogenesis'}}
    >>> file_naming_scheme('mmu04141', parsed_data)
    'Protein_processing_in_endoplasmic_reticulum'
    
    >>> file_naming_scheme('mmu00010', parsed_data)
    'Glycolisis_Gluconeogenesis'
    """
    char_length = 35

    if isinstance(input_data, dict):
        if id is None:
            raise ValueError("When input_data is a dictionary, a valid 'id' must be provided.")
        pathway_name = input_data[id]['name']
    elif isinstance(input_data, str):
        pathway_name = input_data
    else:
        raise ValueError("Unsupported input type. Use a dictionary or a string.")

    truncated_name = pathway_name[:char_length]

    if '/' in truncated_name:
        truncated_name = truncated_name.replace('/', '')
        
    if len(pathway_name) > char_length:
        if truncated_name[-1] != ' ':
            truncated_name = ' '.join(truncated_name.split(' ')[:-1])
            truncated_name = truncated_name.replace(' ', '_')
        else:
            truncated_name = truncated_name[:-1]
            truncated_name = truncated_name.replace(' ', '_')
    else:
        truncated_name = truncated_name.replace(' ', '_')

    return truncated_name

def collect_pathway_info(parsed_output):
    """
    Collect information about genes, KO identifiers, and related data for each pathway in the parsed output.

    Parameters:
    - parsed_output (dict): A dictionary containing parsed information, typically obtained from some data processing.

    Returns:
    dict: A dictionary with pathway identifiers as keys, and information about genes, KO identifiers, and related data
          as values.

    Example:
    >>> parsed_data = {'pathway1': {'name': 'Pathway 1'}, 'pathway2': {'name': 'Pathway 2'}}
    >>> collect_pathway_info(parsed_data)
    {'pathway1': {'gene_symbol_KO': {'Gene1': 'KO123', 'Gene2': 'KO456'},
                  'gene_symbols': ['Gene1', 'Gene2'],
                  'gene_symbol_kegg_id': {'Gene1': 'id123', 'Gene2': 'id456'},
                  'corresponding_KO': 'KO789'},
     'pathway2': {'gene_symbol_KO': {'Gene3': 'KO987', 'Gene4': 'KO654'},
                  'gene_symbols': ['Gene3', 'Gene4'],
                  'gene_symbol_kegg_id': {'Gene3': 'id987', 'Gene4': 'id654'},
                  'corresponding_KO': 'KO321'}}
    """
    pathway_genes = {}
    for pathway in parsed_output.keys():
        gene_symbols = []
        gene_ids = []
        gene_KO = []

        gene_symbol_KO = {}
        gene_KO_symbol = {}
        gene_symbol_kegg_id = {}

        pathway_file = REST.kegg_get(pathway).read()

        current_section = None
        for line in pathway_file.rstrip().split("\n"):
            section = line[:12].strip()
            if not section == "":
                current_section = section
            if current_section == "KO_PATHWAY":
                corresponding_KO = line[12:]
            if current_section == "GENE":
                try:
                    gene_identifiers, gene_description = line[12:].split("; ")
                    gene_id, gene_symbol = gene_identifiers.split()
                except ValueError:
                    continue

                find_KO = gene_description.split('[')
                for substring in find_KO:
                    if ']' in substring and substring.startswith('KO'):
                        gene_KO.append(substring.split(']')[0])
                        gene_symbol_KO[gene_symbol] = substring.split(']')[0].split(':')[1]
                        gene_symbol_KO[gene_symbol]

                        if substring.split(']')[0].split(':')[1] not in gene_KO_symbol:
                            gene_KO_symbol[substring.split(']')[0].split(':')[1]] = gene_symbol

                if not gene_symbol in gene_symbols:
                    gene_symbols.append(gene_symbol)
                if not gene_id in gene_ids:
                    gene_ids.append(gene_id)
                gene_symbol_kegg_id[gene_symbol] = gene_id


        pathway_genes[pathway] = {"gene_symbol_KO": gene_symbol_KO,
                                "gene_symbols": gene_symbols,
                                "gene_symbol_kegg_id": gene_symbol_kegg_id,
                                "corresponding_KO" : corresponding_KO}
        time.sleep(1)

    return pathway_genes

def collect_pathway_info_multiple_interventions(pathway_id):
    """
    Collect information about genes, KO identifiers, and related data for a specific pathway with multiple interventions.

    Parameters:
    - pathway_id (str): The identifier of the pathway for which information needs to be collected.

    Returns:
    dict: A dictionary with the specified pathway identifier as the key and information about genes, KO identifiers,
          and related data as the value.

    Example:
    >>> pathway_id = 'pathway123'
    >>> collect_pathway_info_multiple_interventions(pathway_id)
    {'pathway123': {'gene_symbol_KO': {'Gene1': 'KO123', 'Gene2': 'KO456'},
                    'gene_symbols': ['Gene1', 'Gene2'],
                    'gene_symbol_kegg_id': {'Gene1': 'id123', 'Gene2': 'id456'},
                    'corresponding_KO': 'KO789'}}
    """
    pathway_genes = {}
    gene_symbols = []
    gene_ids = []
    gene_KO = []

    gene_symbol_KO = {}
    gene_KO_symbol = {}
    gene_symbol_kegg_id = {}

    pathway_file = REST.kegg_get(pathway_id).read()

    current_section = None
    for line in pathway_file.rstrip().split("\n"):
        section = line[:12].strip()
        if not section == "":
            current_section = section

        if current_section == "KO_PATHWAY":
            corresponding_KO = line[12:]

        if current_section == "GENE":
            try:
                gene_identifiers, gene_description = line[12:].split("; ")
                gene_id, gene_symbol = gene_identifiers.split()
            except ValueError:
                continue

            find_KO = gene_description.split('[')
            for substring in find_KO:
                if ']' in substring and substring.startswith('KO'):
                    gene_KO.append(substring.split(']')[0])
                    gene_symbol_KO[gene_symbol] = substring.split(']')[0].split(':')[1]
                    gene_symbol_KO[gene_symbol]

                    if substring.split(']')[0].split(':')[1] not in gene_KO_symbol:
                        gene_KO_symbol[substring.split(']')[0].split(':')[1]] = gene_symbol

            if not gene_symbol in gene_symbols:
                gene_symbols.append(gene_symbol)
            if not gene_id in gene_ids:
                gene_ids.append(gene_id)
            gene_symbol_kegg_id[gene_symbol] = gene_id


    pathway_genes[pathway_id] = {"gene_symbol_KO": gene_symbol_KO,
                            "gene_symbols": gene_symbols,
                            "gene_symbol_kegg_id": gene_symbol_kegg_id,
                            "corresponding_KO" : corresponding_KO}
    time.sleep(1)

    return pathway_genes

def create_output_folder(path , folder_extension=None):
    """
    Create an output folder with the specified path, handling existing directories and providing options for overwrite.

    Parameters:
    - path (str): The desired path for the output folder.
    - folder_extension (str, optional): An optional extension to be added to the folder name.

    Returns:
    str: The path of the created or modified output folder.

    Examples:
    >>> create_output_folder('/path/to/output_folder')
    New directory '/path/to/output_folder' created in '/current/working/directory'.
    Returns: '/path/to/output_folder'
    
    >>> create_output_folder('/path/to/output_folder' , folder_extension=None)
    The directory '/path/to/output_folder' already exists in '/current/working/directory'.
    The contents of '/path/to/output_folder' have been removed, and the directory has been recreated.
    Returns: '/path/to/output_folder'

    >>> create_output_folder('/path/to/output_folder', folder_extension='_backup')
    The directory '/path/to/output_folder' already exists in '/current/working/directory'.
    New directory '/path/to/output_folder_backup' created in '/current/working/directory'.
    Returns: '/path/to/output_folder_backup'
    """

    if folder_extension is not None:
        path = path + "_" + folder_extension
        if os.path.exists(path):
            shutil.rmtree(path)
            os.makedirs(path, mode=0o755)
            print(f"The contents of '{path}' have been removed, and the directory has been recreated.", '\n') 
        else:    
            os.makedirs(path, mode=0o755)
            print(f"New directory '{path}' created with extension {folder_extension} ", '\n')

    else:
        if not os.path.exists(path):
            os.makedirs(path, mode=0o755)
            print(f"New directory '{path}' created ", '\n')
        else:
            shutil.rmtree(path)
            os.makedirs(path, mode=0o755)
            print(f"The contents of '{path}' have been removed, and the directory has been recreated.", '\n')                    

    return path

def compile_and_write_output_files(id, pathway_id, output_name, color_legend=None, cmap=None, vmin=None, vmax=None, save_to_eps=False):
    """
    Compile and write output files for a given pathway ID.

    Parameters:
    - id (str): Identifier for the pathway.
    - pathway_id (str): Identifier for the pathway image.
    - output_name (str): Name for the output files.
    - color_legend (dict, optional): Dictionary mapping status labels to colors for legend version.
    - cmap (str, optional): Colormap for colorbar version.
    - vmin (float, optional): Minimum value for color normalization.
    - vmax (float, optional): Maximum value for color normalization.
    - save_to_eps (bool, optional): Flag to save colorbar version to EPS format.

    Returns:
    None

    Raises:
    ValueError: If color_legend is provided along with cmap, vmin, and vmax.

    Note:
    - The function generates PDF files based on the provided parameters, combining pathway image, colorbar, and legend versions.

    Example:
    >>> compile_and_write_output_files(id=id, pathway_id=pathway_id , cmap=cmap , vmin=vmin , vmax=vmax , output_name=output_name , save_to_eps=save_to_eps)

    """
    if cmap is not None and vmin is not None and vmax is not None and color_legend is None:
        # Use the colorbar version
        vertical_figsize = (1, 12)
        horizontal_figsize = (12, 1)
        cbar_width_ratios = [1, 6]
        cbar_height_ratios = [6, 1]

        with PdfPages(pathway_id + "_colorbar.pdf") as pdf:

            fig_v, ax_v = plt.subplots(figsize=vertical_figsize)
            sm_v = ScalarMappable(cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax))
            cb_v = plt.colorbar(sm_v, cax=ax_v, orientation='vertical', shrink=1, aspect=20, pad=0)
            cb_v.ax.set_ylabel(r'$\mathrm{log}_2(\mathrm{FC})$', rotation=270, labelpad=20, fontsize=20)
            cb_v.ax.tick_params(labelsize=16)

            if save_to_eps:
                plt.savefig(f'{id}_{output_name}_colorbar.eps', bbox_inches='tight', pad_inches=0)

            fig_h, ax_h = plt.subplots(figsize=horizontal_figsize)
            sm_h = ScalarMappable(cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax))
            cb_h = plt.colorbar(sm_h, cax=ax_h, orientation='horizontal', shrink=1, aspect=20)
            cb_h.ax.set_xlabel(r'$\mathrm{log}_2(\mathrm{FC})$', labelpad=10, fontsize=20)
            cb_h.ax.tick_params(labelsize=16)

            fig, axs = plt.subplots(nrows=2, ncols=2, gridspec_kw={'width_ratios': cbar_width_ratios, 'height_ratios': cbar_height_ratios})
            axs[0, 1].remove()
            axs[1, 0].remove()
            fig.tight_layout()

            pdf.savefig(fig_v, bbox_inches='tight' , dpi=300)
            pdf.savefig(fig_h, bbox_inches='tight' , dpi=300)

            plt.close(fig)            
            plt.close(fig_v)            
            plt.close(fig_h)            
            plt.close()

        url = f'http://rest.kegg.jp/get/{id}/image'
        response = requests.get(url)
        image = PIL.Image.open(BytesIO(response.content))
        pdf_file_name = f'{id}.pdf'
        pdf_canvas = canvas.Canvas(pdf_file_name, pagesize=image.size)
        temp_image_file = f'{id}.png'
        image.save(temp_image_file)
        pdf_canvas.drawImage(temp_image_file, 0, 0)
        pdf_canvas.save()
        os.remove(temp_image_file)

        pdf1 = open(pathway_id + ".pdf", 'rb')
        pdf_reader1 = PyPDF2.PdfReader(pdf1)
        pdf2 = open(id + ".pdf", 'rb')
        pdf_reader2 = PyPDF2.PdfReader(pdf2)
        pdf3 = open(pathway_id + "_colorbar.pdf", 'rb')
        pdf_reader3 = PyPDF2.PdfReader(pdf3)
        pdf_writer = PyPDF2.PdfWriter()

        for page_num in range(len(pdf_reader1.pages)):
            page = pdf_reader1.pages[page_num]
            pdf_writer.add_page(page)
        for page_num in range(len(pdf_reader2.pages)):
            page = pdf_reader2.pages[page_num]
            pdf_writer.add_page(page)
        for page_num in range(len(pdf_reader3.pages)):
            page = pdf_reader3.pages[page_num]
            pdf_writer.add_page(page)

        output_file = open(id + "_" + output_name + '.pdf', 'wb')
        pdf_writer.write(output_file)

        pdf1.close()
        pdf2.close()
        pdf3.close()
        output_file.close()
        os.remove(pathway_id + ".pdf")
        os.remove(id + ".pdf")
        os.remove(pathway_id + "_colorbar.pdf")

    elif color_legend is not None:
        # Use the legend version
        url = f'http://rest.kegg.jp/get/{id}/image'
        response = requests.get(url)
        image = PIL.Image.open(BytesIO(response.content))
        pdf_file_name = f'{id}.pdf'
        pdf_canvas = canvas.Canvas(pdf_file_name, pagesize=image.size)
        temp_image_file = f'{id}.png'
        image.save(temp_image_file)
        pdf_canvas.drawImage(temp_image_file, 0, 0)
        pdf_canvas.save()
        os.remove(temp_image_file)
        
        handles = [Patch(color=color, label=f'{status_label}') for status_label, color in color_legend.items()]
        fig, ax = plt.subplots(figsize=(6, 4), dpi=300)
        legend = ax.legend(handles=handles, loc='center', bbox_to_anchor=(0, 0, 1, 1))
        ax.axis('off')
        fig.tight_layout()
        plt.savefig('legend.pdf', bbox_inches='tight', pad_inches=0)
        plt.close()

        pdf1 = open(pathway_id + ".pdf", 'rb')
        pdf_reader1 = PyPDF2.PdfReader(pdf1)
        pdf2 = open(id + ".pdf", 'rb')
        pdf_reader2 = PyPDF2.PdfReader(pdf2)
        pdf3 = open('legend.pdf', 'rb')
        pdf_reader3 = PyPDF2.PdfReader(pdf3)
        pdf_writer = PyPDF2.PdfWriter()

        for page_num in range(len(pdf_reader1.pages)):
            page = pdf_reader1.pages[page_num]
            pdf_writer.add_page(page)
        for page_num in range(len(pdf_reader2.pages)):
            page = pdf_reader2.pages[page_num]
            pdf_writer.add_page(page)
        for page_num in range(len(pdf_reader3.pages)):
            page = pdf_reader3.pages[page_num]
            pdf_writer.add_page(page)

        output_file = open(id + "_" + output_name + '.pdf', 'wb')
        pdf_writer.write(output_file)

        pdf1.close()
        pdf2.close()
        pdf3.close()
        output_file.close()
        os.remove(pathway_id + ".pdf")
        os.remove(id + ".pdf")
        os.remove('legend.pdf')

    else:
        raise ValueError("Color legend cannot be used with cmap, vmin, vmax")

def filter_kegg_pathways_genes(filepath, sheet_name_paths, sheet_name_genes, genes_column, log2fc_column, count_threshold, benjamini_threshold , number_interventions = 1 , name_interventions = None):
    """
    Filter KEGG pathways based on specified criteria.

    Parameters:
    - filepath (str): The file path to the Excel file containing pathway and gene or transcript information.
    - sheet_name_paths (str): The sheet name containing pathway information.
    - sheet_name_genes (str): The sheet name containing gene/transcript information.
    - genes_column (str): The column name containing gene symbols in the sheet_name_genes sheet.
    - log2fc_column (str): The column name containing log2 fold change values in the sheet_name_genes sheet.
    - count_threshold (int): The count threshold for pathway inclusion.
    - benjamini_threshold (float or None): The Benjamini threshold for pathway inclusion.
    - number_interventions (int, optional): The number of interventions.
    - name_interventions (str or None, optional): The name of the interventions.

    Returns:
    dict: A dictionary containing filtered KEGG pathways and associated gene information.

    Example:
    >>> filter_kegg_pathways_genes('data.xlsx', 'Pathways', 'Genes', count_threshold=10, benjamini_threshold=0.05, number_interventions=2, name_interventions='Treatment A')
    {'pathway_id1': {'name': 'Pathway1', 'count': 15, 'genes': ['gene1', 'gene2'], 'benjamini': 0.03, 'logFC_dict': {'gene1': 2.5, 'gene2': -1.8}, 'logFC_secondary_dict': {'gene1': [2.5, 1.7, -0.8], 'gene2': [-1.8, 0.5, -2.0]}, 'intervention_number': 2, 'intervention_name': 'Treatment A'},
     'pathway_id2': {'name': 'Pathway2', 'count': 12, 'genes': ['gene3', 'gene4'], 'benjamini': 0.02, 'logFC_dict': {'gene3': 1.2, 'gene4': -3.0}, 'logFC_secondary_dict': {'gene3': [1.2, 0.8, -1.5], 'gene4': [-3.0, -2.5, -3.5]}, 'intervention_number': 2, 'intervention_name': 'Treatment A'}}
    """
    df = pd.read_excel(filepath, sheet_name=sheet_name_paths)
    gene_input = pd.read_excel(filepath, sheet_name=sheet_name_genes)

    if genes_column not in gene_input.columns or log2fc_column not in gene_input.columns:
        raise ValueError(f"Columns '{genes_column}' and '{log2fc_column}' must be present in the sheet '{sheet_name_paths}'")

    kegg_pathways = df[df['Category'] == 'KEGG_PATHWAY']
    results_dict = {}

    for _, row in kegg_pathways.iterrows():
        pathway_id = row['Term'].split(':')[0]
        pathway_name = row['Term'].split(':')[1]
        pathway_count = row['Count']
        pathway_genes = row['Genes'].split(', ')
        pathway_benjamini = row['Benjamini']
        
        if pathway_count >= count_threshold and (benjamini_threshold is None or pathway_benjamini <= benjamini_threshold):
            gene_logFC_dict = {}
            gene_logFC_secondary_dict = {}


            for gene in pathway_genes:               
                logFC_values = list(gene_input.loc[gene_input[genes_column] == gene, log2fc_column])
                max_logFC = max(logFC_values, key=abs)
                gene_logFC_dict[gene] = max_logFC
                gene_logFC_secondary_dict[gene] = logFC_values

            results_dict[pathway_id] = {'name': pathway_name,
                                        'count': pathway_count,
                                        'genes': pathway_genes,
                                        'benjamini': pathway_benjamini,
                                        'logFC_dict': gene_logFC_dict,
                                        'logFC_secondary_dict': gene_logFC_secondary_dict,
                                        'intervention_number': number_interventions,
                                        'intervention_name' : name_interventions}             

    return results_dict

def parse_bulk_kegg_pathway_file(filepath, sheet_name_paths, sheet_name_genes, genes_column, log2fc_column, number_interventions = 1 , name_interventions = None ):
    """
    Parse bulk KEGG pathway file and extract relevant information.
    No David format is required for this analysis. It features a simplified version which will map all the genes provided, on the requested pathways.

    Parameters:
    - filepath (str): The file path to the Excel file containing pathway and gene information.
    - sheet_name_paths (str): The sheet name containing pathway information.
    - sheet_name_genes (str): The sheet name containing gene information.
    - genes_column (str): The column name containing gene symbols.
    - log2fc_column (str): The column name containing log2 fold change values.
    - number_interventions (int, optional): The number of interventions.
    - name_interventions (str or None, optional): The name of the interventions.

    Returns:
    dict: A dictionary containing parsed KEGG pathway information.

    Example:
    >>> parse_bulk_kegg_pathway_file('data.xlsx', 'Pathways', 'Genes', 'GeneSymbol', 'Log2FC', number_interventions=2, name_interventions='Treatment A')
    {'pathway_id1': {'name': 'Pathway1', 'genes': ['gene1', 'gene2'], 'logFC_dict': {'gene1': 2.5, 'gene2': -1.8}, 'intervention_number': 2, 'intervention_name': 'Treatment A'},
     'pathway_id2': {'name': 'Pathway2', 'genes': ['gene3', 'gene4'], 'logFC_dict': {'gene3': 1.2, 'gene4': -3.0}, 'intervention_number': 2, 'intervention_name': 'Treatment A'}}
    """
    kegg_pathways = pd.read_excel(filepath, sheet_name=sheet_name_paths)
    gene_input = pd.read_excel(filepath, sheet_name=sheet_name_genes)

    if genes_column not in gene_input.columns or log2fc_column not in gene_input.columns:
        raise ValueError(f"Columns '{genes_column}' and '{log2fc_column}' must be present in the sheet '{sheet_name_paths}'")

    gene_input = gene_input[[genes_column , log2fc_column]]
    pathway_genes = gene_input[genes_column].to_list()
    results_dict = {}

    for _, row in kegg_pathways.iterrows():
        try:
            pathway_info = row['Term'].split(':')
            pathway_id = pathway_info[0]
            pathway_name = pathway_info[1] if len(pathway_info) > 1 else None
        except ValueError as error:
            print("An exception occurred:", type(error).__name__)
    
        gene_logFC_dict = {}

        for gene in pathway_genes:
            logFC_values = list(gene_input.loc[gene_input[genes_column] == gene, log2fc_column])
            max_logFC = max(logFC_values, key=abs)
            gene_logFC_dict[gene] = max_logFC

            results_dict[pathway_id] = {'name': pathway_name,
                                        'genes': pathway_genes,
                                        'logFC_dict': gene_logFC_dict,
                                        'intervention_number': number_interventions,
                                        'intervention_name' : name_interventions}

    return results_dict

def generate_colorscale_map(log2fc):
    """
    Generate a colorscale map based on log2 fold change values.

    Parameters:
    - log2fc (list): A list of log2 fold change values.

    Returns:
    tuple: A tuple containing the colormap, vmin, and vmax.

    Example:
    >>> generate_colorscale_map(log2fc=log2fc_values)
    """    
    vmin = min(log2fc)
    vmax = max(log2fc)

    num_colors = 10000

    if vmin >= 0 and vmax >= 0:
        vmin = 0
        cmap = LinearSegmentedColormap.from_list("", ["#FFFFFF", "#EE0000"])
    elif vmin < 0 and vmax < 0:
        vmax = 0
        cmap = LinearSegmentedColormap.from_list("", ["#008000", "#FFFFFF"])
    else:
        max_abs_value = max(map(abs, log2fc))
        vmin = -max_abs_value
        vmax = +max_abs_value


    greens_cmap = LinearSegmentedColormap.from_list("", ["#00AA00", "#FFFFFF"])
    reds_cmap = LinearSegmentedColormap.from_list("", ["#FFFFFF", "#EE0000"])
    colors = plt.cm.seismic(np.linspace(0, 1, num_colors))
    midpoint = int(num_colors * abs(vmin) / (abs(vmax) + abs(vmin)))
    greens = greens_cmap(np.linspace(0, 1, midpoint)).reshape(-1, 4)
    colors[:midpoint, :] = greens[:, :4]
    reds = reds_cmap(np.linspace(0, 1, num_colors - midpoint)).reshape(-1, 4)
    colors[midpoint:, :] = reds[:, :4]
    cmap = LinearSegmentedColormap.from_list("", colors)

    return cmap , vmin, vmax

def generate_genes_per_cell_spreadsheet(writer , genes_per_cell , id):
    """
    Generate a spreadsheet with genes per cell information.

    Parameters:
    - writer: Excel writer object.
    - genes_per_cell (dict): Dictionary containing genes per cell information.
    - id (str): Identifier for the sheet in the Excel file.

    Returns:
    None

    Example:
    >>> generate_genes_per_cell_spreadsheet(writer=writer , genes_per_cell=genes_per_cell , id=id)
    """
    genes_per_cell_to_df = {}

    for key, value in genes_per_cell.items():
        if isinstance(value, list):
            new_key = tuple(sorted(value))
        else:
            new_key = (value,)

        if new_key not in genes_per_cell_to_df:
            genes_per_cell_to_df[new_key] = []

        genes_per_cell_to_df[new_key].append(key)


    genes_per_cell_to_df = {tuple(sorted(k)): v for k, v in genes_per_cell_to_df.items()}
    df = pd.DataFrame([(v, k) for k, v in genes_per_cell_to_df.items()], columns=['Genes in Dataset','All Genes in Pathway cell' ])

    df['All Genes in Pathway cell'] = df['All Genes in Pathway cell'].apply(lambda x: ', '.join(x))
    df['Genes in Dataset'] = df['Genes in Dataset'].map(lambda x: ', '.join(map(str, x)))
    df.to_excel(writer, sheet_name=id, index=False, header=True)

def load_metadata(filepath):
    """
    Load metadata from a file specified by the given filepath.

    Parameters:
    - filepath (str): The path to the metadata file.

    Returns:
    pandas.DataFrame: A DataFrame containing the loaded metadata.

    Example:
    >>> data = load_metadata('metadata.csv')
    Opening CSV file: metadata.csv

    >>> data = load_metadata('metadata.tsv')
    Opening TSV file: metadata.tsv

    >>> data = load_metadata('metadata.xlsx')
    Opening Excel file: metadata.xlsx

    >>> data = load_metadata('unsupported_file.txt')
    Unsupported file format: unsupported_file.txt
    """
    file_extension = pathlib.Path(filepath).suffix

    if file_extension in csv_suffixes:
        print('Opening CSV file:', filepath)
        metadata = pd.read_csv(filepath, sep=',')

    elif file_extension in tsv_suffixes:
        print('Opening TSV file:', filepath)
        metadata = pd.read_csv(filepath, sep='\t')

    elif file_extension in excel_suffixes:
        print('Opening Excel file:', filepath)
        try:
            metadata = pd.read_excel(filepath)
        except pd.errors.ParserError:
            print(f"Error reading Excel file: {filepath}")
    else:
        print('Unsupported file format:', filepath)

    return metadata

def evaluate_metadata(metadata_df , pval_column , genes_column):
    """
    Evaluate metadata DataFrame for proper structure and provided columns.

    Parameters:
    - metadata_df (pd.DataFrame): The metadata DataFrame.
    - pval_column (str): The name of the p-value column.
    - genes_column (str): The name of the genes column.

    Raises:
    - ValueError: If metadata_df is None or not a DataFrame.
    - ValueError: If pval_column or genes_column is not provided, not a string, or not found in metadata_df columns.

    Example:
    >>> evaluate_metadata(metadata_df, 'PValue', 'Gene')
    """
    if metadata_df is None:
        raise ValueError(f'Could not load {metadata_df} properly. perhaps check the metadata file for correct structure')
    elif not isinstance(metadata_df, pd.DataFrame):
        raise ValueError('metadata df is not of the correct type (pd.DataFrame)')

    if pval_column is None:
        print('pvalue column is not provided')
    elif not isinstance(pval_column, str):
        print('pvalue column is not of the correct type (str)')
    elif pval_column not in metadata_df.columns:
        raise ValueError(f'{pval_column} not found in metadata dataframe columns')

    if genes_column is None:
        raise ValueError('genes column is not provided')
    elif not isinstance(genes_column, str):
        raise ValueError('genes column is not of the correct type (str)')
    elif genes_column not in metadata_df.columns:
        raise ValueError(f'{genes_column} not found in metadata dataframe columns')











