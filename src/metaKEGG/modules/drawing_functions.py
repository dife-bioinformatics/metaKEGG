from Bio.KEGG import REST
from Bio.Graphics.KGML_vis import KGMLCanvas
import pandas as pd
from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
from pylab import *
from itertools import combinations

from ..config import gray as gray
from ..helpers import helpfunctions as _hf

def make_new_graphic(input_entry):
    """
    Create a new graphic entry based on the provided input entry.

    Parameters:
    - input_entry (KGML_parser.Entry): Input KGML entry containing graphics information.

    Returns:
    KGML_parser.Entry: New KGML entry with graphics information.
    """
    output_entry = KGML_parser.Graphics(input_entry)
    output_entry.graphics = KGML_parser.Entry.add_graphics(input_entry , output_entry)
    output_entry.graphics = input_entry.graphics
    output_entry.name = input_entry.graphics[0].name
    output_entry.x = input_entry.graphics[0].x
    output_entry.y = input_entry.graphics[0].y
    output_entry.type = input_entry.graphics[0].type
    output_entry.width = input_entry.graphics[0].width
    output_entry.height = input_entry.graphics[0].height
    output_entry.fgcolor   = input_entry.graphics[0].fgcolor
    output_entry.bgcolor  = input_entry.graphics[0].bgcolor
        
    return output_entry

def draw_KEGG_pathways_genes(parsed_output , info , save_to_eps):
    """
    Draw KEGG pathways and save to file, with gene expression information.

    Parameters:
    - parsed_output (dict): Parsed output containing log2 fold change information for genes.
    - info (dict): Information about pathways and genes.
    - save_to_eps (bool): Flag to save output to EPS format.

    Returns:
    None

    Example:
    >>> draw_KEGG_pathways_genes(parsed_output, info, save_to_eps=True)
    """    
    writer = pd.ExcelWriter('genes_per_cell.xlsx', engine='xlsxwriter')

    for id , path_data in parsed_output.items():
        genes_per_cell = {}
        pathway_id = info[id]['corresponding_KO']
        pathway = KGML_parser.read(REST.kegg_get(pathway_id, "kgml"))
        log2fc = path_data['logFC_dict']
        gene_symbol_KO = info[id]['gene_symbol_KO']
        output_name = _hf.file_naming_scheme(input_data=parsed_output , id = id)

        log2fc_values = log2fc.values()

        cmap , vmin, vmax = _hf.generate_colorscale_map(log2fc=log2fc_values)

        for entry in pathway.orthologs:
            entry.graphics[0].name = ""
            multiple_kos = [name.split(":")[1] for name in entry.name.split()]
            corresponding_genes = [key for key, values in gene_symbol_KO.items() if values in multiple_kos]

            if not corresponding_genes:
                entry.graphics[0].name = ""
                entry.graphics[0].bgcolor = gray
            if entry.graphics[0].type == 'line':
                continue
            if  len(corresponding_genes) > 1:
                for ko in range(len(corresponding_genes) -1):
                    new_entry = make_new_graphic(entry)

                num_subcells = len(entry.graphics)
                subcell_width = new_entry.width / num_subcells
                left_x = new_entry.graphics[0].x - new_entry.width/2
            
            
                for i, (subcell, gene) in enumerate(zip(entry.graphics, corresponding_genes)):
                    subcell.x = left_x + subcell_width * (i + 0.5)
                    subcell.width = subcell_width
                    subcell.bgcolor = gray
                    subcell.name = ""

                    if gene in log2fc:
                        subcell.bgcolor = cmap((log2fc[gene] - vmin) / (vmax - vmin))
                        subcell.name = gene
                        genes_per_cell[gene] = corresponding_genes

            else:
                for i, (element, gene) in enumerate(zip(entry.graphics, corresponding_genes)):
                    element.bgcolor = gray
                    element.name = gene
                    if gene in log2fc:
                        element.bgcolor = cmap((log2fc[gene] - vmin) / (vmax - vmin))
                        genes_per_cell[gene] = gene

        _hf.generate_genes_per_cell_spreadsheet(writer=writer , genes_per_cell=genes_per_cell , id=id)

        cnvs = KGMLCanvas(pathway, import_imagemap=True , fontsize=9)
        cnvs.draw(pathway_id + ".pdf")

        if save_to_eps:
            cnvs.draw(id + "_" + output_name + ".eps")
        
        _hf.compile_and_write_output_files(id=id, pathway_id=pathway_id , cmap=cmap , vmin=vmin , vmax=vmax , output_name=output_name , save_to_eps=save_to_eps)

    writer.close()

def draw_KEGG_pathways_transcripts(parsed_output , info, save_to_eps):
    """
    Draw KEGG pathways and save to file, with transcript expression information.

    Parameters:
    - parsed_output (dict): Parsed output containing log2 fold change information for genes and transcripts.
    - info (dict): Information about pathways, genes, and transcripts.
    - save_to_eps (bool): Flag to save output to EPS format.

    Returns:
    None

    Example:
    >>> draw_KEGG_pathways_transcripts(parsed_output, info, save_to_eps=True)
    """
    writer = pd.ExcelWriter('genes_per_cell.xlsx', engine='xlsxwriter')

    for id , path_data in parsed_output.items():
        genes_per_cell = {}
        transcripts_per_cell = {}
        pathway_id = info[id]['corresponding_KO']
        pathway = KGML_parser.read(REST.kegg_get(pathway_id, "kgml"))
        log2fc = path_data['logFC_dict']
        log2fc_secondary = path_data['logFC_secondary_dict']
        gene_symbol_KO = info[id]['gene_symbol_KO']
        output_name = _hf.file_naming_scheme(input_data=parsed_output , id = id)
        corresponding_transcripts = [key for key, values in log2fc_secondary.items() if len(values) > 1]

        log2fc_extended = list(log2fc.values()) + [value for values in log2fc_secondary.values() for value in values]

        cmap , vmin, vmax = _hf.generate_colorscale_map(log2fc=log2fc_extended)

        for entry in pathway.orthologs:
            entry.graphics[0].name = ""
            multiple_kos = [name.split(":")[1] for name in entry.name.split()]
            corresponding_genes = [key for key, values in gene_symbol_KO.items() if values in multiple_kos]
            if not corresponding_genes:
                entry.graphics[0].name = ""
                entry.graphics[0].bgcolor = gray
            if entry.graphics[0].type == 'line':
                continue
            if entry.graphics[0].type != 'rectangle':
                continue
            if  len(corresponding_genes) > 1:
                for ko in range(len(corresponding_genes) -1):
                    new_entry = make_new_graphic(entry)

                num_subcells = len(entry.graphics)
                subcell_width = new_entry.width / num_subcells
                left_x = new_entry.graphics[0].x - new_entry.width/2
            
            
                for i, (subcell, gene) in enumerate(zip(entry.graphics, corresponding_genes)):
                    subcell.x = left_x + subcell_width * (i + 0.5)
                    subcell.width = subcell_width
                    subcell.bgcolor =  gray
                    subcell.name = ""
                    if gene in log2fc:
                        subcell.bgcolor = cmap((log2fc[gene] - vmin) / (vmax - vmin))
                        subcell.name = gene
                        genes_per_cell[gene] = corresponding_genes
                    if gene in corresponding_transcripts:
                        for transcript in range(len(log2fc_secondary[gene]) -1):
                            try:
                                new_entry_v = make_new_graphic(subcell)
                            except: pass
                        num_subcells_vertical = len(log2fc_secondary[gene])
                        subcell_height = new_entry_v.height / num_subcells_vertical
                        top_y = new_entry_v.graphics[0].y - new_entry_v.height / 2
                        for j in range(num_subcells_vertical):
                            new_entry_v.y = top_y + subcell_height * (j + 0.5)
                            new_entry_v.x = subcell.x
                            new_entry_v.height = subcell_height
                            new_entry_v.bgcolor = cmap((log2fc_secondary[gene][j] - vmin) / (vmax - vmin))
                            new_entry_v.name = ''
                            transcripts_per_cell[gene] = corresponding_transcripts
            else:
                for i, (element, gene) in enumerate(zip(entry.graphics, corresponding_genes)):
                    if entry.graphics[0].type == 'line':
                        continue
                    element.bgcolor = gray
                    if gene in log2fc:
                        element.bgcolor = cmap((log2fc[gene] - vmin) / (vmax - vmin))
                        element.name = gene
                    if gene in corresponding_transcripts:
                        for transcript in range(len(log2fc_secondary[gene]) -1):
                            try:
                                new_entry_v = make_new_graphic(entry)
                            except: pass
                        num_subcells_vertical = len(log2fc_secondary[gene])
                        subcell_height = new_entry_v.height / num_subcells_vertical
                        top_y = new_entry_v.graphics[0].y - new_entry_v.height / 2
                        for j in range(num_subcells_vertical):
                            new_entry_v.y = top_y + subcell_height * (j + 0.5)
                            new_entry_v.x = element.x
                            new_entry_v.height = subcell_height
                            new_entry_v.bgcolor = cmap((log2fc_secondary[gene][j] - vmin) / (vmax - vmin))
                            new_entry_v.name = ''
                            transcripts_per_cell[gene] = corresponding_transcripts

                            
        _hf.generate_genes_per_cell_spreadsheet(writer=writer , genes_per_cell=genes_per_cell , id=id)

        cnvs = KGMLCanvas(pathway, import_imagemap=True , fontsize=9)
        cnvs.draw(pathway_id + ".pdf")

        if save_to_eps:
            cnvs.draw(id + "_" + output_name + ".eps")
        
        _hf.compile_and_write_output_files(id=id, pathway_id=pathway_id , cmap=cmap , vmin=vmin , vmax=vmax , output_name=output_name , save_to_eps=save_to_eps)

    
    writer.close()

def draw_KEGG_pathways_genes_multiple_interventions(parsed_out_list , intervention_names , colors_list , save_to_eps):
    """
    Draw KEGG pathways with gene expression information for multiple interventions.

    Parameters:
    - parsed_out_list (list of dicts): List of parsed output dictionaries containing log2 fold change information for genes.
    - intervention_names (list of str): List of intervention names.
    - colors_list (list of str): List of colors corresponding to each intervention.
    - save_to_eps (bool): Flag to save output to EPS format.

    Returns:
    None

    Example:
    >>> draw_KEGG_pathways_genes_multiple_interventions(parsed_out_list, intervention_names, colors_list, save_to_eps=True)
    """
    num_interventions = len(intervention_names)
    color_to_intervention = { inter : color for (inter , color) in zip(intervention_names , colors_list)}

    key_sets = [set(d.keys()) for d in parsed_out_list]
    common_keys = list(set.intersection(*key_sets))
    if not common_keys:
        raise ValueError("No common keys found in the list of dictionaries.")

    common_key_dict = {key: [] for key in common_keys}

    for d in parsed_out_list:
        for key in common_keys:
            if key in d:
                common_key_dict[key].append({key: d[key]})


    combination_dict = {common_key: [] for common_key in common_key_dict}
    for common_key , common_list in common_key_dict.items():
        for r in range(2, num_interventions + 1):
            for combination in combinations(common_key_dict[common_key], r):
                combination_key = '+'.join(sorted(interv_dict[common_key]['intervention_name'] for interv_dict in combination))
                genes_lists = [interv_dict[common_key]['genes'] for interv_dict in combination]
                common_genes = list(set.intersection(*map(set, genes_lists)))
                level = len(combination)
                combination_dict[common_key].append({'level': level, 'combination_key': combination_key, 'common_genes': common_genes})

        max_level = max(combination_info['level'] for combinations_list in combination_dict.values() for combination_info in combinations_list)

        for common_key, combinations_list in combination_dict.items():
            for i, combination_info in enumerate(combinations_list):
                level_i_genes = combination_info['common_genes']

                for j in range(combination_info['level'] + 1, max_level + 1):
                    for k in range(i + 1, len(combinations_list)):
                        level_k_genes = combinations_list[k]['common_genes']
                        common_genes_to_remove = set(level_i_genes) & set(level_k_genes)
                        if common_genes_to_remove:
                            combinations_list[i]['common_genes'] = list(set(level_i_genes) - common_genes_to_remove)

        for common_key, combinations_list in combination_dict.items():
            for combination_info in combinations_list:
                level = combination_info['level']
                combination_key = combination_info['combination_key']
                common_genes = combination_info['common_genes']
                
                
        intervention_names_w_combos = intervention_names.copy()

        for common_key, combinations_list in combination_dict.items():
            for combination_info in combinations_list:
                intervention_names_w_combos.append(combination_info['combination_key'])

        color_to_intervention = {inter: color for (inter, color) in zip(intervention_names_w_combos, colors_list)}

    writer = pd.ExcelWriter('genes_per_cell.xlsx', engine='xlsxwriter')

    for common_key , common_list in common_key_dict.items():
        collect_info = _hf.collect_pathway_info_multiple_interventions(common_key)
        genes_per_cell = {}
        pathway_id = collect_info[common_key]['corresponding_KO']
        pathway = KGML_parser.read(REST.kegg_get(pathway_id, "kgml"))
        gene_symbol_KO = collect_info[common_key]['gene_symbol_KO']
        output_name = _hf.file_naming_scheme(input_data=common_list[0][common_key]['name'])
        
        for entry in pathway.orthologs:
            entry.graphics[0].name = ''
            multiple_kos = [name.split(":")[1] for name in entry.name.split()]
            corresponding_genes = [key for key, values in gene_symbol_KO.items() if values in multiple_kos]

            if not corresponding_genes:
                entry.graphics[0].name = ""
                entry.graphics[0].bgcolor = gray
            if entry.graphics[0].type == 'line':
                continue
            if  len(corresponding_genes) > 1:
                for ko in range(len(corresponding_genes) -1):
                    new_entry = make_new_graphic(entry)

                num_subcells = len(entry.graphics)
                subcell_width = new_entry.width / num_subcells
                left_x = new_entry.graphics[0].x - new_entry.width/2
            
            
                for i, (subcell, gene) in enumerate(zip(entry.graphics, corresponding_genes)):
                    subcell.x = left_x + subcell_width * (i + 0.5)
                    subcell.width = subcell_width
                    subcell.bgcolor = gray
                    subcell.name = ""

                    for indiv_pathway in common_list:
                        genes_in_pathway = indiv_pathway[common_key]['genes']
                        name_intervention = indiv_pathway[common_key]['intervention_name']
                        in_pathway = gene in genes_in_pathway
                        key_with_gene = [combo['combination_key'] for combo in combination_dict[common_key] if gene in combo['common_genes']]
                        if len(key_with_gene) >= 1:
                                subcell.bgcolor = color_to_intervention[key_with_gene[0]]
                                subcell.name = gene
                                genes_per_cell[gene] = corresponding_genes
                        elif in_pathway:
                            subcell.bgcolor = color_to_intervention[name_intervention]
                            subcell.name = gene
                            genes_per_cell[gene] = corresponding_genes



            else:
                for i, (element, gene) in enumerate(zip(entry.graphics, corresponding_genes)):
                    element.bgcolor = gray
                    element.name = gene

                    for indiv_pathway in common_list:
                        genes_in_pathway = indiv_pathway[common_key]['genes']
                        name_intervention = indiv_pathway[common_key]['intervention_name']
                        in_pathway = gene in genes_in_pathway
                        key_with_gene = [combo['combination_key'] for combo in combination_dict[common_key] if gene in combo['common_genes']]

                        if len(key_with_gene) >= 1:
                            element.bgcolor = color_to_intervention[key_with_gene[0]]
                            element.name = gene
                            genes_per_cell[gene] = corresponding_genes
                        elif in_pathway:
                            element.bgcolor = color_to_intervention[name_intervention]
                            element.name = gene
                            genes_per_cell[gene] = corresponding_genes

        _hf.generate_genes_per_cell_spreadsheet(writer=writer , genes_per_cell=genes_per_cell , id=common_key)

        cnvs = KGMLCanvas(pathway, import_imagemap=True , fontsize=9)
        cnvs.draw(pathway_id + ".pdf")

        if save_to_eps:
            cnvs.draw(id + "_" + output_name + ".eps")

        _hf.compile_and_write_output_files(id=common_key, pathway_id=pathway_id , color_legend=color_to_intervention , output_name=output_name , save_to_eps=save_to_eps)

    writer.close()

def draw_KEGG_pathways_genes_with_methylation(parsed_output , info , genes_from_MM , color_legend, save_to_eps):
    """
    Draw KEGG pathways with gene expression and methylation information.

    Parameters:
    - parsed_output (dict): Parsed output containing gene expression information.
    - info (dict): Information about pathways, genes, and methylation status.
    - genes_from_MM (list): List of genes that are differentially methylated.
    - color_legend (dict): Dictionary mapping methylation status labels to colors.
    - save_to_eps (bool): Flag to save output to EPS format.

    Returns:
    None

    Example:
    >>> draw_KEGG_pathways_genes_with_methylation(parsed_output, info, genes_from_MM, color_legend, save_to_eps=True)
    """
    writer = pd.ExcelWriter('genes_per_cell.xlsx', engine='xlsxwriter')
    for id , path_data in parsed_output.items():
        genes_per_cell = {}
        pathway_id = info[id]['corresponding_KO']
        pathway = KGML_parser.read(REST.kegg_get(pathway_id, "kgml"))
        genes_in_pathway = path_data['genes']
        gene_symbol_KO = info[id]['gene_symbol_KO']
        output_name = _hf.file_naming_scheme(input_data=parsed_output , id = id)


        for entry in pathway.orthologs:
            entry.graphics[0].bgcolor = gray
            entry.graphics[0].name = ''
            multiple_kos = [name.split(":")[1] for name in entry.name.split()]
            corresponding_genes = [key for key, values in gene_symbol_KO.items() if values in multiple_kos]

            if not corresponding_genes:
                entry.graphics[0].name = ""
                entry.graphics[0].bgcolor = gray
            if entry.graphics[0].type == 'line':
                continue
            if  len(corresponding_genes) > 1:
                for ko in range(len(corresponding_genes) -1):
                    new_entry = make_new_graphic(entry)

                num_subcells = len(entry.graphics)
                subcell_width = new_entry.width / num_subcells
                left_x = new_entry.graphics[0].x - new_entry.width/2
            
            
                for i, (subcell, gene) in enumerate(zip(entry.graphics, corresponding_genes)):
                    subcell.x = left_x + subcell_width * (i + 0.5)
                    subcell.width = subcell_width
                    subcell.bgcolor = gray
                    subcell.name = ""

                    if (gene in genes_in_pathway) and (gene in genes_from_MM):
                        subcell.bgcolor = color_legend['Differentially methylated']
                        subcell.name = gene
                        genes_per_cell[gene] = corresponding_genes
                    elif (gene in genes_in_pathway) and not (gene in genes_from_MM):
                        subcell.bgcolor = color_legend['Not differentially methylated']
                        subcell.name = gene

            else:
                for i, (element, gene) in enumerate(zip(entry.graphics, corresponding_genes)):
                    element.bgcolor = gray
                    element.name = gene
                    if (gene in genes_in_pathway) and (gene in genes_from_MM):
                        element.bgcolor = color_legend['Differentially methylated']
                        genes_per_cell[gene] = gene
                    elif (gene in genes_in_pathway) and not (gene in genes_from_MM):
                        element.bgcolor = color_legend['Not differentially methylated']

        _hf.generate_genes_per_cell_spreadsheet(writer=writer , genes_per_cell=genes_per_cell , id=id)

        cnvs = KGMLCanvas(pathway, import_imagemap=True , fontsize=9)
        cnvs.draw(pathway_id + ".pdf")

        if save_to_eps:
            cnvs.draw(id + "_" + output_name + ".eps")

        _hf.compile_and_write_output_files(id=id, pathway_id=pathway_id , color_legend=color_legend , output_name=output_name , save_to_eps=save_to_eps)

    writer.close()

def draw_KEGG_pathways_genes_with_miRNA(parsed_output , info , genes_from_miRNA , color_legend, save_to_eps):
    """
    Draw KEGG pathways with gene expression and miRNA information.

    Parameters:
    - parsed_output (dict): Parsed output containing gene expression information.
    - info (dict): Information about pathways, genes, and miRNA status.
    - genes_from_miRNA (list): List of genes that are detected by miRNA.
    - color_legend (dict): Dictionary mapping miRNA status labels to colors.
    - save_to_eps (bool): Flag to save output to EPS format.

    Returns:
    None

    Example:
    >>> draw_KEGG_pathways_genes_with_miRNA(parsed_output, info, genes_from_miRNA, color_legend, save_to_eps=True)
    """
    writer = pd.ExcelWriter('genes_per_cell.xlsx', engine='xlsxwriter')
    for id , path_data in parsed_output.items():
        genes_per_cell = {}
        pathway_id = info[id]['corresponding_KO']
        pathway = KGML_parser.read(REST.kegg_get(pathway_id, "kgml"))
        genes_in_pathway = path_data['genes']
        gene_symbol_KO = info[id]['gene_symbol_KO']
        output_name = _hf.file_naming_scheme(input_data=parsed_output , id = id)


        for entry in pathway.orthologs:
            entry.graphics[0].bgcolor = gray
            entry.graphics[0].name = ''
            multiple_kos = [name.split(":")[1] for name in entry.name.split()]
            corresponding_genes = [key for key, values in gene_symbol_KO.items() if values in multiple_kos]

            if not corresponding_genes:
                entry.graphics[0].name = ""
                entry.graphics[0].bgcolor = gray
            if entry.graphics[0].type == 'line':
                continue
            if  len(corresponding_genes) > 1:
                for ko in range(len(corresponding_genes) -1):
                    new_entry = make_new_graphic(entry)

                num_subcells = len(entry.graphics)
                subcell_width = new_entry.width / num_subcells
                left_x = new_entry.graphics[0].x - new_entry.width/2
            
            
                for i, (subcell, gene) in enumerate(zip(entry.graphics, corresponding_genes)):
                    subcell.x = left_x + subcell_width * (i + 0.5)
                    subcell.width = subcell_width
                    subcell.bgcolor = gray
                    subcell.name = ""

                    if (gene in genes_in_pathway) and (gene in genes_from_miRNA):
                        subcell.bgcolor = color_legend['miRNA detected']
                        subcell.name = gene
                        genes_per_cell[gene] = corresponding_genes
                    elif (gene in genes_in_pathway) and not (gene in genes_from_miRNA):
                        subcell.bgcolor = color_legend['miRNA not detected']
                        subcell.name = gene

            else:
                for i, (element, gene) in enumerate(zip(entry.graphics, corresponding_genes)):
                    element.bgcolor = gray
                    element.name = gene
                    if (gene in genes_in_pathway) and (gene in genes_from_miRNA):
                        element.bgcolor = color_legend['miRNA detected']
                        genes_per_cell[gene] = gene
                    elif (gene in genes_in_pathway) and not (gene in genes_from_miRNA):
                        element.bgcolor = color_legend['miRNA not detected']

        _hf.generate_genes_per_cell_spreadsheet(writer=writer , genes_per_cell=genes_per_cell , id=id)

        cnvs = KGMLCanvas(pathway, import_imagemap=True , fontsize=9)
        cnvs.draw(pathway_id + ".pdf")

        if save_to_eps:
            cnvs.draw(id + "_" + output_name + ".eps")

        _hf.compile_and_write_output_files(id=id, pathway_id=pathway_id , color_legend=color_legend , output_name=output_name , save_to_eps=save_to_eps)

    writer.close()

def draw_KEGG_pathways_genes_with_methylation_and_miRNA(parsed_output , info , genes_from_MM , genes_from_miRNA, color_legend, save_to_eps):
    """
    Draw KEGG pathways with gene expression, methylation, and miRNA information.

    Parameters:
    - parsed_output (dict): Parsed output containing gene expression information.
    - info (dict): Information about pathways, genes, and miRNA status.
    - genes_from_MM (list): List of genes that are differentially methylated.
    - genes_from_miRNA (list): List of genes that are detected by miRNA.
    - color_legend (dict): Dictionary mapping gene status labels to colors.
    - save_to_eps (bool): Flag to save output to EPS format.

    Returns:
    None

    Example:
    >>> draw_KEGG_pathways_genes_with_methylation_and_miRNA(parsed_output, info, genes_from_MM, genes_from_miRNA, color_legend, save_to_eps=True)
    """
    writer = pd.ExcelWriter('genes_per_cell.xlsx', engine='xlsxwriter')
    for id , path_data in parsed_output.items():
        genes_per_cell = {}
        pathway_id = info[id]['corresponding_KO']
        pathway = KGML_parser.read(REST.kegg_get(pathway_id, "kgml"))
        genes_in_pathway = path_data['genes']
        gene_symbol_KO = info[id]['gene_symbol_KO']
        output_name = _hf.file_naming_scheme(input_data=parsed_output , id = id)


        for entry in pathway.orthologs:
            entry.graphics[0].bgcolor = gray
            entry.graphics[0].name = ''
            multiple_kos = [name.split(":")[1] for name in entry.name.split()]
            corresponding_genes = [key for key, values in gene_symbol_KO.items() if values in multiple_kos]

            if not corresponding_genes:
                entry.graphics[0].name = ""
                entry.graphics[0].bgcolor = gray
            if entry.graphics[0].type == 'line':
                continue
            if  len(corresponding_genes) > 1:
                for ko in range(len(corresponding_genes) -1):
                    new_entry = make_new_graphic(entry)

                num_subcells = len(entry.graphics)
                subcell_width = new_entry.width / num_subcells
                left_x = new_entry.graphics[0].x - new_entry.width/2
            
            
                for i, (subcell, gene) in enumerate(zip(entry.graphics, corresponding_genes)):
                    subcell.x = left_x + subcell_width * (i + 0.5)
                    subcell.width = subcell_width
                    subcell.bgcolor = gray
                    subcell.name = ""

                    if (gene in genes_in_pathway) and (gene in genes_from_MM) and (gene in genes_from_miRNA):
                        subcell.bgcolor = color_legend['Differentially methylated and miRNA detected']
                        subcell.name = gene
                        genes_per_cell[gene] = corresponding_genes
                    elif (gene in genes_in_pathway) and not (gene in genes_from_MM) and not (gene in genes_from_miRNA):
                        subcell.bgcolor = color_legend['Not differentially methylated and not miRNA detected']
                        subcell.name = gene
                    elif (gene in genes_in_pathway) and not (gene in genes_from_MM) and (gene in genes_from_miRNA):
                        subcell.bgcolor = color_legend['Not differentially methylated and miRNA detected']
                        subcell.name = gene
                    elif (gene in genes_in_pathway) and (gene in genes_from_MM) and not (gene in genes_from_miRNA):
                        subcell.bgcolor = color_legend['Differentially methylated and not miRNA detected']
                        subcell.name = gene

            else:
                for i, (element, gene) in enumerate(zip(entry.graphics, corresponding_genes)):
                    element.bgcolor = gray
                    element.name = gene
                    if (gene in genes_in_pathway) and (gene in genes_from_MM) and (gene in genes_from_miRNA):
                        element.bgcolor = color_legend['Differentially methylated and miRNA detected']
                        genes_per_cell[gene] = gene
                    elif (gene in genes_in_pathway) and not (gene in genes_from_MM) and not (gene in genes_from_miRNA):
                        element.bgcolor = color_legend['Not differentially methylated and not miRNA detected']
                    elif (gene in genes_in_pathway) and not (gene in genes_from_MM) and (gene in genes_from_miRNA):
                        element.bgcolor = color_legend['Not differentially methylated and miRNA detected']
                    elif (gene in genes_in_pathway) and (gene in genes_from_MM) and not (gene in genes_from_miRNA):
                        element.bgcolor = color_legend['Differentially methylated and not miRNA detected']


        _hf.generate_genes_per_cell_spreadsheet(writer=writer , genes_per_cell=genes_per_cell , id=id)

        cnvs = KGMLCanvas(pathway, import_imagemap=True , fontsize=9)
        cnvs.draw(pathway_id + ".pdf")

        if save_to_eps:
            cnvs.draw(id + "_" + output_name + ".eps")

        _hf.compile_and_write_output_files(id=id, pathway_id=pathway_id , color_legend=color_legend , output_name=output_name , save_to_eps=save_to_eps)

    writer.close()