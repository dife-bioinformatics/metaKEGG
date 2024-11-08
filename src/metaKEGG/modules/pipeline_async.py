import os
import datetime
import sys
import shutil
import asyncio
from typing import Optional, List, Union
from pathlib import Path
from ..helpers import helpfunctions as _hf
from ..modules import drawing_functions as _df
from ..modules import colorscale as _cs
from ..config import analysis_types_to_execute as analysis_types_to_execute

class PipelineAsync:
    """
    Class for executing different analyses on KEGG pathways with various data inputs.
    """
    def __init__(
        self,
        input_file_path: Union[str, Path, List[str], List[Path]],
        input_label: Optional[Union[str, List[str]]], 
        sheet_name_paths: str = "pathways",
        sheet_name_genes: str = "gene_metrics",
        genes_column: str = "gene_symbol",
        log2fc_column: str = "logFC",
        output_folder_name: Optional[str] = None,
        folder_extension: Optional[str] = None,
        compounds_list: Optional[List[str]] = None,
        save_to_eps: bool = False

    ) -> None:
        
        self.input_file_path = input_file_path
        self.input_label = input_label

        self.sheet_name_paths = sheet_name_paths
        self.sheet_name_genes = sheet_name_genes
        self.genes_column = genes_column
        self.log2fc_column = log2fc_column

        self.compounds_list = compounds_list if compounds_list is not None else []
        self.output_folder_name = output_folder_name
        self.folder_extension = folder_extension
        self.save_to_eps = save_to_eps


    def find_file_folder(self):
        """
        Find the folder containing the input file(s) and set the current working directory to that location.

        Returns:
        str: The path to the folder containing the input file(s).

        Raises:
        FileNotFoundError: If the specified file or the first file in the list (for multiple files) does not exist.
        """
        if isinstance(self.input_file_path, list):
            if os.path.exists(self.input_file_path[0]):
                folder_path = os.path.dirname(self.input_file_path[0])
                os.chdir(folder_path)
                return folder_path
            else:
                raise FileNotFoundError(f"The file '{self.input_file_path[0]}' does not exist.")
        
        else:
            if os.path.exists(self.input_file_path):
                folder_path = os.path.dirname(self.input_file_path)
                os.chdir(folder_path)
                return folder_path
            else:
                raise FileNotFoundError(f"The file '{self.input_file_path}' does not exist.")
        
    def make_output_folder(self , folder_path , analysis_extension):
        """
        Create a unique output folder based on the current date, analysis extension, and an optional folder extension.

        Args:
        folder_path (str): The path to the folder where the output folder will be created.
        analysis_extension (str): A string representing the type of analysis being performed.

        Returns:
        str: The path to the created output folder.

        Notes:
        - The output folder is named in the format: "draw_KEGG_<current_date>_<analysis_extension>_<folder_extension>"
        - The folder_extension is optional and can be set during the class initialization.
        """
        if self.output_folder_name and self.output_folder_name is not None:
            create_folder = os.path.join(folder_path, self.output_folder_name)
            print(f'Will create folder : {create_folder}')
            output_folder = _hf.create_output_folder(create_folder, self.folder_extension)
        else:
            today = datetime.date.today().strftime("%Y-%m-%d")
            folder_today = f"draw_KEGG_{today}_{analysis_extension}"
            create_folder = os.path.join(folder_path, folder_today)
            print(f'Will create folder : {create_folder}')
            output_folder = _hf.create_output_folder(create_folder, self.folder_extension)
        return output_folder

    async def single_input_genes(self, benjamini_threshold: Optional[float] = None, count_threshold: Optional[int] = 2 , pathway_pvalue: Optional[float] = None):
        """
        Perform the Single Input Analysis for Gene IDs.

        Raises:
        TypeError: If the input_file_path is a list, as this analysis expects a single input file.

        Prints:
        - Execution message.
        - Output folder path.
        - Parsing and collecting pathway information messages.
        - Completion message.

        Notes:
        - Calls helper functions to filter KEGG pathways for genes, parse the input file, and draw KEGG pathways.
        - The output files are located in the created output folder.
        """
        if isinstance(self.input_file_path , list):
            raise TypeError('Please provide a single input to perform \'Single input analysis (Genes)')

        print("Executing analysis: Single input (Gene IDs)...\n")
        entry_dir = os.getcwd()
        folder_of_input = self.find_file_folder()

        analysis_extension = 'genes'
        output_folder = self.make_output_folder(folder_of_input , analysis_extension)
        print(f'Output folder is {output_folder}\n')
        print('Parsing input file...\n')
        parsed_out, all_genes = _hf.filter_kegg_pathways_genes(filepath=self.input_file_path,
                                                    sheet_name_paths=self.sheet_name_paths,
                                                    sheet_name_genes=self.sheet_name_genes,
                                                    genes_column= self.genes_column,
                                                    log2fc_column=self.log2fc_column,
                                                    count_threshold = count_threshold,
                                                    benjamini_threshold=benjamini_threshold,
                                                    raw_pvalue_threshold=pathway_pvalue)
        if len(parsed_out) == 0:
            raise ValueError("Could not detect pathways in the input file with the selected default & user settings. Check your input file and/or settings.")
        print('Finished parsing input file\n')
        print('Collecting pathway info...')
        pathway_info = await _hf.collect_pathway_info_async(parsed_output=parsed_out)
        print('Finished collecting pathway info\n')
        os.chdir(output_folder)
        _hf.generate_pathways_per_gene_spreadsheet(gene_list=all_genes, pathway_dict=parsed_out , name_extension=None)
        print('Mapping pathways...\n')
        _df.draw_KEGG_pathways_genes(parsed_output=parsed_out , info=pathway_info , compounds_list=self.compounds_list ,save_to_eps=self.save_to_eps)
        os.chdir(entry_dir)
        print(f'\nDone! \nOutput files are located in {output_folder}')

    async def single_input_transcripts(self, benjamini_threshold: Optional[float] = None, count_threshold: Optional[int] = 2 , pathway_pvalue: Optional[float] = None):
        """
        Perform the Single Input Analysis for Transcript IDs.

        Raises:
        TypeError: If the input_file_path is a list, as this analysis expects a single input file.

        Prints:
        - Execution message.
        - Output folder path.
        - Parsing and collecting pathway information messages.
        - Completion message.

        Notes:
        - Calls helper functions to filter KEGG pathways for genes, parse the input file, and draw KEGG pathways for transcripts.
        - The output files are located in the created output folder.
        """    
        if isinstance(self.input_file_path , list):
            raise TypeError('Please provide a single input to perform \'Single input analysis (Transcripts)')
        
        print("Executing analysis: Single input (Transcript IDs)...\n")
        entry_dir = os.getcwd()
        folder_of_input = self.find_file_folder()
        analysis_extension = 'transcripts'
        output_folder = self.make_output_folder(folder_of_input , analysis_extension)
        print(f'Output folder is {output_folder}\n')
        print('Parsing input file...\n')
        parsed_out, all_genes = _hf.filter_kegg_pathways_genes(filepath=self.input_file_path,
                                                    sheet_name_paths=self.sheet_name_paths,
                                                    sheet_name_genes=self.sheet_name_genes,
                                                    genes_column=self.genes_column,
                                                    log2fc_column=self.log2fc_column,
                                                    count_threshold = count_threshold,
                                                    benjamini_threshold=benjamini_threshold,
                                                    raw_pvalue_threshold=pathway_pvalue)
        if len(parsed_out) == 0:
            raise ValueError("Could not detect pathways in the input file with the selected default & user settings. Check your input file and/or settings.")
        print('Finished parsing input file\n')
        print('Collecting pathway info...')
        pathway_info = await _hf.collect_pathway_info_async(parsed_output=parsed_out)
        print('Finished collecting pathway info\n')
        os.chdir(output_folder)
        _hf.generate_pathways_per_gene_spreadsheet(gene_list=all_genes, pathway_dict=parsed_out , name_extension=None)
        print('Mapping pathways...\n')
        _df.draw_KEGG_pathways_transcripts(parsed_output=parsed_out , info=pathway_info , compounds_list=self.compounds_list , save_to_eps=self.save_to_eps)
        os.chdir(entry_dir)
        print(f'Done! \nOutput files are located in {output_folder}')

    async def multiple_inputs(self, benjamini_threshold: Optional[float] = None, count_threshold: Optional[int] = 2 , pathway_pvalue: Optional[float] = None):
        """
        Perform the Multiple Inputs Analysis.

        Raises:
        TypeError: If the input_file_path is not a list, if the input_label is not a list,
                or if the number of input files does not match the number of labels.

        Prints:
        - Execution message.
        - Output folder path.
        - Information about the number of inputs to be mapped.
        - Parsing and collecting pathway information messages.
        - Completion message.

        Notes:
        - Calls helper functions to filter KEGG pathways for genes, parse input files for multiple interventions,
        and draw KEGG pathways for genes with multiple interventions.
        - The output files are located in the created output folder.
        """
        if not isinstance(self.input_file_path , list):
            raise TypeError('Please provide a list of inputs to perform \'Multiple inputs analysis')
        elif not isinstance(self.input_label , list):
            raise TypeError('Please provide a list with a label for each input file.')
        elif isinstance(self.input_label , list) and isinstance(self.input_file_path , list) and (len(self.input_file_path) != len(self.input_label)):
            raise TypeError('Please make sure that every input file has a corresponding label.')
        print("Executing analysis : Multiple inputs...\n")
        
        entry_dir = os.getcwd()        
        how_many =  len(self.input_label)
        analysis_extension = f'{how_many}_inputs'
        
        folder_of_input = self.find_file_folder()
        output_folder = self.make_output_folder(folder_of_input , analysis_extension)
        print(f'Output folder is {output_folder}\n')
        if how_many > 1:
            print(f"You want to map {how_many} inputs in total.")
        else:
            raise TypeError("Please provide more than one input files to perform this analysis")

        parsed_out_list = []
        all_genes_list = []
        file_counter = 1
        
        print('Parsing input file...\n')
        for (file, inter_name), file_counter in zip(zip(self.input_file_path, self.input_label), range(1, len(self.input_file_path) + 1)):
            print(f"File Counter: {file_counter}, File: {file}, with name {inter_name}\n")
            parsed_out_counter = 'parsed_out_' + str(file_counter)
            all_genes_counter = 'all_genes_' + str(file_counter)
            globals()[parsed_out_counter], globals()[all_genes_counter] =  _hf.filter_kegg_pathways_genes(filepath=file,
                                                                            sheet_name_paths=self.sheet_name_paths,
                                                                            sheet_name_genes=self.sheet_name_genes,
                                                                            genes_column=self.genes_column,
                                                                            log2fc_column=self.log2fc_column,
                                                                            count_threshold = count_threshold , benjamini_threshold=benjamini_threshold,
                                                                            number_interventions=file_counter , name_interventions=inter_name, raw_pvalue_threshold=pathway_pvalue)

            parsed_out_list.append(globals()[parsed_out_counter])
            all_genes_list.append(globals()[all_genes_counter])

            if len(globals()[parsed_out_counter]) == 0:
                print(f"Input file {file_counter} located in {file}, with name {inter_name} did not return any pathways with the selected default & user settings. Check your input file and/or settings.")
            
            file_counter += 1
    
        if len(parsed_out_list) == 0:
            raise ValueError("Could not detect pathways in the input file with the selected default & user settings. Check your input file and/or settings.")
        print('Finished parsing input file\n')

        os.chdir(output_folder)

        for list_counter , (parsed_out_i, all_genes_i) in enumerate(zip(parsed_out_list , all_genes_list)):
            _hf.generate_pathways_per_gene_spreadsheet(gene_list=all_genes_i, pathway_dict=parsed_out_i , name_extension="input"+str(list_counter+1))

        print('Collecting pathway info & mapping pathways...\n')
        await _df.draw_KEGG_pathways_genes_multiple_interventions_async(parsed_out_list=parsed_out_list , intervention_names=self.input_label , colors_list=_cs.colors_list , compounds_list=self.compounds_list , save_to_eps=self.save_to_eps)
        os.chdir(entry_dir)
        print(f'Done! \nOutput files are located in {output_folder}\n')

    async def single_input_with_methylation(self, methylation_path: Union[str, Path], methylation_pvalue: Optional[str] = None, methylation_pvalue_thresh: Optional[float] = None, methylation_genes: str = None,
                                      benjamini_threshold: Optional[float] = None, count_threshold: Optional[int] = 2 , pathway_pvalue: Optional[float] = None):
        """
        Perform Single Input Analysis with Methylation.

        Raises:
        TypeError: If the input_file_path is a list.

        ValueError: If the methylation file path is not provided or is invalid, or if there are no genes with a methylation profile.

        Prints:
        - Execution message.
        - Output folder path.
        - Parsing and collecting pathway information messages.
        - Completion message.

        Notes:
        - Calls helper functions to load and evaluate methylation metadata, filter KEGG pathways for genes with methylation,
        and draw KEGG pathways for genes with methylation.
        - The output files are located in the created output folder.
        """
        if isinstance(self.input_file_path , list):
            raise TypeError('Please provide a single input to perform \'Single input analysis w Methylation')
        
        print("Executing analysis : Single input w Methylation...\n")

        entry_dir = os.getcwd()        
        folder_of_input = self.find_file_folder()
        analysis_extension = 'methylation'
        output_folder = self.make_output_folder(folder_of_input , analysis_extension)
        print(f'Output folder is {output_folder}\n')
        
        if methylation_path is not None or isinstance(methylation_path, (str , os.PathLike)):
            try:
                methylation_df = _hf.load_metadata(methylation_path)
            except ValueError:
                raise ValueError(f'Please provide a proper methylation file path')
            
        _hf.evaluate_metadata(methylation_df , methylation_pvalue , methylation_genes)

        if methylation_pvalue_thresh is None or not isinstance(methylation_pvalue_thresh, (int, float)) or methylation_pvalue is None:
            genes_from_MM = methylation_df[methylation_genes].unique().tolist()
        else:
            if methylation_pvalue is not None and methylation_pvalue not in methylation_df.columns:
                raise KeyError(f'Column {methylation_pvalue} not found in the methylation dataframe.')
            
            try:
                genes_from_MM = methylation_df.loc[methylation_df[methylation_pvalue] < methylation_pvalue_thresh][methylation_genes].unique().tolist()
            except ValueError:
                raise ValueError(f'Invalid value provided for pvalue_thresh. It should be a number.')

        if len(genes_from_MM) == 0:
            raise ValueError('There are no genes with a methylation profile')
        
        methylation_options = ['DEG with DMP' , 'DEG without DMP']
        color_to_methylation = { meth : color for (meth , color) in zip(methylation_options , _cs.colors_list)}
        print('Parsing input file...\n')
        parsed_out, all_genes = _hf.filter_kegg_pathways_genes(filepath=self.input_file_path,
                                                    sheet_name_paths=self.sheet_name_paths,
                                                    sheet_name_genes=self.sheet_name_genes,
                                                    genes_column=self.genes_column,
                                                    log2fc_column=self.log2fc_column,
                                                    count_threshold = count_threshold,
                                                    benjamini_threshold=benjamini_threshold,
                                                    raw_pvalue_threshold=pathway_pvalue)
        if len(parsed_out) == 0:
            raise ValueError("Could not detect pathways in the input file with the selected default & user settings. Check your input file and/or settings.")
        print('Finished parsing input file\n')
        print('Collecting pathway info...')
        pathway_info = await _hf.collect_pathway_info_async(parsed_output=parsed_out)
        print('Finished collecting pathway info\n')
        os.chdir(output_folder)
        _hf.generate_pathways_per_gene_spreadsheet(gene_list=all_genes, pathway_dict=parsed_out , name_extension=None)
        print('Mapping pathways...\n')
        _df.draw_KEGG_pathways_genes_with_methylation(parsed_output=parsed_out , info=pathway_info , genes_from_MM=genes_from_MM , color_legend=color_to_methylation , compounds_list=self.compounds_list , save_to_eps=self.save_to_eps)
        os.chdir(entry_dir)
        print(f'Done! \nOutput files are located in {output_folder}')

    async def single_input_with_miRNA(self, miRNA_path: Union[str, Path] = None, miRNA_pvalue: Optional[str] = None, miRNA_pvalue_thresh: Optional[float] = None, miRNA_genes: str = None, 
                                benjamini_threshold: Optional[float] = None, count_threshold: Optional[int] = 2 , pathway_pvalue: Optional[float] = None):
        """
        Perform Single Input Analysis with miRNA.

        Raises:
        TypeError: If the input_file_path is a list.

        ValueError: If the miRNA file path is not provided or is invalid, or if there are no genes with a miRNA profile.

        Prints:
        - Execution message.
        - Output folder path.
        - Parsing and collecting pathway information messages.
        - Completion message.

        Notes:
        - Calls helper functions to load and evaluate miRNA metadata, filter KEGG pathways for genes with miRNA,
        and draw KEGG pathways for genes with miRNA.
        - The output files are located in the created output folder.
        """
        if isinstance(self.input_file_path , list):
            raise TypeError('Please provide a single input to perform \'Single input analysis w miRNA')
        
        print("Executing analysis : Single input w miRNA...\n")

        entry_dir = os.getcwd()        
        folder_of_input = self.find_file_folder()    
        analysis_extension = 'miRNA'
        output_folder = self.make_output_folder(folder_of_input , analysis_extension)
        print(f'Output folder is {output_folder}\n')

        if miRNA_path is not None or isinstance(miRNA_path, (str , os.PathLike)):
            try:
                miRNA_df = _hf.load_metadata(miRNA_path)
            except ValueError:
                raise ValueError(f'Please provide a proper miRNA file path')

        _hf.evaluate_metadata(miRNA_df , miRNA_pvalue , miRNA_genes)

        if miRNA_pvalue_thresh is None or not isinstance(miRNA_pvalue_thresh, (int, float)) or miRNA_pvalue is None:
            genes_from_miRNA = miRNA_df[miRNA_genes].unique().tolist()
        else:
            if miRNA_pvalue is not None and miRNA_pvalue not in miRNA_df.columns:
                raise KeyError(f'Column {miRNA_pvalue} not found in the miRNA dataframe.')
            
            try:
                genes_from_miRNA = miRNA_df.loc[miRNA_df[miRNA_pvalue] < miRNA_pvalue_thresh][miRNA_genes].unique().tolist()
            except ValueError:
                raise ValueError(f'Invalid value provided for pvalue_thresh. It should be a number.')

        if len(genes_from_miRNA) == 0:
            raise ValueError('There are no genes with a miRNA profile')

        miRNA_options = ['DEmiR target gene' , 'No DEmiR target gene']
        color_to_miRNA = {miRNA : color for (miRNA , color) in zip(miRNA_options , _cs.colors_list)}
        print('Parsing input file...\n')
        parsed_out, all_genes = _hf.filter_kegg_pathways_genes(filepath=self.input_file_path,
                                                    sheet_name_paths=self.sheet_name_paths,
                                                    sheet_name_genes=self.sheet_name_genes,
                                                    genes_column=self.genes_column,
                                                    log2fc_column=self.log2fc_column,
                                                    count_threshold = count_threshold,
                                                    benjamini_threshold=benjamini_threshold,
                                                    raw_pvalue_threshold=pathway_pvalue)
        if len(parsed_out) == 0:
            raise ValueError("Could not detect pathways in the input file with the selected default & user settings. Check your input file and/or settings.")        
        print('Finished parsing input file\n')
        print('Collecting pathway info...')
        pathway_info = await _hf.collect_pathway_info_async(parsed_output=parsed_out)
        print('Finished collecting pathway info\n')
        os.chdir(output_folder)
        _hf.generate_pathways_per_gene_spreadsheet(gene_list=all_genes, pathway_dict=parsed_out , name_extension=None)
        print('Mapping pathways...\n')
        _df.draw_KEGG_pathways_genes_with_miRNA(parsed_output=parsed_out , info=pathway_info , genes_from_miRNA=genes_from_miRNA , color_legend=color_to_miRNA , compounds_list=self.compounds_list , save_to_eps=self.save_to_eps)
        os.chdir(entry_dir)
        print(f'Done! \nOutput files are located in {output_folder}')

    async def single_input_with_methylation_and_miRNA(self, methylation_path: Union[str, Path], methylation_pvalue: Optional[str] = None, methylation_pvalue_thresh: Optional[float] = None, methylation_genes: str = None,
                                                miRNA_path: Union[str, Path] = None, miRNA_pvalue: Optional[str] = None, miRNA_pvalue_thresh: Optional[float] = None, miRNA_genes: str = None,
                                                benjamini_threshold: Optional[float] = None, count_threshold: Optional[int] = 2 , pathway_pvalue: Optional[float] = None):
        """
        Perform a single input analysis with Methylation and miRNA data.

        Raises:
            TypeError: If input_file_path is a list.
            ValueError: If there are issues with loading methylation or miRNA metadata.
                        If invalid values are provided for pvalue_thresh.
            KeyError: If a specified column is not found in the metadata dataframe.

        Prints:
            Execution message.
            Output folder location.

        Returns:
            None. Results are saved in the output folder.
        """
        if isinstance(self.input_file_path , list):
            raise TypeError('Please provide a single input to perform \'Single input analysis w Methylation & miRNA')
        
        print("Executing analysis : Single input w Methylation & miRNA...\n")

        entry_dir = os.getcwd()        
        folder_of_input = self.find_file_folder()    
        analysis_extension = 'methylation_and_miRNA'
        output_folder = self.make_output_folder(folder_of_input , analysis_extension)
        print(f'Output folder is {output_folder}\n')

        if methylation_path is not None or isinstance(methylation_path, (str , os.PathLike)):
            try:
                methylation_df = _hf.load_metadata(methylation_path)
            except ValueError:
                raise ValueError(f'Please provide a proper methylation file path')


        _hf.evaluate_metadata(methylation_df , methylation_pvalue , methylation_genes)

        if methylation_pvalue_thresh is None or not isinstance(methylation_pvalue_thresh, (int, float)) or methylation_pvalue is None:
            genes_from_MM = methylation_df[methylation_genes].unique().tolist()
        else:
            if methylation_pvalue is not None and methylation_pvalue not in methylation_df.columns:
                raise KeyError(f'Column {methylation_pvalue} not found in the methylation dataframe.')
            
            try:
                genes_from_MM = methylation_df.loc[methylation_df[methylation_pvalue] < methylation_pvalue_thresh][methylation_genes].unique().tolist()
            except ValueError:
                raise ValueError(f'Invalid value provided for pvalue_thresh. It should be a number.')

        if len(genes_from_MM) == 0:
            raise ValueError('There are no genes with a methylation profile')


        if miRNA_path is not None or isinstance(miRNA_path, (str , os.PathLike)):
            try:
                miRNA_df = _hf.load_metadata(miRNA_path)
            except ValueError:
                raise ValueError(f'Please provide a proper miRNA file path')
            

        _hf.evaluate_metadata(miRNA_df , miRNA_pvalue , miRNA_genes)

        if miRNA_pvalue_thresh is None or not isinstance(miRNA_pvalue_thresh, (int, float)) or miRNA_pvalue is None:
            genes_from_miRNA = miRNA_df[miRNA_genes].unique().tolist()
        else:
            if miRNA_pvalue is not None and miRNA_pvalue not in miRNA_df.columns:
                raise KeyError(f'Column {miRNA_pvalue} not found in the miRNA dataframe.')
            
            try:
                genes_from_miRNA = miRNA_df.loc[miRNA_df[miRNA_pvalue] < miRNA_pvalue_thresh][miRNA_genes].unique().tolist()
            except ValueError:
                raise ValueError(f'Invalid value provided for pvalue_thresh. It should be a number.')

        if len(genes_from_miRNA) == 0:
            raise ValueError('There are no genes with a miRNA profile')
        

        methylation_w_miRNA_options = ['DEG with DMP and DEmiR target gene', 'DEG without DMP and DEmiR target gene',
                                       'DEG with DMP and no DEmiR target gene', 'DEG without DMP and no DEmiR target gene']
        color_to_methylation_w_miRNA = { meth_miRNA : color for (meth_miRNA , color) in zip(methylation_w_miRNA_options , _cs.colors_list)}
        print('Parsing input file...\n')
        parsed_out, all_genes = _hf.filter_kegg_pathways_genes(filepath=self.input_file_path,
                                                    sheet_name_paths=self.sheet_name_paths,
                                                    sheet_name_genes=self.sheet_name_genes,
                                                    genes_column=self.genes_column,
                                                    log2fc_column=self.log2fc_column,
                                                    count_threshold = count_threshold,
                                                    benjamini_threshold=benjamini_threshold,
                                                    raw_pvalue_threshold=pathway_pvalue)
        if len(parsed_out) == 0:
            raise ValueError("Could not detect pathways in the input file with the selected default & user settings. Check your input file and/or settings.")        
        print('Finished parsing input file\n')
        print('Collecting pathway info...')
        pathway_info = await _hf.collect_pathway_info_async(parsed_output=parsed_out)
        print('Finished collecting pathway info\n')
        os.chdir(output_folder)
        _hf.generate_pathways_per_gene_spreadsheet(gene_list=all_genes, pathway_dict=parsed_out , name_extension=None)
        print('Mapping pathways...\n')
        _df.draw_KEGG_pathways_genes_with_methylation_and_miRNA(parsed_output=parsed_out , info=pathway_info ,
                                                                genes_from_MM=genes_from_MM , genes_from_miRNA=genes_from_miRNA,
                                                                color_legend=color_to_methylation_w_miRNA , compounds_list=self.compounds_list , save_to_eps=self.save_to_eps)
        os.chdir(entry_dir)
        print(f'Done! \nOutput files are located in {output_folder}')

    async def single_input_genes_bulk_mapping(self):
        """
        Perform a single input analysis with bulk mapping for genes.

        Raises:
            TypeError: If input_file_path is a list.

        Prints:
            Execution message.
            Output folder location.

        Returns:
            None. Results are saved in the output folder.
        """
        if isinstance(self.input_file_path , list):
            raise TypeError('Please provide a single input to perform \'Single input (Bulk mapping)')

        print("Executing analysis : Single input (Bulk mapping)...\n")

        entry_dir = os.getcwd()        
        folder_of_input = self.find_file_folder()    
        analysis_extension = 'bulk'
        output_folder = self.make_output_folder(folder_of_input , analysis_extension)
        print(f'Output folder is {output_folder}')
        print('Parsing input file...\n')
        parsed_out, all_genes = _hf.parse_bulk_kegg_pathway_file(filepath=self.input_file_path,
                                                    sheet_name_paths=self.sheet_name_paths,
                                                    sheet_name_genes=self.sheet_name_genes,
                                                    genes_column=self.genes_column,
                                                    log2fc_column=self.log2fc_column)
        if len(parsed_out) == 0:
            raise ValueError("Could not detect pathways in the input file with the selected default & user settings. Check your input file and/or settings.")        
        print('Finished parsing input file\n')
        print('Collecting pathway info...')
        pathway_info = await _hf.collect_pathway_info_async(parsed_output=parsed_out)
        print('Finished collecting pathway info\n')
        os.chdir(output_folder)
        _hf.generate_pathways_per_gene_spreadsheet(gene_list=all_genes, pathway_dict=parsed_out , name_extension=None)
        print('Mapping pathways...\n')
        _df.draw_KEGG_pathways_genes(parsed_output=parsed_out , info=pathway_info , compounds_list=self.compounds_list , save_to_eps=self.save_to_eps)
        os.chdir(entry_dir)
        print(f'Done! \nOutput files are located in {output_folder}')

    async def single_input_with_miRNA_quantification(self, miRNA_path: Union[str, Path] = None, miRNA_pvalue: Optional[str] = None, miRNA_pvalue_thresh: Optional[float] = None, miRNA_genes: str = None, miRNA_ID_column: Optional[str] = None,
                                               benjamini_threshold: Optional[float] = None, count_threshold: Optional[int] = 2 , pathway_pvalue: Optional[float] = None):
        """
        Perform Single Input Analysis with miRNA.

        Raises:
        TypeError: If the input_file_path is a list.

        ValueError: If the miRNA file path is not provided or is invalid, or if there are no genes with a miRNA profile.

        Prints:
        - Execution message.
        - Output folder path.
        - Parsing and collecting pathway information messages.
        - Completion message.

        Notes:
        - Calls helper functions to load and evaluate miRNA metadata, filter KEGG pathways for genes with miRNA,
        and draw KEGG pathways for genes with miRNA.
        - The output files are located in the created output folder.
        """
        if isinstance(self.input_file_path , list):
            raise TypeError('Please provide a single input to perform \'Single input analysis w miRNA quantification')
        
        print("Executing analysis : Single input w miRNA quantification...\n")

        entry_dir = os.getcwd()        
        folder_of_input = self.find_file_folder()    
        analysis_extension = 'miRNA_quantification'
        output_folder = self.make_output_folder(folder_of_input , analysis_extension)
        print(f'Output folder is {output_folder}\n')

        if miRNA_path is not None or isinstance(miRNA_path, (str , os.PathLike)):
            try:
                miRNA_df = _hf.load_metadata(miRNA_path)
            except ValueError:
                raise ValueError(f'Please provide a proper miRNA file path')

        _hf.evaluate_metadata(miRNA_df , miRNA_pvalue , miRNA_genes)

        if miRNA_pvalue_thresh is None or not isinstance(miRNA_pvalue_thresh, (int, float)) or miRNA_pvalue is None:
            genes_from_miRNA = miRNA_df[miRNA_genes].unique().tolist()
        else:
            if miRNA_pvalue is not None and miRNA_pvalue not in miRNA_df.columns:
                raise KeyError(f'Column {miRNA_pvalue} not found in the miRNA dataframe.')
            
            try:
                genes_from_miRNA = miRNA_df.loc[miRNA_df[miRNA_pvalue] < miRNA_pvalue_thresh][miRNA_genes].unique().tolist()
            except ValueError:
                raise ValueError(f'Invalid value provided for pvalue_thresh. It should be a number.')

        if len(genes_from_miRNA) == 0:
            raise ValueError('There are no genes with a miRNA profile')

        if miRNA_ID_column is None:
            raise KeyError(f'Please provide the column name with the miRNA IDs')
        elif miRNA_ID_column is not None and miRNA_ID_column not in miRNA_df.columns:
            raise KeyError(f'Column {miRNA_ID_column} not found in the miRNAs dataframe.')


        miRNA_options = ['DEmiR target gene' , 'No DEmiR target gene']
        color_to_miRNA = {miRNA : color for (miRNA , color) in zip(miRNA_options , _cs.colors_list)}
        print('Parsing input file...\n')
        parsed_out, all_genes = _hf.filter_kegg_pathways_genes(filepath=self.input_file_path,
                                                    sheet_name_paths=self.sheet_name_paths,
                                                    sheet_name_genes=self.sheet_name_genes,
                                                    genes_column=self.genes_column,
                                                    log2fc_column=self.log2fc_column,
                                                    count_threshold = count_threshold,
                                                    benjamini_threshold=benjamini_threshold,
                                                    raw_pvalue_threshold=pathway_pvalue)
        if len(parsed_out) == 0:
            raise ValueError("Could not detect pathways in the input file with the selected default & user settings. Check your input file and/or settings.")        
        print('Finished parsing input file\n')
        print('Collecting pathway info...')
        pathway_info = await _hf.collect_pathway_info_async(parsed_output=parsed_out)
        print('Finished collecting pathway info\n')
        os.chdir(output_folder)
        _hf.generate_pathways_per_gene_spreadsheet(gene_list=all_genes, pathway_dict=parsed_out , name_extension=None)
        print('Mapping pathways...\n')
        _df.draw_KEGG_pathways_genes_with_miRNA_quantification(parsed_output=parsed_out , info=pathway_info , genes_from_miRNA=genes_from_miRNA , miRNA_df=miRNA_df , miRNA_genes_col = miRNA_genes , miRNA_id_col=miRNA_ID_column , compounds_list=self.compounds_list ,save_to_eps=self.save_to_eps)
        os.chdir(entry_dir)
        print(f'Done! \nOutput files are located in {output_folder}')

    async def single_input_with_methylation_quantification(self, methylation_path: Union[str, Path], methylation_pvalue: Optional[str] = None, methylation_pvalue_thresh: Optional[float] = None, methylation_genes: str = None,
                                                     probes_to_cgs: Optional[bool] = False, methylation_probe_column: Optional[str] = None,
                                                     benjamini_threshold: Optional[float] = None, count_threshold: Optional[int] = 2 , pathway_pvalue: Optional[float] = None):
        """
        Perform Single Input Analysis with methylation quantification.

        Raises:
        TypeError: If the input_file_path is a list.

        ValueError: If the methylation file path is not provided or is invalid, or if there are no genes with a methylation profile.

        Prints:
        - Execution message.
        - Output folder path.
        - Parsing and collecting pathway information messages.
        - Completion message.

        Notes:
        - Calls helper functions to load and evaluate methylation metadata, filter KEGG pathways for genes with methylation,
        and draw KEGG pathways for genes with methylation.
        - The output files are located in the created output folder.
        """
        if isinstance(self.input_file_path , list):
            raise TypeError('Please provide a single input to perform \'Single input analysis w methylation quantification')
        
        print("Executing analysis : Single input w methylation quantification...\n")

        entry_dir = os.getcwd()        
        folder_of_input = self.find_file_folder()    
        if probes_to_cgs:
            analysis_extension = 'methylation_quantification_probe_correction'
            print("Selected setting : Probe correction\n")
        else:
            analysis_extension = 'methylation_quantification'

        output_folder = self.make_output_folder(folder_of_input , analysis_extension)
        print(f'Output folder is {output_folder}\n')

        if methylation_path is not None or isinstance(methylation_path, (str , os.PathLike)):
            try:
                methylation_df = _hf.load_metadata(methylation_path)
            except ValueError:
                raise ValueError(f'Please provide a proper methylation file path')

        _hf.evaluate_metadata(methylation_df , methylation_pvalue , methylation_genes)

        if methylation_pvalue_thresh is None or not isinstance(methylation_pvalue_thresh, (int, float)) or methylation_pvalue is None:
            genes_from_methylation = methylation_df[methylation_genes].unique().tolist()
        else:
            if methylation_pvalue is not None and methylation_pvalue not in methylation_df.columns:
                raise KeyError(f'Column {methylation_pvalue} not found in the methylation dataframe.')
            
            try:
                genes_from_methylation = methylation_df.loc[methylation_df[methylation_pvalue] < methylation_pvalue_thresh][methylation_genes].unique().tolist()
            except ValueError:
                raise ValueError(f'Invalid value provided for pvalue_thresh. It should be a number.')

        if len(genes_from_methylation) == 0:
            raise ValueError('There are no genes with a genes_from_methylation profile')

        if methylation_probe_column is None:
            raise KeyError(f'Please provide the column name with the probe IDs.')
        elif methylation_probe_column is not None and methylation_probe_column not in methylation_df.columns:
            raise KeyError(f'Column {methylation_probe_column} not found in the methylation dataframe.')
        else:
            metadata_id_col = methylation_probe_column
            methylation_df = methylation_df.drop_duplicates(subset=[metadata_id_col , methylation_genes] , keep='first')
            if probes_to_cgs is True:
                insert_new_column = 'unique_CG_quantification'
                if insert_new_column in methylation_df.columns:
                    raise KeyError(f'Could not insert unique probe column ({insert_new_column}). It already exists in the dataframe.')
                else:
                    methylation_df[insert_new_column] = methylation_df[methylation_probe_column].str.split("_").str[0]
                    methylation_df = methylation_df.drop_duplicates(subset=[insert_new_column, methylation_genes] , keep='first')
                    metadata_id_col = insert_new_column
            elif probes_to_cgs is False:
                print('Will not perform probe correction...\n')

        methylation_options = ['DEG with DMP' , 'DEG without DMP']
        color_to_methylation = { meth : color for (meth , color) in zip(methylation_options , _cs.colors_list)}
        print('Parsing input file...\n')
        parsed_out, all_genes = _hf.filter_kegg_pathways_genes(filepath=self.input_file_path,
                                                    sheet_name_paths=self.sheet_name_paths,
                                                    sheet_name_genes=self.sheet_name_genes,
                                                    genes_column=self.genes_column,
                                                    log2fc_column=self.log2fc_column,
                                                    count_threshold = count_threshold,
                                                    benjamini_threshold=benjamini_threshold,
                                                    raw_pvalue_threshold=pathway_pvalue)
        if len(parsed_out) == 0:
            raise ValueError("Could not detect pathways in the input file with the selected default & user settings. Check your input file and/or settings.")        
        print('Finished parsing input file\n')
        print('Collecting pathway info...')
        pathway_info = await _hf.collect_pathway_info_async(parsed_output=parsed_out)
        print('Finished collecting pathway info\n')
        os.chdir(output_folder)
        _hf.generate_pathways_per_gene_spreadsheet(gene_list=all_genes, pathway_dict=parsed_out , name_extension=None)
        print('Mapping pathways...\n')
        _df.draw_KEGG_pathways_genes_with_methylation_quantification(parsed_output=parsed_out , info=pathway_info , genes_from_MM=genes_from_methylation , MM_df=methylation_df , MM_genes_col = methylation_genes , MM_id_col=metadata_id_col , compounds_list=self.compounds_list , save_to_eps=self.save_to_eps)
        os.chdir(entry_dir)
        print(f'Done! \nOutput files are located in {output_folder}')
