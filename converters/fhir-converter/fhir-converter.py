import uuid
import os
import pathlib
import ijson
import glob
from typing import Any
from typing import Tuple
from typing import List
from typing import Dict
import pandas as pd
import requests
import hashlib
import sqlite3
from oakvar import BaseConverter
import oakvar as ov 




class Converter(BaseConverter):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.format_name = 'FHIR' or 'fhir'
        self.resources = None
        self.filename = pathlib.Path(self.input_path)

    def setup(self, *args, **kwargs):
        fhir_dir = ov.lib.module.local.get_module_dir("fhir-converter")   
        mapping_file = str(fhir_dir) + "\GRCh38_RefSeq2UCSC.txt"

        gencode_dir = ov.lib.module.local.get_module_dir("gencode")
        gen_data = str(gencode_dir) + "\data"
        path_data = pathlib.Path(gen_data)
        file_list = glob.glob(os.path.join(path_data, f"*.sqlite"))
        self.gencode_db = ""
        self.highest_gen_version = 0
        def version_finder(a_file):
            parts = a_file.split(".")
            prefix = parts[0]
            numbers = prefix.split("_")
            return numbers

        for file in file_list:
            name = os.path.basename(file)
            version =version_finder(name)
            if int(version[1]) > int(self.highest_gen_version):
                self.gencode_db = file
                self.highest_gen_version = version[1]


        self.chr_dict = self.refseq_to_USCS(mapping_file)

        self.start_terms = ['resourceType','name','resource','variant']
        self.var_dict = {
             "chrom": None ,
             "pos":  None,
             "ref_base": None ,
             "alt_base": None,
             "sample_id": None,
         }
        self.ref_or_alt = None 
        self.allele_codes = ["69547-8","69551-0"]
        self.HGVS_codes = ["48013-7","51958-7",'48004-6']
        self.pos_codes = ["81254-5"]
        self.in_hgvs_comp = False
        self.in_pos_comp = False
        self.hgvs_systems = ["http://varnomen.hgvs.org",'https://www.ncbi.nlm.nih.gov/dbvar/','https://api.ncbi.nlm.nih.gov/variation/v0/',
                        "http://www.ensembl.org"]
        json_file = open(self.input_path, 'r')
        self.parser = ijson.parse(json_file)
     


    def check_format(self,input_path, *__args__, **__kwargs__)-> bool:
        """
        Detect the format of an input file.

        Arguments:
            f: a file handle to an input file
        Returns:
            bool: True if the input file is for this converter,
                  False if not.

        The example below checks if the input file's first line indicates
        VCF file format.
        """
        filename = pathlib.Path(input_path)
        print(input_path)   
        with open(filename, 'r') as fn:
            parser = ijson.parse(fn)
            state = None
            for prefix,event,value in parser: 
                if event == 'map_key': 
                    state = value 
                    if state == 'resourceType':
                        return True 


    # If your converter module needs something else than
    # the standard way of opening a text input file,
    # read line by line, and coverting each line into
    # a list of dictionaries of variants,
    # you may want to start with modifying
    # convert_file method. In that case, uncomment
    # the below convert_file method and add your implementation.
    def refseq_to_USCS(self,chromosome_mappings):
        refseq_2_uscs_dict = {}
        with open(chromosome_mappings,'r') as file:
            for line in file: 
                key, value  = line.split('c')
                updated_key = key.split(".")[0]
                update_value = value.split("_")[0]
                refseq_2_uscs_dict[updated_key.strip()] = update_value.strip()

        for item in refseq_2_uscs_dict.keys():
            refseq_2_uscs_dict[item] = "c" + refseq_2_uscs_dict[item]

        return refseq_2_uscs_dict
    def chrom_finder(self,db,transcript,dict) -> str: 
        #trim transcript if not in correct format
        str_tr = transcript.split(':')[0]

        final_str = str_tr.split(".")[0]

        if final_str[0] == "E":
            conn = sqlite3.connect(db)
            curs = conn.cursor()

            curs.execute(f'SELECT "chromid" from "transcript" WHERE "name" = "{str_tr}"')

            chrom_id = curs.fetchall()[0][0]

            curs.execute(f'SELECT "chrom" from "chroms" WHERE "chromid" = {chrom_id}')
            chrom = curs.fetchall()[0][0]

            conn.close()
            return chrom
        else:
            try: dict[final_str] 
            except: final_str
    def dict_checker(self,dict):
        for key in dict:
            if key is not 'sample_id' and dict[key] is None:
                return False
        return True



    def get_variant_lines(self,input_path, num_pool, start_line_no, batch_size):
        lines = []
        component = []
        immature_exit = True
        in_pos = False
        var_dict = {
            "chrom": None ,
            "pos":  None,
            "ref_base": None ,
            "alt_base": None,
            "sample_id": None,
        }
        with open(input_path, 'r') as json_file:
                parser = ijson.parse(json_file)
                for prefix,event,value in parser:
                    if self.dict_checker(var_dict):
                        lines.append(var_dict)
                        var_dict = {
                            "chrom": None ,
                            "pos":  None,
                            "ref_base": None ,
                            "alt_base": None,
                            "sample_id": None,
                        }


                    #check to see if in
                    if event == 'map_key':
                        current_key = value
                    #    if current_key in self.start_terms:


                    #look at component value
                    elif event =='string' or event == 'number':
                        if self.in_pos_comp is True:
                            if len(component) > 2: 
                                var_dict["pos"] = component[1]['value']
                                component = []
                                self.in_pos_comp = False
                            else:
                                component.append({current_key:value})

                        #check to see if in hgvs system
                        if current_key == 'system':
                            if value in self.hgvs_systems:
                                self.in_hgvs_comp = True

                        # check to see if in a position component to extract pos 
                        if self.in_pos_comp and current_key == 'value':
                            self.in_pos_comp = False
                            var_dict['pos'] = value

                        if current_key == 'code':

                            if self.in_hgvs_comp:
                                if self.chrom_finder(self.gencode_db,value,self.chr_dict) is None: 
                                    self.in_hgvs_comp = False
                                try:
                                    var_dict['chrom'] = self.chrom_finder(self.gencode_db,value,self.chr_dict)
                                except: var_dict['chrom'] = value

                            if value in self.HGVS_codes:
                                self.in_hgvs_comp = True
                            if value in self.pos_codes:
                                self.in_pos_comp = True 

                            if value == "69547-8":
                                self.ref_or_alt = 'Ref'
                            if value == "69551-0":
                                self.ref_or_alt = 'Alt'

                        if current_key == 'valueString' and self.ref_or_alt == "Ref":
                            var_dict['ref_base'] = value

                        if current_key == 'valueString' and self.ref_or_alt == 'Alt':

                            var_dict['alt_base'] = value

                    elif event == 'end_map':
                        pass
        return lines
    def convert_line(self, line):
        return line
