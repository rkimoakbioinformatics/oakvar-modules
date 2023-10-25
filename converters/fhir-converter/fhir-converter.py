import uuid
import os
import pathlib
from typing import Any
from typing import Tuple
from typing import List
from typing import Dict
import ijson
import pandas as pd
import requests
import hashlib
import sqlite3
from oakvar import BaseConverter




class Converter(BaseConverter):
    def __init__(self, input_path, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.format_name = 'FHIR' or 'fhir'
        self.resources = None
        self.filename = pathlib.Path(input_path)
     


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
        with open(filename, 'r') as fn:
            parser = ijson.parse(fn)
            state = None
            for prefix,event,value in parser: 
                if event == 'map_key': 
                    state = value 
                    if state == 'resourceType':
                        return True 


    def refseq_to_USCS(self,chromosome_mappings):
        refseq_2_uscs_dict = {}
        with open(chromosome_mappings,'r') as file:
            for line in file: 
                key, value  = line.split('c')
                updated_key = key.split(".")[0]
                refseq_2_uscs_dict[updated_key.strip()] = value.strip()

        for item in refseq_2_uscs_dict.keys():
            refseq_2_uscs_dict[item] = "c" + refseq_2_uscs_dict[item]

        return refseq_2_uscs_dict
    def chrom_finder(self,db,transcript,dict) -> str: 
        #trim transcript if not in correct format
        str_tr = transcript.split(':')[0]

        final_str = str_tr.split(".")[0]

        if final_str[0] != "N":
            conn = sqlite3.connect(db)
            curs = conn.cursor()

            curs.execute(f'SELECT "chromid" from "transcript" WHERE "name" = "{str_tr}"')


            chrom_id = curs.fetchall()[0][0]
            #print(chrom_id)

            curs.execute(f'SELECT "chrom" from "chroms" WHERE "chromid" = {chrom_id}')
            chrom = curs.fetchall()[0][0]

            conn.close()
            return chrom
        else: 
            return dict[final_str]



    def convert_file(self, file):
        chr_dict = self.refseq_to_USCS("GRCh38_RefSeq2UCSC.txt")
        gencode_db = ""
        start_terms = ['resourceType','name','resource','variant']
        var_dicts = []
        var_dict = {
            "chrom": None ,
            "pos":  None,
            "ref_base": None ,
            "alt_base": None,
            "sample_id": None,
        }
        ref_or_alt = None 


        allele_codes = ["69547-8","69551-0"]
        HGVS_codes = ["48013-7","51958-7",'48004-6']
        pos_codes = ["81254-5"]


        in_hgvs_comp = False
        in_pos_comp = False

        hgvs_systems = ["http://varnomen.hgvs.org",'https://www.ncbi.nlm.nih.gov/dbvar/','https://api.ncbi.nlm.nih.gov/variation/v0/',
                        "http://www.ensembl.org"]




        with open(file, 'r') as json_file:
            parser = ijson.parse(json_file)
            for prefix,event,value in parser:


                #check to see if in
                if event == 'map_key':
                    current_key = value
                    if current_key in start_terms:
                            var_dict = {
                                "chrom": None ,
                                "pos":  None,
                                "ref_base": None ,
                                "alt_base": None,
                                "sample_id": None,
                            }

                #look at component value
                elif event =='string' or event == 'number':

                    #check to see if in hgvs system
                    if current_key == 'system':
                        if value in hgvs_systems:
                            #print(value)
                            in_hgvs_comp = True

                    # check to see if in a position component to extract pos 
                    if in_pos_comp and current_key == 'value':
                        in_pos_comp = False
                        var_dict['pos'] = value

                    if current_key == 'code':

                        if in_hgvs_comp:
                            #print(value)
                            in_hgvs_comp = False
                            try:
                                var_dict['chrom'] = self.chrom_finder(gencode_db,value,chr_dict)
                            except: var_dict['chrom'] = value

                        if value in HGVS_codes:
                            in_hgvs_comp = True
                        if value in pos_codes:
                            in_pos_comp = True 

                        if value == "69547-8":
                            ref_or_alt = 'Ref'
                        if value == "69551-0":
                            ref_or_alt = 'Alt'

                    if current_key == 'valueString' and ref_or_alt == "Ref":
                        var_dict['ref_base'] = value

                    if current_key == 'valueString' and ref_or_alt == 'Alt':

                        var_dict['alt_base'] = value

                        var_dicts.append(var_dict)



                elif event == 'end_map':
                    pass
            return var_dicts