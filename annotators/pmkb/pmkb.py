from typing import Optional
from oakvar import BaseAnnotator
import re
import sqlite3


class Annotator(BaseAnnotator):
    
    def annotate(self, input_data: dict, secondary_data: Optional[dict] = None):
        assert input_data is not None
        if not self.cursor:
            return None
        #get data for querying the pmkb database from input_data
        gene = input_data['hugo']
        transcript = input_data['transcript']
        exonno = input_data['exonno']
        if input_data['achange'] != None or input_data['achange'] != '':
            achange = input_data['achange']
        qr = []
        query = """ 
                    SELECT 
                        interpretations_final.gene_name, interpretations_final.tumor_type,interpretations_final.tissue_type,
                        interpretations_final.pmkb_url_interpretations,
                        interpretations_final.interpretations, interpretations_final.citations,
                        variant.achange, variant.pmkb_url_variants 
                    FROM 
                        interpretations_final
                    LEFT JOIN
                        variant
                        ON variant.gene = interpretations_final.gene_name
                    WHERE 
                        achange = :achange_pmkb AND variant.gene = :gene
                """
        # exonno = input_data['exonno']
        variant_type = input_data['so']
        #Handle missense variants
        if variant_type == 'MIS' or variant_type == 'missense_variant':
            #search for the description of missense variant
            input_ref_alt_pos = re.search(r'(\w{3})(\d+)(\w{3}|\?)',achange)
            if input_ref_alt_pos is not None:
                ref_alt_pos_catch = input_ref_alt_pos.groups()
                input_pos = ref_alt_pos_catch[1]
                input_ref_allele = seq1(ref_alt_pos_catch[0])
                input_alt_allele = seq1(ref_alt_pos_catch[2])
                #Query the database for the possible variants
                self.cursor.execute("""
                    SELECT 
                        achange
                    FROM 
                        variant
                    WHERE 
                        gene = :gene 
                """, {"gene": gene})
                #save results in pmkb_variants variable    
                pmkb_variants = self.cursor.fetchall()
                #loop through results to get a match
                for line in pmkb_variants:
                    pmkb_achange = line[0].split(':')
                    #get pos:ref_allele_alt_allele
                    pmkb_pos = pmkb_achange[1]
                    pmkb_ref_allele = pmkb_achange[-2]
                    pmkb_alt_allele = pmkb_achange[-1]
                    pmkb_so = pmkb_achange[-3]
                    pmkb_kind = pmkb_achange[0]
                    #match results of pmkb achange with input data
                    if input_pos in pmkb_pos:
                        if  pmkb_ref_allele == input_ref_allele and pmkb_alt_allele == input_alt_allele:
                            self.cursor.execute(f'{query}', {"achange_pmkb": line[0],"gene": gene})
                            qr.append(self.cursor.fetchone())
                        #if positions meet 
                        elif pmkb_ref_allele == '_any' and pmkb_alt_allele == '_any':
                            self.cursor.execute(query, {"achange_pmkb": line[0],"gene": gene})
                            qr.append(self.cursor.fetchone())
                        if pmkb_kind == '_exon':
                            if pmkb_pos == exonno:
                                self.cursor.execute(query,{"achange_pmkb": line[0],"gene": gene})
                                qr.append(self.cursor.fetchone())
                    #handle _any kind
                    if pmkb_so == '_any' or pmkb_so == 'any':
                        if pmkb_pos == '_any':
                            self.cursor.execute(query, {"achange_pmkb": line[0],"gene": gene})
                            qr.append(self.cursor.fetchone())
        
        #handle frameshift mutations
        elif  variant_type == 'FSD' or variant_type  == 'FSI' or variant_type == 'frameshift_insertion' or variant_type == 'frameshift_deletion': 
            #get input data for querying the pmkb database
            input_ref_alt_pos = re.search(r'(\w{3})(\d+)',achange)
            if input_ref_alt_pos is not None:
                ref_alt_pos_catch = input_ref_alt_pos.groups()
                input_pos = ref_alt_pos_catch[1]
                input_ref_allele = seq1(ref_alt_pos_catch[0])
                #query the pmkb database
                self.cursor.execute("""
                    SELECT 
                        achange
                    FROM 
                        variant
                    WHERE gene = :gene
                """, {"gene": gene})
                #fetch the achanges
                pmkb_variants = self.cursor.fetchall()
                #get pmkb pos:ref_allele:alt_allele
                for line in pmkb_variants:
                    #get pos:ref_allele_alt_allele
                    pmkb_achange = line[0].split(':')
                    pmkb_pos = pmkb_achange[1]
                    pmkb_ref_allele = pmkb_achange[-2]
                    pmkb_alt_allele = pmkb_achange[-1]
                    pmkb_so = pmkb_achange[-3]
                    pmkb_kind = pmkb_achange[0]
                    if input_pos in pmkb_pos:
                        if  pmkb_ref_allele == input_ref_allele:
                            self.cursor.execute(query, {"achange_pmkb": line[0],"gene": gene})
                            qr.append(self.cursor.fetchone())
                        #if positions meet 
                        elif pmkb_ref_allele == '_any' and pmkb_alt_allele == '_any':
                            self.cursor.execute(query, {"achange_pmkb": line[0],"gene": gene})
                            qr.append(self.cursor.fetchone())
                        if pmkb_kind == '_exon':
                            if pmkb_pos == exonno:
                                self.cursor.execute(query,{"achange_pmkb": line[0],"gene": gene})
                                qr.append(self.cursor.fetchone())
                    #handle _any kind
                    if pmkb_so == '_any' or pmkb_so == 'any':
                        if pmkb_pos == '_any':
                            self.cursor.execute(query, {"achange_pmkb": line[0],"gene": gene})
                            qr.append(self.cursor.fetchone())
        #handle insertion mutation cases
        elif variant_type == 'INI' or variant_type == 'inframe_insertion':
            input_ref_alt_pos = re.search(r'(\w{3})(\d+)_?(\w{3})(\d+)ins(\w+)',achange)
            if input_ref_alt_pos is not None:
                ref_alt_pos_catch = input_ref_alt_pos.groups()
                input_start_pos = ref_alt_pos_catch[1]
                input_end_pos = ref_alt_pos_catch[3]
                input_alt_allele = ref_alt_pos_catch[-1]
            
                #Query pmkb database
                self.cursor.execute("""
                    SELECT 
                        achange
                    FROM 
                        variant
                    WHERE 
                        gene = :gene
                """, {'gene': gene})
                pmkb_variants = self.cursor.fetchall()
                #loop through results to get a match
                for line in pmkb_variants:
                    pmkb_achange = line[0].split(':')
                    pmkb_pos = pmkb_achange[1]
                    pmkb_alt_allele = pmkb_achange[-1]
                    input_alt_allele_seq1 = ''
                    pmkb_so = pmkb_achange[-3]
                    pmkb_kind = pmkb_achange[0]
                    
                    #convert alt_allele to one letter amino acid notation
                    if input_alt_allele != None:
                        for i in range(0,len(input_alt_allele)+1,3):
                            if i < len(input_alt_allele):
                                input_alt_allele_seq1 += seq1(alt_allele[i:i+3])
                    if pmkb_pos == input_start_pos:    
                        if pmkb_alt_allele == input_alt_allele_seq1:
                            self.cursor.execute(query, {"achange_pmkb": line[0],"gene": gene})
                            qr.append(self.cursor.fetchone())
                        elif pmkb_ref_allele == '_any' and pmkb_alt_allele == '_any':
                            self.cursor.execute(query, {"achange_pmkb": line[0],"gene": gene})
                            qr.append(self.cursor.fetchone())
                        if pmkb_kind == '_exon':
                            if pmkb_pos == exonno:
                                self.cursor.execute(query,{"achange_pmkb": line[0],"gene": gene})
                                qr.append(self.cursor.fetchone())
                    #handle _any kind
                    if pmkb_so == '_any' or pmkb_so == 'any':
                        if pmkb_pos == '_any':
                            self.cursor.execute(query, {"achange_pmkb": line[0],"gene": gene})
                            qr.append(self.cursor.fetchone())

        elif variant_type == 'CSS' or variant_type == 'complex_substitution':
        
            input_ref_alt_pos = re.search(r'(\w{3})(\d+)_?(\w{3})?(\d+)?delins(\w+)',achange)
            if input_ref_alt_pos is not None:
                ref_alt_pos_catch = input_ref_alt_pos.groups()
                input_start_pos = ref_alt_pos_catch[1]
                input_end_pos = ref_alt_pos_catch[3] if ref_alt_pos_catch[3] != None else '' 
                input_alt_allele = ref_alt_pos_catch[-1]
                self.cursor.execute("""
                    SELECT 
                        achange
                    FROM 
                        variant
                    WHERE 
                        gene = :gene
                    """, {'gene': gene})
                
                pmkb_variants = self.cursor.fetchall()
                for line in pmkb_variants:
                    pmkb_achange = line[0].split(':')
                    #get pos:ref_allele_alt_allele
                    pmkb_start_pos = pmkb_achange[1]
                    #get pos:alt_allele
                    pmkb_end_pos = pmkb_achange[2]
                    pmkb_alt_allele = pmkb_achange[-1]
                    input_alt_allele_seq1 = ''
                    pmkb_so = pmkb_achange[-3]
                    pmkb_kind = pmkb_achange[0]
                    #convert alt_allele to one letter amino acid notation
                    if input_alt_allele != None:
                        for i in range(0,len(input_alt_allele)+1,3):
                            if i < len(input_alt_allele):
                                input_alt_allele_seq1 += seq1(alt_allele[i:i+3])
                    #if variant description has start and end pos - check matching with the pmkb database
                    if input_end_pos != '':
                        if pmkb_start_pos == input_start_pos:
                            if pmkb_end_pos == input_end_pos and pmkb_alt_allele == input_alt_allele_seq1:
                                self.cursor.execute(query, {"achange_pmkb": line[0],"gene": gene})
                                qr.append(self.cursor.fetchall())
                    # if there is no end pos for the input variant ex: p.Asp618delinsGluLeuHis 
                    else:
                        if pmkb_start_pos == input_start_pos and pmkb_alt_allele == input_alt_allele_seq1:
                            self.cursor.execute(query, {"achange_pmkb": line[0],"gene": gene})
                            qr.append(self.cursor.fetchall())
                        elif pmkb_ref_allele == '_any' and pmkb_alt_allele == '_any':
                            self.cursor.execute(query, {"achange_pmkb": line[0],"gene": gene})
                            qr.append(self.cursor.fetchone())
                        if pmkb_kind == '_exon':
                            if pmkb_pos == exonno:
                                self.cursor.execute(query,{"achange_pmkb": line[0],"gene": gene})
                                qr.append(self.cursor.fetchone())
                    #handle _any so
                    if pmkb_so == '_any' or pmkb_so == 'any':
                        if pmkb_pos == '_any':
                            self.cursor.execute(query, {"achange_pmkb": line[0],"gene": gene})
                            qr.append(self.cursor.fetchone())
            #Handle deletion cases
        
        #Handle deletion category
        elif variant_type == 'IND' or variant_type == 'inframe_deletion':
            input_ref_alt_pos = re.search(r'(\w{3})(\d+)_?(\w{3})?(\d+)?del', achange)
            if input_ref_alt_pos is not None:
                ref_alt_pos_catch = input_ref_alt_pos.groups()
                input_ref_allele = seq1(ref_alt_pos_catch[0])
                input_alt_allele = '' if ref_alt_pos_catch[2] == None else seq1(ref_alt_pos_catch[2])
                input_start_pos = ref_alt_pos_catch[1] 
                input_end_pos = ref_alt_pos_catch[3]
                input_ref_alt = input_ref_allele + 'X' + input_alt_allele
                self.cursor.execute("""
                    SELECT 
                        achange
                    FROM 
                        variant
                    WHERE 
                        gene = :gene
                    """, {'gene': gene})
                pmkb_variants = self.cursor.fetchall()
                for line in pmkb_variants:
                    pmkb_achange = line[0].split(':')
                    pmkb_start_pos = pmkb_achange[1]
                    pmkb_ref_allele = pmkb_achange[-2]
                    pmkb_alt_allele = pmkb_achange[-1]
                    pmkb_end_pos = pmkb_achange[2] if type(pmkb_achange[2]) is int else ''
                    pmkb_so = pmkb_achange[-3]
                    pmkb_kind = pmkb_achange[0]
                    #convert input ref_allele to one letter notation
                    #if the deletion is in an interval
                    if input_end_pos != '' and pmkb_end_pos != '' and input_alt_allele != '':
                        if pmkb_start_pos == input_start_pos and pmkb_ref_allele == input_ref_alt and pmkb_end_pos == input_end_pos:
                            self.cursor.execute(query, {"achange_pmkb": line[0],"gene": gene})
                            qr.append(self.cursor.fetchall())
                    #else the deletion is in one position only
                    else:
                        if pmkb_ref_allele == input_ref_allele and pmkb_start_pos == input_start_pos:
                            self.cursor.execute(query, {"achange_pmkb": line[0],"gene": gene})
                            qr.append(self.cursor.fetchone())
                        elif pmkb_ref_allele == '_any' and pmkb_alt_allele == '_any':
                            self.cursor.execute(query, {"achange_pmkb": line[0],"gene": gene})
                            qr.append(self.cursor.fetchone())
                        if pmkb_kind == '_exon':
                            if pmkb_pos == exonno:
                                self.cursor.execute(query,{"achange_pmkb": line[0],"gene": gene})
                                qr.append(self.cursor.fetchone())
                    #handle _any so
                    if pmkb_so == '_any' or pmkb_so == 'any':
                        if pmkb_start_pos == '_any':
                            self.cursor.execute(query, {"achange_pmkb": line[0],"gene": gene})
                            qr.append(self.cursor.fetchone())
        if qr is not None and qr != [] and None not in qr:
            return{
              "gene_name": qr[0][0],
              "tumor_type":qr[0][1],
              "tissue_type":qr[0][2],
              "pmkb_url_interpretation":qr[0][3],
              "interpretations": qr[0][4],
              "citations":qr[0][5],
              "achange": qr[0][6],
              "pmkb_url_variants": qr[0][7]
            }
        _ = secondary_data

def seq1(seq, custom_map = None):
    protein_letters_3to1 = {
    "Ala": "A",
    "Cys": "C",
    "Asp": "D",
    "Glu": "E",
    "Phe": "F",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Lys": "K",
    "Leu": "L",
    "Met": "M",
    "Asn": "N",
    "Pro": "P",
    "Gln": "Q",
    "Arg": "R",
    "Ser": "S",
    "Thr": "T",
    "Val": "V",
    "Trp": "W",
    "Tyr": "Y",
    }
    if custom_map is None:
        custom_map = {"Ter": "*"}
    
    onecode = {k.upper(): v for k,v in protein_letters_3to1.items()}
    onecode.update((k.upper(), v) for k, v in custom_map.items())
    seqlist = [seq[3 * i : 3 * (i+1)] for i in range(len(seq)//3)]
    return "".join(onecode.get(aa.upper()) for aa in seqlist)



