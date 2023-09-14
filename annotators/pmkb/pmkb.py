from typing import Optional
from oakvar import BaseAnnotator
from Bio.SeqUtils import seq1
class Annotator(BaseAnnotator):
    def annotate(self, input_data: dict, secondary_data: Optional[dict] = None):
        assert input_data is not None
        if not self.cursor:
            return None
        #get data for querying the pmkb database from input_data
        gene = input_data['hugo']
        transcript = input_data['transcript']
        achange = input_data['achange'].split('.')[0]
        exonno = input_data['exonno']
        variant_type = input_date['all_mappings'].split(':')[3]
        #Handle missense variants
        if 'missense' in variant_type:
            input_ref_alt_pos = re.search(r'(\w{3})(\d+)(\w{3})',achange)
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
                    Gene = :gene AND transcript_id = :transcript
            """, {"gene": gene, "transcript": transcript})
            #save results in pmkb_variants variable    
            pmkb_variants = self.cursor.fetchall()
            #loop through results to get a match
            for line in pmkb_variants:
                pmkb_achange = line[0].split(':')
                #for exon cases 
                #             
                #get pos:ref_allele_alt_allele
                pmkb_pos = pmkb_achange[1]
                pmkb_ref_allele = pmkb_achange[-2]
                pmkb_alt_allele = pmkb_achange[-1]
                #match results of pmkb achange with input data
                if pmkb_pos == input_pos and pmkb_ref_allele == input_ref_allele and pmkb_alt_allele == input_alt_allele:
                    self.cursor.execute("""
                        SELECT 
                            gene_name, tumor_type,tissue_type,tier,variant.pmkb_url,interpretations,
                            citations,achange, interpretations_final.pmkb_url
                        FROM 
                            variant
                        LEFT JOIN 
                            interpretations_final
                            ON variant.achange = interpretations_final.variants
                        WHERE 
                            achange = :achange_pmkb
                    """, {"achange_pmkb": pmkb_achange})
                    qr = self.cursor.fetchall()
        #handle frameshift mutations
        elif 'frameshift' in variant_type or 'FSD' in variant_type or 'FSI' in variant_type: 
            #get input data for querying the pmkb database
            input_ref_alt_pos = re.search(r'(\w{3})(\d+)',achange)
            ref_alt_pos_catch = input_ref_alt_pos.groups()
            input_pos = ref_alt_pos_catch[1]
            input_ref_allele = seq1(ref_alt_pos_catch[0])
            #query the pmkb database
            self.cursor.execute("""
                SELECT 
                    achange
                FROM 
                    variant
                WHERE Gene = :gene AND transcript_id = :transcript
            """, {"gene": gene, "transcript": transcript})
            #fetch the achanges
            pmkb_variants = self.cursor.fetchall()
            #get pmkb pos:ref_allele:alt_allele
            for line in pmkb_variants:
                #get pos:ref_allele_alt_allele
                pmkb_pos = achange_pmkb.split(':')[1]
                pmkb_ref_allele = achange_pmkb.split(':')[-2]
                if pmkb_pos == input_pos and pmkb_ref_allele == input_ref_allele:
                    self.cursor.execute("""
                        SELECT 
                            gene_name, tumor_type,tissue_type,tier,variant.pmkb_url,interpretations,
                            citations,achange, interpretations_final.pmkb_url
                        FROM 
                            variant
                        LEFT JOIN 
                            interpretations_final
                            ON variant.achange = interpretations_final.variants
                        WHERE 
                            achange = :achange_pmkb""", {"achange_pmkb": achange_pmkb})
                    qr = self.cursor.fetchall()
        #handle insertion mutation cases
        elif 'insertion' in variant_type or 'inframe_insertion' in variant_type:
            input_ref_alt_pos = re.search(r'(\w{3})(\d+)_?(\w{3})(\d+)ins(\w+)',achange)
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
                    Gene = :gene AND transcript_id = :transcript
            """, {'gene': gene, 'transcript': transcript})
            pmkb_variants = self.cursor.fetchall()
            #loop through results to get a match
            for line in pmkb_variants:
                pmkb_achange = line[0].split(':')
                #for exon cases 
                #
                #get pos:ref_allele_alt_allele
                pmkb_pos = pmkb_achange[1]
                # pmkb_ref_allele = achange_pmkb.split(':')[-2]
                pmkb_alt_allele = pmkb_achange[-1]
                input_alt_allele_seq1 = ''
                #convert alt_allele to one letter amino acid notation
                for i in range(0,len(input_alt_allele)+1,3):
                    if i < len(input_alt_allele):
                        input_alt_allele_seq1 += seq1(alt_allele[i:i+3])
                    
                if pmkb_pos == input_start_pos and pmkb_alt_allele == input_alt_allele_seq1:
                    self.cursor.execute("""
                        SELECT 
                            gene_name, tumor_type,tissue_type,tier,variant.pmkb_url,interpretations,
                            citations,achange, interpretations_fina.pmkb_url
                        FROM 
                            variant
                        LEFT JOIN 
                            interpretations_final
                            ON variant.achange = interpretations_final.variants
                        WHERE 
                            achange = :achange_pmkb
                    """, {"achange_pmkb": pmkb_achange})
                    qr = self.cursor.fetchall()
            
        elif 'indel' in variant_type or 'complex_substitution' in variant_type:
        
            ref_alt_pos_catch = re.search(r'(\w{3})(\d+)_?(\w{3})?(\d+)?delins(\w+)',achange).groups()
            input_start_pos = ref_alt_pos_catch[1]
            input_end_pos = ref_alt_pos_catch[3] if ref_alt_pos_catch[3] != None else '' 
            input_alt_allele = ref_alt_pos_catch[-1]
            self.cursor.execute("""
                SELECT 
                    achange
                FROM 
                    variant
                WHERE 
                    Gene = :gene AND transcript_id = :transcript
                """, {'gene': gene, 'transcript': transcript})
            
            pmkb_variants = self.cursor.fetchall()

            for line in pmkb_variants:
                pmkb_achange = line[0].split(':')
                #for exon cases 
                #
                #get pos:ref_allele_alt_allele
                pmkb_start_pos = pmkb_achange[1]
                #get pos:alt_allele
                pmkb_end_pos = pmkb_achange[2]
                # pmkb_ref_allele = achange_pmkb.split(':')[-2]
                pmkb_alt_allele = pmkb_achange[-1]
                input_alt_allele_seq1 = ''
                #convert alt_allele to one letter amino acid notation
                for i in range(0,len(input_alt_allele)+1,3):
                    if i < len(input_alt_allele):
                        input_alt_allele_seq1 += seq1(alt_allele[i:i+3])
                #if variant description has start and end pos - check matching with the pmkb database
                if input_end_pos != '':
                    if pmkb_start_pos == input_start_pos and pmkb_end_pos == input_end_pos and pmkb_alt_allele == input_alt_allele_seq1:
                        self.cursor.execute("""
                            SELECT 
                                gene_name, tumor_type,tissue_type,tier,variant.pmkb_url,interpretations,
                                citations,achange, interpretations_fina.pmkb_url
                            FROM 
                                variant
                            LEFT JOIN 
                                interpretations_final
                                ON variant.achange = interpretations_final.variants
                            WHERE 
                                achange = :achange_pmkb
                        """, {"achange_pmkb": pmkb_achange})
                        qr = self.cursor.fetchall()
                # if there is no end pos for the input variant ex: p.Asp618delinsGluLeuHis 
                else:
                    if pmkb_start_pos == input_start_pos and pmkb_alt_allele == input_alt_allele_seq1:
                        self.cursor.execute("""
                            SELECT 
                                gene_name, tumor_type,tissue_type,tier,variant.pmkb_url,interpretations,
                                citations,achange, interpretations_fina.pmkb_url
                            FROM 
                                variant
                            LEFT JOIN 
                                interpretations_final
                                ON variant.achange = interpretations_final.variants
                            WHERE 
                                achange = :achange_pmkb
                        """, {"achange_pmkb": pmkb_achange})
                        qr = self.cursor.fetchall()
            #Handle deletion cases
        
        #Handle deletion category
        elif 'deletion' in variant_type:
            ref_alt_pos_catch = re.search(r'(\w{3})(\d+)_?(\w{3})?(\d+)?del', achange).groups()
            input_ref_allele = seq1(ref_alt_pos_catch[0])
            input_alt_allele == seq1(ref_alt_pos_catch[2]) if ref_alt_pos_catch != None else ''
            input_start_pos = ref_alt_pos_catch[1] 
            input_end_pos = ref_alt_pos_catch[3] if ref_alt_pos_catch[3] != None else ''
            input_ref_alt = input_ref_allele + 'X' + input_alt_allele
            self.cursor.execute("""
                SELECT 
                    achange
                FROM 
                    variant
                WHERE 
                    Gene = :gene AND transcript_id = :transcript
                """, {'gene': gene, 'transcript': transcript})
            pmkb_variants = self.cursor.fetchall()

            for line in pmkb_variants:
                pmkb_achange = line[0].split(':')
                pmkb_start_pos = pmkb_achange[1]
                pmkb_ref_allele = pmkb_achange[-2]
                pmkb_end_pos = pmkb_achange[2] if type(pmkb_achang[2]) is int else ''

                #convert input ref_allele to one letter notation
                #if the deletion is in an interval
                if input_end_pos != '' and pmkb_end_pos != '' and input_alt_allele != '':
                    if pmkb_start_pos == input_start_pos and pmkb_ref_allele == input_ref_alt and pmkb_end_pos == input_end_pos:
                        self.cursor.execute("""
                            SELECT 
                                gene_name, tumor_type,tissue_type,tier,variant.pmkb_url,interpretations,
                                citations,achange, interpretations_fina.pmkb_url
                            FROM 
                                variant
                            LEFT JOIN 
                                interpretations_final
                                ON variant.achange = interpretations_final.variants
                            WHERE 
                                achange = :achange_pmkb
                        """, {"achange_pmkb": pmkb_achange})
                        qr = self.cursor.fetchall()
                #else the deletion is in one position only
                else:
                    if pmkb_ref_allele == input_ref_allele and pmkb_start_pos == input_start_pos:
                        self.cursor.execute("""
                            SELECT 
                                gene_name, tumor_type,tissue_type,tier,variant.pmkb_url,interpretations,
                                citations,achange, interpretations_fina.pmkb_url
                            FROM 
                                variant
                            LEFT JOIN 
                                interpretations_final
                                ON variant.achange = interpretations_final.variants
                            WHERE 
                                achange = :achange_pmkb
                        """, {"achange_pmkb": pmkb_achange})
                        qr = self.cursor.fetchall()

        if qr is not None:
            return{
              "gene_name": qr[0],
              "tumor_type":qr[1],
              "tissue_type":qr[2],
              "tier": qr[3],
              "pmkb_url":qr[4],
              "interpretations": qr[5],
              "citations":qr[6],
              "achange": qr[7],
              "pmkb_url": qr[8]
            }
        _ = secondary_data
        # out = {}
        # return out
