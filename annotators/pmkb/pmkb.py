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
        #Handle missense variants
        if 'missense' in input_data['all_mappings'].split(':')[3]:
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
                pmkb_achange = line[0]
                #for exon cases 
                #             
                #get pos:ref_allele_alt_allele
                pmkb_pos = achange_pmkb.split(':')[1]
                pmkb_ref_allele = achange_pmkb.split(':')[-2]
                pmkb_alt_allele = achange_pmkb.split(':')[-1]
                #match results of pmkb achange with input data
                if pmkb_pos == input_pos and pmkb_ref_allele == input_ref_allele and pmkb_alt_allele == input_alt_allele:
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
                    """, {"achange_pmkb": achange_pmkb})
                    qr = self.cursor.fetchall()
        #handle frameshift mutations
        elif 'frameshift' in input_date['all_mappings'].split(':')[3]:
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
                    varian
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
                            citations,achange, interpretations_fina.pmkb_url
                        FROM 
                            variant
                        LEFT JOIN 
                            interpretations_final
                            ON variant.achange = interpretations_final.variants
                        WHERE 
                            achange = :achange_pmkb""", {"achange_pmkb": achange_pmkb})
                    qr = self.cursor.fetchall()
        elif 'insertion' in input_date['all_mappings'].split(':')[3]:
            pass    
        elif 'indel' in input_data['all_mappings'].split(':')[3]:
            pass
        elif 'deletion' in input_data['all_mappings'].split(':')[3]:
            pass
        
        #handle frameshift mutation

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
