import duckdb
from typing import Dict, Any
from oakvar import BaseAnnotator
import json
import re

class Annotator(BaseAnnotator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
            
    def annotate(self, input_data: Dict[str, Any]) -> Dict[str, Any]:    
        missense_variants = self.parse_missense_variants(input_data)
        if not missense_variants:
            return None

        for variant in missense_variants:
            uniprot_id = variant["Uniprot"]
            position = variant["pos"]
            var_aa = variant["var"]
            hit = self.query_score(uniprot_id, position, var_aa)
            if hit not in (None, 0):
                return {"score": hit}

        return None

    def query_score(self, uniprot_id, position, var_aa):
        # List of variant amino acids in the order they are stored in the arrays
        var_aas = ['K', 'R', 'H', 'E', 'D', 'N', 'Q', 'T', 'S', 'C', 'G', 'A', 'V', 'L', 'I', 'M', 'P', 'Y', 'F', 'W']
        
        # Find the index of the requested variant amino acid
        try:
            aa_index = var_aas.index(var_aa)
        except ValueError:
            return None
        
        # SQL query to retrieve the array of scores for the specified position
        query = f'''
        SELECT Scores[{position}][{aa_index + 1}] AS Score
        FROM "{uniprot_id}"
        WHERE Positions[{position}] = {position}
        '''
        
        # Execute the query and fetch the result using self.cursor
        try:
            self.cursor.execute(query)
            result = self.cursor.fetchall()
            if result:
                score = result[0][0] 
                return round(score, 3)
            else:
                return None
        except duckdb.CatalogException:
            return None
        
    def parse_missense_variants(self, input_data):
        result = []

        # Parse the all_mappings JSON string into a dictionary
        all_mappings = json.loads(input_data['all_mappings'])

        # Amino acid three-letter to single-letter code dictionary
        aa_dict = {
            "Ala": "A", "Cys": "C", "Asp": "D", "Glu": "E", "Phe": "F", "Gly": "G", "His": "H",
            "Ile": "I", "Lys": "K", "Leu": "L", "Met": "M", "Asn": "N", "Pro": "P", "Gln": "Q",
            "Arg": "R", "Ser": "S", "Thr": "T", "Val": "V", "Trp": "W", "Tyr": "Y"
        }

        for gene_symbol in all_mappings:
            mappings = all_mappings[gene_symbol]
            for mapping in mappings:
                uniprot_id, aa_change, molecular_consequence, transcript_id, cdna_change, exon_number = mapping

                # Check for "MIS" as the molecular consequence
                if "MIS" in molecular_consequence and aa_change.startswith("p."):
                    # Extract the parts of the amino acid change
                    match = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', aa_change)
                    if match:
                        ref_aa, pos, var_aa = match.groups()
                        if var_aa != "Ter":  # Ensure it is not a termination mutation
                            # Convert the variant amino acid from three-letter to single-letter code
                            var_aa_single = aa_dict.get(var_aa, var_aa)

                            result.append({
                                "Uniprot": uniprot_id,
                                "pos": int(pos),  # Ensure pos is an integer
                                "var": var_aa_single
                            })

        return result