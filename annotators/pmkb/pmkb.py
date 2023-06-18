from typing import Optional
from oakvar import BaseAnnotator

class Annotator(BaseAnnotator):
    def annotate(self, input_data: dict, secondary_data: Optional[dict] = None):
        assert input_data is not None
        if not self.cursor:
            return None
        gene_name = input_data["Description"]
        #query interpretations database
        self.cursor.execute(
            """SELECT 
                interpretations.Gene, "TumorType(s)", "TissueType(s)", "Variant(s)", Tier,
                interpretations.PMKBURL, Interpretations, Citations, Description, variants.PMKBURL
            FROM	
	            "PMKB_Interpretations_Complete_20230526-1011" as interpretations
            LEFT JOIN
                "PMKB_Variants_Complete_20230526-1011" as variants ON
                interpretations.Gene = variants.Gene
            WHERE
                interpretations.Gene = ?
            """, (gene_name,)
        )
        qr = self.cursor.fetchone()
        if qr is not None:
            return{
              "gene_name": qr[0],  
              "tumor_type":qr[1],
              "tissue_type":qr[2],
              "tier": qr[3],
              "pmkb_url":qr[4],
              "interpretations": qr[5],
              "citations":qr[6],
              "gene_description": qr[7],
              "pmkb_url": qr[8]
            }
        _ = secondary_data
        # out = {}
        # return out
