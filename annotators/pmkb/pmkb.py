from typing import Optional
from oakvar import BaseAnnotator

class Annotator(BaseAnnotator):
    def annotate(self, input_data: dict, secondary_data: Optional[dict] = None):
        assert input_data is not None
        if not self.cursor:
            return None
        gene_name = input_data["Gene"]
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
              "Gene": qr[0],  
              "TumorType(s)":qr[1],
              "TissueType(s)":qr[2],
              "Variant(s)": qr[3],
              "Tier":qr[4],
              "PMKBURL": qr[5],
              "Interpretations":qr[6],
              "Citations": qr[7],
              "Description": qr[8],
              "PMKBURL": qr[9]
            }
        _ = secondary_data
        # out = {}
        # return out
