from typing import Optional
from oakvar import BaseAnnotator

class Annotator(BaseAnnotator):
    def annotate(self, input_data: dict, secondary_data: Optional[dict] = None):
        assert input_data is not None
        if not self.cursor:
            return None
        gene_name = input_data["hugo"]
        a_change = input_data["achange"]
        #query interpretations database
        self.cursor.execute(
            """SELECT 
	                gene_name, tumor_type, tissue_type, tier, interpretations.pmkb_url,
	                interpretations, citations, achange, variants.pmkb_url
                FROM
	                interpretations
                LEFT JOIN
	                variants as variants ON
	                interpretations.gene_name = variants.gene
                WHERE
	                (CASE
		                WHEN variants.achange IS NULL THEN interpretations.gene_name = :gene
		                ELSE interpretations.gene_name = :gene AND variants.achange = :achange
		            END)
            """, {"gene": gene_name, "achange":a_change}
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
              "achange": qr[7],
              "pmkb_url": qr[8]
            }
        _ = secondary_data
        # out = {}
        # return out
