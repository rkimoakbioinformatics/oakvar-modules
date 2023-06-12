from typing import Optional
from oakvar import BaseAnnotator

class Annotator(BaseAnnotator):
    def annotate(self, input_data: dict, secondary_data: Optional[dict] = None):
        assert input_data is not None
        if not self.cursor:
            return None
        #query interpretations database
        self.cursor.execute(
            """SELECT Gene, 'TUMORType(s)','TissueType(s)','Variant(s)',Tier, PMKBURL, Interpretations, Citations 
            FROM Interpretations""" 
        )
        
        _ = secondary_data
        out = {}
        return out
