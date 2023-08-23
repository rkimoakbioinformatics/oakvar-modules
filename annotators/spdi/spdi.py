"""Module for SPDI annotator"""
import oakvar as ov
from oakvar import BaseAnnotator
from scripts.Variant import Variant



class Annotator(BaseAnnotator):
    """
    Annotator class that inherits from BaseAnnotator. It is responsible for calling the Variant class and returning the SPDI notation.
    The value in the SPDI notation then populates the input file as a new column.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.wgs_reader = ov.get_wgs_reader(assembly="hg38")

    def annotate(self, input_data):
        chrom = input_data["chrom"]
        pos = int(input_data["pos"])
        ref_base = input_data["ref_base"]
        alt_base = input_data["alt_base"]

        variant_instance = Variant.classify_variant(chrom, pos, ref_base, alt_base)
        spdi_notation = variant_instance.to_spdi()
        return {"spdi": spdi_notation}
