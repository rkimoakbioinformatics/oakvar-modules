from scripts.Variant import Variant

class SNP(Variant):
    
    """
    A class representing a single nucleotide polymorphism (SNP).
    """
    
    def __init__(self, chrom, pos, ref_base, alt_base):
        super().__init__(chrom, pos, ref_base, alt_base)
        self.pos = pos - 1  # Adjusting for 0-based indexing
         

    @classmethod
    def matches_criteria(cls, chrom, pos, ref_base, alt_base):
        return len(ref_base) == 1 and len(alt_base) == 1 and ref_base != "-" and alt_base != "-"

    def to_spdi(self):
        spdi_notation = f"{self.chrom}:{self.pos}:{self.ref_base}:{self.alt_base}"
        return spdi_notation
    
    def _left_align(self):
        raise NotImplementedError("Left alignment not applicable for SNP")

    def _get_left_context(self):
        raise NotImplementedError("Getting left context not applicable for SNP")

    def _get_right_context(self):
        raise NotImplementedError("Getting right context not applicable for SNP")

    def _construct_contextual_allele(self):
        raise NotImplementedError("Constructing contextual allele not applicable for SNP")