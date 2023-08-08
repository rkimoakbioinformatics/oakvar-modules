import oakvar as ov
from oakvar import BaseAnnotator

class Variant:
    """Base class for variants.
    
    Attributes:
        chrom (str): The chromosome where the variant is located.
        pos (int): The position of the variant on the chromosome.
        ref_base (str): The reference base at the variant position.
        alt_base (str): The alternate base at the variant position.
    """
    def __init__(self, chrom, pos, ref_base, alt_base):
        
        self.pos = pos
        self.chrom = chrom
        self.ref_base = ref_base
        self.alt_base = alt_base

    def to_spdi(self):
        """Converts the variant to SPDI format. Subclases must implement this method.
        
        Returns:
            str: The SPDI representation of the variant.
        """
        raise NotImplementedError("Subclasses must implement this method.")

class SNP(Variant):
    """
    A class representing a single nucleotide polymorphism (SNP).
    """
    def __init__(self, chrom, pos, ref_base, alt_base):
        super().__init__(chrom, pos, ref_base, alt_base)
        self.pos = pos - 1  # Adjusting for 0-based indexing 

    def to_spdi(self):
        spdi_notation = f"{self.chrom}:{self.pos}:{self.ref_base}:{self.alt_base}"
        return spdi_notation

class SnInsertion(Variant):
    """
    Class for representing single nucleotide insertions.
    """

    def __init__(self, chrom, pos, ref_base, alt_base, wgs_reader):
        super().__init__(chrom, pos, ref_base, alt_base)
        self.pos = pos - 1  # Adjusting for 0-based indexing
        self.wgs_reader = wgs_reader
        self.left_align()

    def left_align(self):
        """
        Adjusts the insertion position to the leftmost possible position.
        """
        while self.pos > 0:
            base_before = self.wgs_reader.get_bases(self.chrom, self.pos - 1)
            if base_before == self.alt_base:
                self.pos -= 1
            else:
                break

    def to_spdi(self):
        spdi_notation = f"{self.chrom}:{self.pos}:-:{self.alt_base}"
        return spdi_notation

class SnDeletion(Variant):
    
    """Class for single nucleotide deletions."""
    
    def __init__(self, chrom, pos, ref_base, alt_base, wgs_reader):
        super().__init__(chrom, pos, ref_base, alt_base)
        self.pos = pos - 1  # Adjusting for 0-based indexing
        self.wgs_reader = wgs_reader
        self.left_align()
        
    def left_align(self):
        while self.pos > 0:
            base_before = self.wgs_reader.get_bases(self.chrom, self.pos - 1)
            if base_before == self.ref_base:
                self.pos -= 1
            else:
                break
        
    def to_spdi(self):
        spdi_notation = f"{self.chrom}:{self.pos}:{self.ref_base}:-"
        return spdi_notation
        
class MNP(Variant):
    """Class for multi-nucleotide polymorphisms (MNPs).
    """
    def __init__(self, chrom, pos, ref_base, alt_base, wgs_reader):
        assert len(ref_base) == len(alt_base), "Ref_base and alt_base must be the same length for MNPs"
        super().__init__(chrom, pos, ref_base, alt_base)
        self.pos = pos - 1  # Adjusting for 0-based indexing
        self.wgs_reader = wgs_reader
        self.left_align()

    def left_align(self):
        """Left-aligns the MNP in case it is part of a repeat region."""
        while self.pos > 0:
            base_before = self.wgs_reader.get_bases(self.chrom, self.pos - len(self.ref_base))
            if base_before == self.ref_base:
                self.pos -= 1
            else:
                break

    def to_spdi(self):
        spdi_notation = f"{self.chrom}:{self.pos}:{self.ref_base}:{self.alt_base}"
        return spdi_notation

class MnInsertion(Variant):
    """Class for multi-nucleotide insertions."""
    def __init__(self, chrom, pos, ref_base, alt_base, wgs_reader):
        super().__init__(chrom, pos, ref_base, alt_base)
        self.pos = pos - 1  # Adjusting for 0-based indexing
        self.wgs_reader = wgs_reader
        self.left_align()
        
    def left_align(self):
        while self.pos > 0:
            current_bases = self.wgs_reader.get_bases(self.chrom, self.pos + 1, self.pos + len(self.alt_base))
            if current_bases == self.alt_base:
                self.pos -= 1
            else:
                break

    def to_spdi(self):
        spdi_notation = f"{self.chrom}:{self.pos}:-:{self.alt_base}"
        return spdi_notation

class MnDeletion(Variant):
    """Class for multi-nucleotide deletions."""
    def __init__(self, chrom, pos, ref_base, alt_base, wgs_reader):
        super().__init__(chrom, pos, ref_base, alt_base)
        self.pos = pos - 1  # Adjusting for 0-based indexing
        self.wgs_reader = wgs_reader
        self.left_align()
        
    def left_align(self):
        while self.pos > 0:
            current_bases = self.wgs_reader.get_bases(self.chrom, self.pos + 1, self.pos + len(self.ref_base))
            if current_bases == self.ref_base:
                self.pos -= 1
            else:
                break

    def to_spdi(self):
        spdi_notation = f"{self.chrom}:{self.pos}:{self.ref_base}:-"
        return spdi_notation
    
    
class MNV(Variant):
    """Class for multi-nucleotide variations."""
    def __init__(self, chrom, pos, ref_base, alt_base, wgs_reader):
        super().__init__(chrom, pos, ref_base, alt_base)
        self.pos = pos - 1  # Adjusting for 0-based indexing
        self.wgs_reader = wgs_reader
        self.left_align()

    def left_align(self):
        while self.pos > 0:
            bases_before = self.wgs_reader.get_bases(self.chrom, self.pos - len(self.ref_base), self.pos - 1)
            if bases_before == self.ref_base:
                self.pos -= 1
            else:
                break

    def to_spdi(self):
        spdi_notation = f"{self.chrom}:{self.pos}:{self.ref_base}:{self.alt_base}"
        return spdi_notation


class Annotator(BaseAnnotator):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)  # Pass any additional arguments to the superclass - Error with kwargs
        self.wgs_reader = ov.get_wgs_reader(assembly="hg38")
        
    def convert_variant_to_spdi(self, chrom, pos, ref_base, alt_base):

        #SNP Variant Type
        if len(ref_base) == 1 and len(alt_base) == 1 and ref_base != "-" and alt_base != "-":
            variant = SNP(chrom, int(pos), ref_base, alt_base)
            
        #MNP Variant Type
        elif len(ref_base) > 1 and len(alt_base) > 1 and len(ref_base) == len(alt_base):
            variant = MNP(chrom, int(pos), ref_base, alt_base, self.wgs_reader)
            
        #MNV Variant Type
        elif len(ref_base) > 1 and len(alt_base) > 1 and len(ref_base) != len(alt_base):
            variant = MNV(chrom, int(pos), ref_base, alt_base, self.wgs_reader)
            
        #SN Insertion Variant Type
        elif ref_base == "-" and len(alt_base) == 1 and alt_base != "-":
            variant = SnInsertion(chrom, int(pos), ref_base, alt_base, self.wgs_reader)
            
        #MN Insertion Variant Type
        elif ref_base == "-" and len(alt_base) > 1:
            variant = MnInsertion(chrom, int(pos), ref_base, alt_base, self.wgs_reader)
            
        #SN Deletion Variant Type
        elif len(ref_base) == 1 and alt_base == "-" and ref_base != "-":
            variant = SnDeletion(chrom, int(pos), ref_base, alt_base, self.wgs_reader)
            
        #MN Deletion Variant Type
        elif len(ref_base) > 1 and alt_base == "-":
            variant = MnDeletion(chrom, int(pos), ref_base, alt_base, self.wgs_reader)
            
        return variant.to_spdi()
         
    def annotate(self, input_data):
        chrom = input_data["chrom"]
        pos = int(input_data["pos"])
        ref_base = input_data["ref_base"]
        alt_base = input_data["alt_base"]
        spdi_notation = self.convert_variant_to_spdi(chrom, pos, ref_base, alt_base)
        return {"spdi": spdi_notation}