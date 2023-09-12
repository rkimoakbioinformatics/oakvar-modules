from scripts.Variant import Variant

class SnDeletion(Variant):
    """
    Class for representing single nucleotide insertions.
    """
    Variant.initialize_wgs_reader()  # set up the _wgs_reader attribute for the class.

    def __init__(self, chrom, pos, ref_base, alt_base):
        super().__init__(chrom, pos, ref_base, alt_base)
        self.left_context = self._get_left_context()
        self.right_context = self._get_right_context()
        self._construct_contextual_allele()
        self._left_align()
        self.to_spdi()

    @classmethod
    def matches_criteria(cls, chrom, pos, ref_base, alt_base):
        # Logic to determine if it's an SnDeletion.
        return alt_base == "-" and len(ref_base) == 1
    
    def _get_left_context(self):
        left_context = ""
        left_pos = self.pos - 1

        while True:
            base_at_pos = self._wgs_reader.get_bases(self.chrom, left_pos)
            # Update the left_context if the base matches the last character of ref_base
            if base_at_pos == self.ref_base[-1]:
                left_context = base_at_pos + left_context
                left_pos -= 1
            else:
                break

        return left_context

    
    def _get_right_context(self):
        right_context = ""
        right_pos = self.pos + 1

        while True:
            base_at_pos = self._wgs_reader.get_bases(self.chrom, right_pos)
            # Update the right_context if the base matches the first character of ref_base
            if base_at_pos == self.ref_base[0]:
                right_context = right_context + base_at_pos
                right_pos += 1
            else:
                break

        return right_context
    
    def _construct_contextual_allele(self):
        "Construct the reference and alternate bases with the contextual bases."
        if not self.left_context and not self.right_context:
            self.ref_base = self.ref_base
            self.alt_base = "-"
        
        elif self.left_context and not self.right_context:
            self.ref_base = self.left_context + self.ref_base
            self.alt_base = self.left_context
            
        elif self.right_context and not self.left_context:
            self.ref_base = self.ref_base + self.right_context
            self.alt_base = self.right_context
               
        else :
            self.ref_base = self.left_context + self.ref_base + self.right_context
            self.alt_base = self.left_context + self.right_context
            
            
        print(f"Left Context: {self.left_context}")
        print(f"Right Context: {self.right_context}")
        print(f"Ref Base: {self.ref_base}")
        print(f"Alt Base: {self.alt_base}")


    def _left_align(self):
        shift = len(self.left_context)  # This determines how many bases we can shift to the left
        if shift > 0:
            self.pos -= shift
        print(len(self.left_context))

    
    def to_spdi(self):
        adjusted_pos = self.pos - 1  # Convert from 1-based to 0-based indexing
        return f"{self.chrom}:{adjusted_pos}:{self.ref_base}:{self.alt_base}"