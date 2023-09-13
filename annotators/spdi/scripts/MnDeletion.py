from scripts.Variant import Variant

class MnDeletion(Variant):
    """
    Class for representing single nucleotide deletions.
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
        """Logic to determine if it's an MnDeletion."""
        return len(ref_base) > 1 and alt_base == "-"
    
    def _get_left_context(self):
        left_context = ""
        left_pos = self.pos - 1  # Start at the immediate left base of the position.
        ref_index = len(self.ref_base) - 1  # Start with the last base of ref_base.

        while True:
            base_at_pos = self._wgs_reader.get_bases(self.chrom, left_pos, left_pos, to_upper = True)
            if base_at_pos == self.ref_base[ref_index]:
                left_context = base_at_pos + left_context
                ref_index -= 1
                left_pos -= 1
                # If we have checked all bases of ref_base, restart from the last base of ref_base
                if ref_index < 0:
                    ref_index = len(self.ref_base) - 1
            else:
                break

        return left_context

    def _get_right_context(self):
        right_context = ""
        right_pos = self.pos + len(self.ref_base)  # Start immediately after the deletion
        ref_index = 0  # Start with the first base of ref_base.

        while True:
            base_at_pos = self._wgs_reader.get_bases(self.chrom, right_pos, right_pos, to_upper = True)
            if base_at_pos == self.ref_base[ref_index]:
                right_context += base_at_pos
                ref_index += 1
                right_pos += 1
                # If we have checked all bases of ref_base, restart from the first base of ref_base
                if ref_index >= len(self.ref_base):
                    ref_index = 0
            else:
                break

        return right_context

    def _construct_contextual_allele(self):
        if not self.left_context and not self.right_context:
            self.ref_base = self.ref_base
            self.alt_base = "-"
        elif self.left_context and not self.right_context:
            self.ref_base = self.left_context + self.ref_base
            self.alt_base = self.left_context
        elif not self.left_context and self.right_context:
            self.ref_base = self.ref_base + self.right_context
            self.alt_base = self.right_context
        else:  # Both contexts exist
            self.ref_base = self.left_context + self.ref_base + self.right_context
            self.alt_base = self.left_context + self.right_context
            
    def _left_align(self):
        shift = len(self.left_context)  # Calculate the shift based on the length of the left context
        self.pos -= shift  # Subtract the shift from the current position

    def to_spdi(self):
        adjusted_pos = self.pos - 1  # Convert from 1-based to 0-based indexing
        return f"{self.chrom}:{adjusted_pos}:{self.ref_base}:{self.alt_base}"