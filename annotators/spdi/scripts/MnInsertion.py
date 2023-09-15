from scripts.Variant import Variant

class MnInsertion(Variant):
    """
    Class for processing multi nucleotide insertions.
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
        """Logic to determine if it's an MnInsertion."""
        return ref_base == "-" and len(alt_base) > 1
    
    def _get_left_context(self):
        left_context = ""
        left_pos = self.pos - 1  # Start at the immediate left base of the position.
        alt_index = len(self.alt_base) - 1  # Start with the last base of alt_base.

        while True:
            base_at_pos = self._wgs_reader.get_bases(self.chrom, left_pos, left_pos, to_upper = True)
            if base_at_pos == self.alt_base[alt_index]:
                left_context = base_at_pos + left_context
                alt_index -= 1
                left_pos -= 1
                # If we have checked all bases of alt_base, restart from the last base of alt_base
                if alt_index < 0:
                    alt_index = len(self.alt_base) - 1
            else:
                break

        return left_context

    def _get_right_context(self):
        right_context = ""
        right_pos = self.pos  # Start at the immediate right base of the position.
        alt_index = 0  # Start with the first base of alt_base.

        while True:
            base_at_pos = self._wgs_reader.get_bases(self.chrom, right_pos, right_pos, to_upper = True)
            if base_at_pos == self.alt_base[alt_index]:
                right_context += base_at_pos
                alt_index += 1
                right_pos += 1
                # If we have checked all bases of alt_base, restart from the first base of alt_base
                if alt_index >= len(self.alt_base):
                    alt_index = 0
            else:
                break

        return right_context

    def _construct_contextual_allele(self):
        """Construct the reference and alternate alleles with context."""
        if not self.left_context and not self.right_context:
            self.ref_base = "-"
            self.alt_base = self.alt_base
        elif self.left_context and not self.right_context:
            self.ref_base = self.left_context
            self.alt_base = self.left_context + self.alt_base
        elif not self.left_context and self.right_context:
            self.ref_base = self.right_context
            self.alt_base = self.alt_base + self.right_context
        else:  # Both contexts exist
            self.ref_base = self.left_context + self.right_context
            self.alt_base = self.left_context + self.alt_base + self.right_context

    def _left_align(self):
        shift = len(self.left_context)  # Calculate the shift based on the length of the left context
        self.pos -= shift  # Subtract the shift from the current position

    def to_spdi(self):
        adjusted_pos = self.pos - 1  # Convert from 1-based to 0-based indexing
        return f"{self.chrom}:{adjusted_pos}:{self.ref_base}:{self.alt_base}"