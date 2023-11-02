from scripts.Variant import Variant

class SnInsertion(Variant):
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
        # Logic to determine if it's an SnInsertion.
        return ref_base == '-' and len(alt_base) == 1

    def _get_left_context(self):
        left_context = ''
        left_pos = self.pos - 1

        while True:
            base_at_pos = self._wgs_reader.get_bases(self.chrom, left_pos, to_upper=True)
            # Update the left_context if the base matches the last character of alt_base
            if base_at_pos == self.alt_base[-1]:
                left_context = base_at_pos + left_context
                left_pos -= 1
            else:
                break

        return left_context

    def _get_right_context(self):
        """Retrieve the repeated sequence to the right of the insertion."""
        right_context = ''
        right_pos = self.pos + 1  # Start from the position next to the insertion

        base_at_insertion_pos = self._wgs_reader.get_bases(self.chrom, self.pos, to_upper=True)
        if base_at_insertion_pos != self.alt_base[-1]:  # If the base at insertion position doesn't match the inserted base
            return right_context  # Return empty context

        while True:
            base_at_pos = self._wgs_reader.get_bases(self.chrom, right_pos, to_upper=True)
            if base_at_pos == self.alt_base[-1]:
                # If the base at the insertion position matches and it's the first time in the loop
                if not right_context:
                    right_context += base_at_insertion_pos
                right_context += base_at_pos
                right_pos += 1
            else:
                break
            
        return right_context
    
    def _construct_contextual_allele(self):
        """Construct the reference and alternate alleles with context."""
        if not self.left_context and not self.right_context:
            self.ref_base = '-'
            self.alt_base = self.alt_base[-1]
        elif self.left_context and not self.right_context:
            self.ref_base = self.left_context
            self.alt_base = self.left_context + self.alt_base[-1]
        elif not self.left_context and self.right_context:
            self.ref_base = self.right_context
            self.alt_base = self.alt_base[-1] + self.right_context
        else:  # Both contexts exist
            self.ref_base = self.left_context + self.right_context
            self.alt_base = self.left_context + self.alt_base[-1] + self.right_context
            
    def _left_align(self):
        shift = len(self.left_context)
        if shift > 0:
            self.pos -= shift

    def to_spdi(self):
        adjusted_pos = self.pos - 1  # Convert from 1-based to 0-based indexing
        return f'{self.chrom}:{adjusted_pos}:{self.ref_base}:{self.alt_base}'