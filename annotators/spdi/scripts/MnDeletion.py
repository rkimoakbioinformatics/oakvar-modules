from scripts.Variant import Variant

class MnDeletion(Variant):
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
        """Logic to determine if it's an MnDeletion."""
        return len(ref_base) > 1 and alt_base == "-"
    
    def _get_left_context(self):
        left_context = ""
        left_pos = self.pos

        while True:
            # Define the range to check, which is the length of ref_base to the left of the current position
            start_pos = left_pos - len(self.ref_base)
            end_pos = left_pos - 1
            bases_at_pos = self._wgs_reader.get_bases(self.chrom, start_pos, end_pos)

            # Check if the bases match the ref_base; if so, update the left context
            if bases_at_pos == self.ref_base:
                left_context = bases_at_pos + left_context  # Keep track of reference bases being deleted
                left_pos -= len(self.ref_base)
            else:
                # If the bases don't match the ref_base, break the loop
                break

        return left_context
    
    
    def _get_right_context(self):
        right_context = ""
        right_pos = self.pos + len(self.ref_base)  # Start after the reference base

        while True:
            # Define the range to check, which is the length of ref_base to the right of the current position
            start_pos = right_pos
            end_pos = right_pos + len(self.ref_base) - 1
            bases_at_pos = self._wgs_reader.get_bases(self.chrom, start_pos, end_pos)
            
            # Update the right_context if the bases match the ref_base
            if bases_at_pos == self.ref_base:
                right_context += bases_at_pos
                right_pos += len(self.ref_base)
            else:
                break
        print(f"Final right_context: {right_context}") # Debugging print
        return right_context

    def _construct_contextual_allele(self):
        """Construct the reference and alternate alleles with context."""
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

        print(f"Left Context: {self.left_context}")
        print(f"Right Context: {self.right_context}")
        print(f"Ref Base: {self.ref_base}")
        print(f"Alt Base: {self.alt_base}")


    def _left_align(self):
        shift = len(self.left_context)  # Calculate the shift based on the length of the left context
        self.pos -= shift  # Subtract the shift from the current position
        print(f"Shifted position: {self.pos}")  # Print the updated position for debugging

    def to_spdi(self):
        adjusted_pos = self.pos - 1  # Convert from 1-based to 0-based indexing
        return f"{self.chrom}:{adjusted_pos}:{self.ref_base}:{self.alt_base}"