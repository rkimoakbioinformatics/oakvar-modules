from scripts.SNP import SNP

class TestSNP:
    
    def test_spdi(self):
        # Initialize the SNP object with correct SNP arguments
        snp = SNP("chr1", 1000, "A", "G")
        
        # Call the to_spdi method
        result = snp.to_spdi()
        
        # Assert that the result is as expected using pytest's assert
        assert result == "chr1:999:A:G"
    
    def test_matches_criteria(self):
        # Test for valid SNP criteria
        assert SNP.matches_criteria("chr1", 1000, "A", "G")
        
        # Test for invalid SNP criteria due to ref_base being "-"
        assert not SNP.matches_criteria("chr1", 1000, "-", "G")
        
        # Test for invalid SNP criteria due to alt_base being "-"
        assert not SNP.matches_criteria("chr1", 1000, "A", "-")
        
        # Test for invalid SNP criteria due to both ref_base and alt_base being multi-nucleotide
        assert not SNP.matches_criteria("chr1", 1000, "AT", "GC")