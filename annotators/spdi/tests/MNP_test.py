from scripts.MNP import MNP

class TestMNP:
    
    def test_spdi(self):
        # Initialize the MNP object with correct MNP arguments
        mnp = MNP("chr1", 1000, "AT", "GC")
        
        # Call the to_spdi method
        result = mnp.to_spdi()
        
        # Assert that the result is as expected using pytest's assert
        assert result == "chr1:999:AT:GC"  # Adjust according to your SPDI logic
    
    def test_matches_criteria(self):
        # Test for valid MNP criteria
        assert MNP.matches_criteria("chr1", 1000, "AT", "GC")
        
        # Test for invalid MNP criteria due to ref_base being "-"
        assert not MNP.matches_criteria("chr1", 1000, "-", "GC")
        
        # Test for invalid MNP criteria due to alt_base being "-"
        assert not MNP.matches_criteria("chr1", 1000, "AT", "-")