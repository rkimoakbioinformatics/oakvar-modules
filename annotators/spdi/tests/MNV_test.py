from scripts.MNV import MNV

class TestMNV:
    
    def test_spdi(self):
        # Initialize the MNV object with correct MNV arguments
        mnv = MNV("chr1", 1000, "ATC", "GC")
        
        # Call the to_spdi method
        result = mnv.to_spdi()
        
        assert result == "chr1:999:ATC:GC"
    
    def test_matches_criteria(self):
        # Test for valid MNV criteria
        assert MNV.matches_criteria("chr1", 1000, "ATC", "GC")
        
        # Test for invalid MNV criteria due to ref_base being "-"
        assert not MNV.matches_criteria("chr1", 1000, "-", "GC")
        
        # Test for invalid MNV criteria due to alt_base being "-"
        assert not MNV.matches_criteria("chr1", 1000, "AT", "-")
        
        # Test for invalid MNV criteria due to ref_base and alt_base being "-"
        assert not MNV.matches_criteria("chr1", 1000, "-", "-")
        
        # Test for invalid MNV criteria due to ref_base and alt_base being the same length
        assert not MNV.matches_criteria("chr1", 1000, "AT", "GC")