from scripts.MnInsertion import MnInsertion
import pytest


class TestMnInsertion:
    """Test class for MnInsertion.py"""
    
    @pytest.mark.parametrize(
    "chrom, position, ref_base, alt_base, is_valid_MnInsertion",
    [
        # Valid MnInsertion scenarios...
        ("chr1", 1000, "-", "GG", True),
        ("chr1", 1000, "-", "GA", True),
        ("chr1", 1000, "-", "GATC", True),
        # Non MnInsertion scenarios...
        ("chr1", 1000, "A", "G", False),
        ("chr1", 1000, "A", "-", False),
        ("chr1", 1000, "AT", "GC", False),
        ("chr1", 1000, "-", "C", False)
    ]
    )
    
    def test_matches_criteria(self, chrom, position, ref_base, alt_base, is_valid_MnInsertion):
        """test to see if MnInsertion.py correctly identifies MnInsertions"""
        assert MnInsertion.matches_criteria(chrom, position, ref_base, alt_base) == is_valid_MnInsertion
        
        
    @pytest.mark.parametrize(
    "initial_position, alt_bases, mock_data, expected_left_context",
    [
        (1000, "GGG", {1000: "A", 999: "G", 998: "G", 997: "G", 996: "A"}, "GGG"),
        (1000, "G", {1000: "A", 999: "A"}, ""),
        (1000, "T", {1000: "T", 999: "T", 998: "T", 997: "T", 996: "A"}, "TTT"),
        (1000, "GTC", {1000: "A", 999: "C", 998: "T", 997: "G", 996: "C", 995: "T", 994: "G", 993: "A"}, "GTCGTC"),
    ]
)
    def test_get_left_context(self, initial_position, alt_bases, mock_data, expected_left_context, mocked_wgs_reader):
    
        def mock_get_bases(chrom, start, end):
            bases = ''
            for pos in range(start, end + 1):
                base = mock_data.get(pos, '')
                bases += base
            return bases

        mocked_wgs_reader.get_bases.side_effect = mock_get_bases

        MnInsertion._wgs_reader = mocked_wgs_reader
        mninsertion = MnInsertion("chr1", initial_position, "-", alt_bases)
        print(f"Result during assertion: {mninsertion.left_context}")
        assert mninsertion.left_context == expected_left_context
        
        
    @pytest.mark.parametrize(
    "initial_position, alt_bases, mock_data, expected_right_context",
    [
        (1000, "GGG", {1000:"G", 1001: "G", 1002: "G", 1003: "A"}, "GGG"),
        (1000, "GA", {1000: "G", 1001: "A", 1002: "A"}, "GA"),
        (1000, "TA", {1000: "T", 1001: "A", 1002: "T", 1003: "A"}, "TATA"),
        (1000, "ACT", {1000: "C", 1001: "T", 1002: "T", 1003: "T", 1004: "A"}, ""),
        (1000, "GATC", {1000: "G", 1001: "A", 1002: "T", 1003: "C", 1004: "G", 1005: "A", 1006: "T", 1007: "C"}, "GATCGATC")
    ]
)
    def test_get_right_context(self, initial_position, alt_bases, mock_data, expected_right_context, mocked_wgs_reader):
        
        def mock_get_bases(chrom, start, end):
            bases = ''
            for pos in range(start, end + 1):
                base = mock_data.get(pos, '')
                bases += base
            return bases

        mocked_wgs_reader.get_bases.side_effect = mock_get_bases

        MnInsertion._wgs_reader = mocked_wgs_reader
        mninsertion = MnInsertion("chr1", initial_position, "-", alt_bases)
        print(f"Result during assertion: {mninsertion.right_context}")
        assert mninsertion.right_context == expected_right_context
        
        
    @pytest.mark.parametrize(
    "initial_position, alt_base, mock_data, expected_ref_base, expected_alt_base",
    [
        # No left or right context
        (1000, "GG", {1000: "G"}, "-", "GG"),
        
        # Only left context
        (1000, "GG", {1000: "G", 999: "G", 998: "G", 997: "G"}, "GG", "GGGG"),
        
        # Only right context
        (1000, "TT", {1000: "T", 1001: "T", 1002: "T", 1003: "T"}, "TTTT", "TTTTTT"),
        
        # Both left and right context
        (1000, "AA", {1000: "A", 999: "A", 998: "A", 1001: "A", 1002: "A", 1003: "G"}, "AAAA", "AAAAAA"),
    ]
)
    def test_construct_contextual_allele(self, initial_position, alt_base, mock_data, expected_ref_base, expected_alt_base, mocked_wgs_reader):
        
        # Mocking method
        def mock_get_bases(chrom, start, end):
            bases = ''
            for pos in range(start, end + 1):
                base = mock_data.get(pos, '')
                bases += base
            return bases
            
        mocked_wgs_reader.get_bases.side_effect = mock_get_bases
        MnInsertion._wgs_reader = mocked_wgs_reader
        
        mninsertion = MnInsertion("chr1", initial_position, "-", alt_base)
        assert mninsertion.ref_base == expected_ref_base
        assert mninsertion.alt_base == expected_alt_base

    @pytest.mark.parametrize(
        "initial_position, alt_base, mock_data, expected_position_after_left_align",
        [
            # Case 1: "AA" repeated, with some on the left
            (1000, "AA", {1000: "A", 999: "A", 998: "A", 997: "A"}, 998),

            # Case 2: "ACG" not repeated on the left
            (1000, "ACG", {1000: "T", 999: "C", 998: "A", 997: "G"}, 1000),

            # Case 3: differing multinucleotides, with varying repetitions on the left
            (1000, "ACG", {1000: "G", 999: "G", 998: "C", 997: "A", 996: "C", 995: "A"}, 997),
            (1000, "TT", {1000: "T", 999: "T", 998: "T", 997: "T", 996: "T"}, 996),
        ]
    )

    
    def test_left_align(self, initial_position, alt_base, mock_data, expected_position_after_left_align, mocked_wgs_reader):
        
        def mock_get_bases(chrom, start, end):
            bases = ''
            for pos in range(start, end + 1):
                base = mock_data.get(pos, '')
                bases += base
            return bases
                
        mocked_wgs_reader.get_bases.side_effect = mock_get_bases
        
        MnInsertion._wgs_reader = mocked_wgs_reader
        mninsertion = MnInsertion("chr1", initial_position, "-", alt_base)
        assert mninsertion.pos == expected_position_after_left_align
        
        
    @pytest.mark.parametrize(
    "chrom, pos, ref_base, alt_base, mock_data, expected_spdi",
    [
        # No left or right context
        ("chr1", 1000, "-", "GAG", {1000: "A", 999: "T"}, "chr1:999:-:GAG"),
        
        # Only left context (assuming 3 G's to the left of the position 1000)
        ("chr1", 1000, "-", "TT", {1000: "A", 999: "T", 998: "T"}, "chr1:997:TT:TTTT"),
        
        # Only right context (assuming 4 T's to the right of the position 1000)
        ("chr2", 1000, "-", "TAT", {1000: "T", 1001: "A", 1002: "T", 1003: "G"}, "chr2:999:TAT:TATTAT"),
        
        # Both left and right context (assuming 2 A's to the left and 3 A's to the right of the position 1000)
        ("chr3", 1000, "-", "ACGT", {999: "T", 998: "G", 997: "C", 996: "A", 995: "T", 994: "G", 993: "C", 992: "A", 1000: "A", 1001: "C", 1002: "G", 1003: "T", 1004: "A", 1005: "C", 1006: "G", 1007: "T"}, "chr3:991:ACGTACGTACGTACGT:ACGTACGTACGTACGTACGT"),        
        # Additional cases can be added here
    ]
)
    def test_to_spdi(self, chrom, pos, ref_base, alt_base, mock_data, expected_spdi, mocked_wgs_reader):

        def mock_get_bases(chrom, start, end):
            bases = ''
            for position in range(start, end + 1):
                base = mock_data.get(position, '')
                bases += base
            return bases

        mocked_wgs_reader.get_bases.side_effect = mock_get_bases
        MnInsertion._wgs_reader = mocked_wgs_reader
        mninsertion = MnInsertion(chrom, pos, ref_base, alt_base)
        assert mninsertion.to_spdi() == expected_spdi