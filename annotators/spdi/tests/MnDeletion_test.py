from scripts.MnDeletion import MnDeletion
import pytest


class TestMnDeletion:
    """Test class for MnDeletion.py"""
    
    @pytest.mark.parametrize(
    "chrom, position, ref_base, alt_base, is_valid_MnDeletion",
    [
        # Valid MnDeletion scenarios...
        ("chr1", 1000, "GATC", "-", True),
        ("chr1", 1000, "GATCGA", "-", True),
        ("chr1", 1000, "GATCGATC", "-", True),
        # Non MnDeletion scenarios...
        ("chr1", 1000, "A", "G", False),
        ("chr1", 1000, "C", "-", False),
        ("chr1", 1000, "AT", "GC", False),
        ("chr1", 1000, "-", "C", False)
    ]
    )
    
    def test_matches_criteria(self, chrom, position, ref_base, alt_base, is_valid_MnDeletion):
        """test to see if MnDeletion.py correctly identifies MnDeletions"""
        assert MnDeletion.matches_criteria(chrom, position, ref_base, alt_base) == is_valid_MnDeletion
        
    @pytest.mark.parametrize(
    "initial_position, ref_bases, mock_data, expected_left_context",
    [
        (1000, "GGG", {1000: "G", 999: "G", 998: "G", 997: "G", 996: "A"}, "GGG"),
        (1000, "GAG", {1000: "G", 999: "G", }, "G"),
        (1000, "T", {1000: "T", 999: "T", 998: "T", 997: "T", 996: "A"}, "TTT"),
        (1000, "GTC", {1000: "G", 999: "C", 998: "T", 997: "G", 996: "C", 995: "T", 994: "G", 993: "A"}, "GTCGTC"),
    ]
)
    def test_get_left_context(self, initial_position, ref_bases, mock_data, expected_left_context, mocked_wgs_reader):

        def mock_get_bases(chrom, start, end, to_upper = True):
            bases = ''
            for pos in range(start, end + 1):
                base = mock_data.get(pos, '')
                bases += base
            return bases

        mocked_wgs_reader.get_bases.side_effect = mock_get_bases

        MnDeletion._wgs_reader = mocked_wgs_reader
        mndeletion = MnDeletion("chr1", initial_position, ref_bases, "-")
        assert mndeletion.left_context == expected_left_context
        
    @pytest.mark.parametrize(
    "initial_position, ref_bases, mock_data, expected_right_context",
    [
        (1000, "GGG", {1003:"G", 1004: "G", 1005: "G", 1006: "A"}, "GGG"),
        (1000, "GA", {1002: "G", 1003: "A", 1004: "A"}, "GA"),
        (1000, "TA", {1000: "T", 1001: "A", 1002: "T", 1003: "A"}, "TA"),
        (1000, "ACT", {1000: "A", 1001: "C", 1002: "T", 1003: "T", 1004: "A"}, ""),
        (1000, "GATC", {1004: "G", 1005: "A", 1006: "T", 1007: "C", 1008: "G", 1009: "A", 1010: "T"}, "GATCGAT"),
    ]
)
    def test_get_right_context(self, initial_position, ref_bases, mock_data, expected_right_context, mocked_wgs_reader):

        def mock_get_bases(chrom, start, end, to_upper = True):
            bases = ''
            for pos in range(start, end + 1):
                base = mock_data.get(pos, '')
                bases += base
            return bases

        mocked_wgs_reader.get_bases.side_effect = mock_get_bases

        MnDeletion._wgs_reader = mocked_wgs_reader
        mndeletion = MnDeletion("chr1", initial_position, ref_bases, "-")
        assert mndeletion.right_context == expected_right_context
        
        
    @pytest.mark.parametrize(
    "initial_position, ref_base, mock_data, expected_ref_base, expected_alt_base",
    [
        # No left or right context
        (1000, "GG", {1000: "G", 1001: "G"}, "GG", "-"),
        
        # Only left context
        (1000, "GTC", {999: "C", 998: "T", 997: "G"}, "GTCGTC", "GTC"),
        
        # Only right context
        (1000, "TAC", {1003: "T", 1004: "A"}, "TACTA", "TA"),
        
        # Both left and right context
        (1000, "AA", {999: "A", 998: "A", 1002: "A", 1003: "A"}, "AAAAAA", "AAAA"),
    ]
)
    def test_construct_contextual_allele(self, initial_position, ref_base, mock_data, expected_ref_base, expected_alt_base, mocked_wgs_reader):

        # Mocking method
        def mock_get_bases(chrom, start, end, to_upper = True):
            bases = ''
            for pos in range(start, end + 1):
                base = mock_data.get(pos, '')
                bases += base
            return bases
            
        mocked_wgs_reader.get_bases.side_effect = mock_get_bases
        MnDeletion._wgs_reader = mocked_wgs_reader
        
        mndeletion = MnDeletion("chr1", initial_position, ref_base, "-")
        assert mndeletion.ref_base == expected_ref_base
        assert mndeletion.alt_base == expected_alt_base
        
        
    @pytest.mark.parametrize(
        "initial_position, ref_base, mock_data, expected_position_after_left_align",
        [
            # Case 1: "AA" repeated, with some on the left
            (1000, "AA", {1000: "A", 999: "A", 998: "A", 997: "A"}, 997),

            # Case 2: "ACG" not repeated on the left
            (1000, "ACG", {1000: "T", 999: "C", 998: "A", 997: "G"}, 1000),

            # Case 3: differing multinucleotides, with varying repetitions on the left
            (1000, "ACG", {1000: "G", 999: "G", 998: "C", 997: "A", 996: "C", 995: "A"}, 997),
            (1000, "TT", {1000: "T", 999: "T", 998: "T", 997: "T", 996: "T"}, 996),
        ]
    )

    
    def test_left_align(self, initial_position, ref_base, mock_data, expected_position_after_left_align, mocked_wgs_reader):
        
        def mock_get_bases(chrom, start, end, to_upper = True):
            bases = ''
            for pos in range(start, end + 1):
                base = mock_data.get(pos, '')
                bases += base
            return bases
                
        mocked_wgs_reader.get_bases.side_effect = mock_get_bases
        
        MnDeletion._wgs_reader = mocked_wgs_reader
        mndeletion = MnDeletion("chr1", initial_position, ref_base, "-")
        assert mndeletion.pos == expected_position_after_left_align
        
    @pytest.mark.parametrize(
    "chrom, pos, ref_base, mock_data, expected_spdi",
    [
        # No left or right context
        ("chr1", 1000, "GAG", {1002: "G", 1001: "A", 1000: "G", 999: "A", 998: "T"}, "chr1:999:GAG:-"),

        # Only left context 
        ("chr1", 1000, "TT", {1000: "T", 999: "T", 998: "T"}, "chr1:997:TTTT:TT"),

        # Only right context
        ("chr2", 1000, "TAT", {1003:"T", 1004: "A", 1005: "T", 1006: "T", 1007: "T"}, "chr2:999:TATTATT:TATT"),

        # Both left and right context 
        ("chr3", 1000, "ACGT", {999: "T", 998: "G", 997: "C", 996: "A",  1004: "A", 1005: "C", 1006: "G", 1007: "T"}, "chr3:995:ACGTACGTACGT:ACGTACGT"),
        # Additional cases can be added here
    ]
)
    def test_to_spdi(self, chrom, pos, ref_base, mock_data, expected_spdi, mocked_wgs_reader):

        def mock_get_bases(chrom, start, end, to_upper = True):
            bases = ''
            for position in range(start, end + 1):
                base = mock_data.get(position, '')
                bases += base
            return bases

        mocked_wgs_reader.get_bases.side_effect = mock_get_bases
        MnDeletion._wgs_reader = mocked_wgs_reader
        mndeletion = MnDeletion(chrom, pos, ref_base, "-")
        assert mndeletion.to_spdi() == expected_spdi