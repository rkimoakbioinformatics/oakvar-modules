from scripts.SnDeletion import SnDeletion
import pytest

class TestSnDeletion:
    
    @pytest.mark.parametrize(
    "chrom, position, ref_base, alt_base, is_valid_SnDeletion, mock_data, expected_spdi",
    [
        # Valid SnDeletion scenarios...
        ("chr1", 1000, "G", "-", True, {999: "A", 1000: "G"}, "chr1:999:G:-"),
        ("chr1", 1000, "G", "-", True, {999: "G", 1000: "G"}, "chr1:999:GG:G"),
        ("chr1", 1000, "G", "-", True, {999: "G", 1000: "G", 1001: "G"}, "chr1:999:GGG:GG"),
        # Non SnDeletion scenarios...
        ("chr1", 1000, "A", "G", False, {999: "A", 1000: "A"}, "chr1:999:A:G"),
        ("chr1", 1000, "-", "A", False, {999: "A", 1000: "A"}, "chr1:999:A:-"),
        ("chr1", 1000, "AT", "GC", False, {999: "A", 1000: "A"}, "chr1:999:AT:GC"),
    ]
    )
    
    def test_matches_criteria(self, chrom, position, ref_base, alt_base, is_valid_SnDeletion, mock_data, expected_spdi):
        assert SnDeletion.matches_criteria(chrom, position, ref_base, alt_base) == is_valid_SnDeletion
        
    @pytest.mark.parametrize(
    "initial_position, ref_base, mock_data, expected_left_context",
    [
        (1000, "G", {1000: "G", 999: "G"}, "G"),
        (1000, "G", {1000: "G", 999: "A"}, ""),
        (1000, "T", {1000: "T", 999: "T", 998: "T", 997: "T", 996: "A"}, "TTT"),
    ]
)
    def test_get_left_context(self, initial_position, ref_base, mock_data, expected_left_context, mocked_wgs_reader):
            
        def mock_get_base(chrom, pos, to_upper = True, data=mock_data):  # Bind mock_data to mock_get_base
            if pos not in data:
                print(f"Position {pos} not found in mock data!")
                print(f"Current mock data: {data}")
            return data.get(pos)
        
        mocked_wgs_reader.get_bases.side_effect = mock_get_base
    
        SnDeletion._wgs_reader = mocked_wgs_reader
        sndeletion = SnDeletion("chr1", initial_position, ref_base, "-")
        #print(f"Result during assertion: {Sndeletion.left_context}")
        assert sndeletion.left_context == expected_left_context
        
    @pytest.mark.parametrize(
        "initial_position, ref_base, mock_data, expected_right_context",
        [
            (1000, "G", {1000: "G", 1001: "G"}, "G"),
            (1000, "G", {1000: "G", 1001: "A"}, ""),
            (1000, "T", {1000: "T", 1001: "T", 1002: "T", 1003: "T", 1004: "A"}, "TTT"),
        ]
    )   
    def test_get_right_context(self, initial_position, ref_base, mock_data, expected_right_context, mocked_wgs_reader):
            
        def mock_get_base(chrom, pos, to_upper = True, data=mock_data):  # Bind mock_data to mock_get_base
            if pos not in data:
                print(f"Position {pos} not found in mock data!")
                print(f"Current mock data: {data}")
            return data.get(pos)
        
        mocked_wgs_reader.get_bases.side_effect = mock_get_base
        
        SnDeletion._wgs_reader = mocked_wgs_reader
        
        sndeletion = SnDeletion("chr1", initial_position, ref_base, "-")
        #print(f"Result during assertion: {Sndeletion.right_context}")
        assert sndeletion.right_context == expected_right_context
        
    @pytest.mark.parametrize(
        "initial_position, ref_base, mock_data, expected_ref_base, expected_alt_base",
        [
        # No left or right context
        (1000, "G", {1000: "G"}, "G", "-" ),
        
        # Only left context
        (1000, "G", {1000: "G", 999: "G", 998: "G", 997: "G"}, "GGGG", "GGG"),
        
        # Only right context
        (1000, "T", {1000: "T", 1001: "T", 1002: "T", 1003: "T"}, "TTTT", "TTT"),
        
        # Both left and right context
        (1000, "A", {1000: "A", 999: "A", 998: "A", 1001: "A", 1002: "A", 1003: "G"}, "AAAAA", "AAAA"),
    ]
        
    )
    
    def test_construct_contextual_allele(self, initial_position, ref_base, mock_data, expected_ref_base, expected_alt_base, mocked_wgs_reader):
        
        # Mocking method
        def mock_get_base(chrom, pos, to_upper = True, data=mock_data):
            return data.get(pos)
        
        mocked_wgs_reader.get_bases.side_effect = mock_get_base
        SnDeletion._wgs_reader = mocked_wgs_reader
        
        sndeletion = SnDeletion("chr1", initial_position, ref_base, "-")
        assert sndeletion.ref_base == expected_ref_base
        assert sndeletion.alt_base == expected_alt_base
        
    @pytest.mark.parametrize(
        "initial_position, mock_data, expected_position_after_left_align",
        
        [
            (1000, {999: "G", 1000: "G", 998: "G", 997: "G"}, 997),
            (1000, {999: "G", 1000: "G", 998: "G", 997: "A"}, 998),
           
        ]
    )
    
    def test_left_align(self, initial_position, mock_data, expected_position_after_left_align, mocked_wgs_reader):

        def mock_get_base(chrom, pos, to_upper = True):
            return mock_data.get(pos, "A")
                
        mocked_wgs_reader.get_bases.side_effect = mock_get_base

        SnDeletion._wgs_reader = mocked_wgs_reader
        sndeletion = SnDeletion("chr1", initial_position, "G", "-")
        assert sndeletion.pos == expected_position_after_left_align
        
    @pytest.mark.parametrize(
    "chrom, pos, ref_base, alt_base, mock_data, expected_spdi",
    [
        # No left or right context
        ("chr1", 1000, "G", "-", {1000: "G", 999: "T"}, "chr1:999:G:-"),
        
        # Only left context (assuming 3 G's to the left of the position 1000)
        ("chr1", 1000, "G", "-", {1000: "G", 999: "G", 998: "G", 997: "G"}, "chr1:996:GGGG:GGG"),
        
        # Only right context (assuming 4 T's to the right of the position 1000)
        ("chr2", 1000, "T", "-", {1000: "T", 1001: "T", 1002: "T", 1003: "T"}, "chr2:999:TTTT:TTT"),
        
        # Both left and right context (assuming 2 A's to the left and 3 A's to the right of the position 1000)
        ("chr3", 1000, "A", "-", {1000: "A", 999: "A", 998: "A", 1001: "A", 1002: "A", 1003: "G"}, "chr3:997:AAAAA:AAAA"),
    ]
)

    def test_to_spdi(self, chrom, pos, ref_base, alt_base, mock_data, expected_spdi, mocked_wgs_reader):
        # Mocking method
        def mock_get_base(chrom, pos, to_upper = True, data=mock_data):
            return data.get(pos)

        # Set up the mock
        mocked_wgs_reader.get_bases.side_effect = mock_get_base
        SnDeletion._wgs_reader = mocked_wgs_reader

        sndeletion = SnDeletion(chrom, pos, ref_base, alt_base)
        assert sndeletion.to_spdi() == expected_spdi      