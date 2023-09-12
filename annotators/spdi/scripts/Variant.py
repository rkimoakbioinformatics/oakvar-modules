"""The purpose of this module is to provide a base class for all variant types for the SPDI annotation module."""

from abc import ABC, abstractmethod
import oakvar as ov
import importlib


class Variant(ABC):
    """Base class for variants."""

    _wgs_reader = None

    def __init__(self, chrom, pos, ref_base, alt_base):
        self.pos = pos
        self.chrom = chrom
        self.ref_base = ref_base
        self.alt_base = alt_base

    @classmethod
    @abstractmethod
    def matches_criteria(cls, chrom, pos, ref_base, alt_base):
        """Check if the provided data matches the criteria for this variant type."""

    @classmethod
    def classify_variant(cls, chrom, pos, ref_base, alt_base):
        """Classify the input .crv data into a variant type by checking the matches_criteria method of each subclass."""
        print(f"Classifying variant: chrom={chrom}, pos={pos}, ref_base={ref_base}, alt_base={alt_base}") # Print input values
        importlib.import_module('scripts.SNP')
        importlib.import_module('scripts.SnInsertion')
        importlib.import_module('scripts.SnDeletion')
        importlib.import_module('scripts.MNP')
        importlib.import_module('scripts.MNV')
        importlib.import_module('scripts.MnInsertion')
        importlib.import_module('scripts.MnDeletion')
        
        
        for subclass in cls.__subclasses__():
            print(f"Checking subclass: {subclass.__name__}") # Print subclass name
            match = subclass.matches_criteria(chrom, pos, ref_base, alt_base)
            print(f"Result of matches_criteria for {subclass.__name__}: {match}") # Print result of matches_criteria

            if match:
                return subclass(chrom, pos, ref_base, alt_base)

        return UnclassifiedVariant(chrom, pos, ref_base, alt_base)

    @abstractmethod
    def to_spdi(self):
        """Converts the variant to SPDI format."""

    @abstractmethod
    def _left_align(self):
        """Helper method for left alignment."""

    @abstractmethod
    def _get_left_context(self):
        """Abstract method to capture the repeated sequence to the left."""

    @abstractmethod
    def _get_right_context(self):
        """Abstract method to capture the repeated sequence to the right."""

    @abstractmethod
    def _construct_contextual_allele(self):
        """Helper method for constructing contextual allele."""

    @staticmethod
    def initialize_wgs_reader(assembly="hg38"):
        """
        Initialize the WGS reader for the module.
        """
        Variant._wgs_reader = ov.get_wgs_reader(assembly)


class UnclassifiedVariant(Variant):
    """
    The purpose of this class is to populate all the cells of input data that are unclassified.
    """

    def _construct_contextual_allele(self):
        pass

    def _get_left_context(self):
        pass

    def _get_right_context(self):
        pass

    def _left_align(self):
        pass

    @classmethod
    def matches_criteria(cls, chrom, pos, ref_base, alt_base):
        pass

    def to_spdi(self):
        return "Unable to classify the given variant data."