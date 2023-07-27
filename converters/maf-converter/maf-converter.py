from typing import List
from typing import Dict
from oakvar import BaseConverter


# TODO the official documentations has 12 steps for validating a file. Maybe implement in the future?


class Converter(BaseConverter):
    MAF_STANDARD_COLS = (
        "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome", "Start_Position", "End_Position",
        "Strand", "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1",
        "Tumor_Seq_Allele2", "dbSNP_RS", "dbSNP_Val_Status", "Tumor_Sample_Barcode",
        "Matched_Norm_Sample_Barcode", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
        "Tumor_Validation_Allele1", "Tumor_Validation_Allele2", "Match_Norm_Validation_Allele1",
        "Match_Norm_Validation_Allele2", "Verification_Status", "Validation_Status", "Mutation_Status",
        "Sequencing_Phase", "Sequence_Source", "Validation_Method",
        "Score", "BAM_File", "Sequencer", "Tumor_Sample_UUID"
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.format_name = "maf"

    def check_format(self, f) -> bool:
        from pathlib import Path

        if not Path(f).exists():
            return False

        file = open(f)
        line = file.readline()
        file_cols = line.split('\t')

        # if any of the standard columns are not in the file headers, return False
        for col in self.MAF_STANDARD_COLS:
            if col not in file_cols:
                return False

        return True

    def convert_line(self, line) -> List[Dict]:
        line_list = line.split('\t')
        var_dicts = []

        # inferred from column 5 "Chromosome"
        maf_chrom = line_list[4]

        # inferred from column 6 "Start_Position"
        maf_pos = line_list[5]

        # inferred from column 10 "Reference_Allele"
        maf_ref = line_list[11]

        # inferred from column 11 "Tumor_Seq_Allele1"
        maf_alt = line_list[10]

        # inferred from column 15 "Tumor_Sample_Barcode"
        maf_sample = line_list[15]

        var_dict = {
            'chrom': maf_chrom,
            'pos': maf_pos,
            'ref_base': maf_ref,
            'alt_base': maf_alt,
            'sample_id': maf_sample,
            "var_no": '',
        }

        var_dicts.append(var_dict)

        return var_dicts
