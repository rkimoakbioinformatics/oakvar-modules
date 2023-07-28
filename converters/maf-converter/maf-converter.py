import re
from typing import List
from typing import Dict
from oakvar import BaseConverter


# TODO the official documentations has 12 steps for validating a file. Maybe implement in the future?


class Converter(BaseConverter):
    MAF_STANDARD_COLS = (
        "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele1"
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

        # if any of the standard columns are not in the file headers, return False
        for col in self.MAF_STANDARD_COLS:
            if not re.search(col, line, flags=re.IGNORECASE):
                return False

        return True

    def convert_line(self, line) -> List[Dict]:
        line_list = line.split('\t')
        var_dicts = []

        if line.startswith("Hugo_Symbol"):
            return self.IGNORE

        # inferred from column 5 "Chromosome"
        maf_chrom = line_list[4]

        # inferred from column 6 "Start_Position"
        maf_pos = line_list[5]

        # inferred from column 10 "Reference_Allele"
        maf_ref = line_list[10]

        # inferred from column 11 "Tumor_Seq_Allele1"
        maf_alt = line_list[11]

        # inferred from column 15 "Tumor_Sample_Barcode"
        maf_sample = line_list[15]

        var_dict = {
            'chrom': maf_chrom,
            'pos': maf_pos,
            'ref_base': maf_ref,
            'alt_base': maf_alt,
            'sample_id': maf_sample,
            'var_no': self.line_no
        }

        var_dicts.append(var_dict)

        return var_dicts
