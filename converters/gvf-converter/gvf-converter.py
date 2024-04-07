from typing import List
from typing import Dict
from oakvar import BaseConverter



class Converter(BaseConverter):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.format_name = 'gvf'

    def check_format(self, f) -> bool:
        from pathlib import Path

        if not Path(f).exists():
            return False

        if f.endswith(('gff', 'gvf', 'gff3')):
            # reset format_name to the value of the given file
            self.format_name = f.split('.')[-1]
            return True
        else:
            return False

    def convert_line(self, line) -> List[Dict]:
        """
        Converts a line from the file into a structured dictionary.
        """
        var_dicts = []

        if line.startswith('#') or not line.strip():
            return self.IGNORE

        delimiter = "\t" if self.format_name in ("gff", "gff3") else " "
        line_values = line.split(delimiter)
        if len(line_values) < 8:  # Basic validation of line structure
            return None
        
        chrom_val = line_values[0] if 'chr' in line_values[0] else 'chr' + line_values[0]
        pos_val = line_values[3]
        end_pos = line_values[4]
        ref_base, alt_base, sample = '', '', 'unknown'  # Initialize with default values


        if self.format_name == 'gvf':
            value_attrs = line_values[-1].split(';')
            for val in value_attrs:
                if 'Reference_seq' in val:
                    ref_base = val.split('=')[-1]
                elif 'Variant_seq' in val:
                    alt_base = val.split('=')[-1]

            if "," in alt_base:
                alt_bases = []
                alt_bases = alt_base.split(',')
                for alt in alt_bases:
                    var_dict = {
                        "chrom": chrom_val,
                        "pos": pos_val,
                        'end_pos': end_pos,
                        "ref_base": ref_base,
                        "alt_base": alt,
                        "sample_id": sample,
                    }
                    var_dicts.append(var_dict)
            else:
                var_dict = {
                    "chrom": chrom_val,
                    "pos": pos_val,
                    'end_pos': end_pos,
                    "ref_base": ref_base,
                    "alt_base": alt_base,
                    "sample_id": sample,
                }
                var_dicts.append(var_dict)
            return var_dicts
        elif self.format_name in ("gff", "gff3"):
            var_dict = {
                "chrom": chrom_val,
                "pos": pos_val,
                'end_pos': end_pos,
                "ref_base": "!",
                "alt_base": "!",
                "sample_id": sample,
            }
            var_dicts.append(var_dict)
            return var_dicts