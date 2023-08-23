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
        ref_base = ''

        if line.startswith('#'):
            return self.IGNORE

        line_values = line.split()

        # chrom is first value
        chrom_val = line_values[0]

        if 'chr' not in chrom_val:
            chrom_val = 'chr' + chrom_val

        # pos
        pos_val = line_values[3]

        # ref_base
        if self.format_name == 'gvf':
            value_attrs = line_values[-1].split(';')

            for val in value_attrs:
                if 'Reference_seq' in val:
                    ref_base = val.split('-')[-1]
                else:
                    ref_base = '-'

        var_dicts = []
        var_dict = {
            "chrom": chrom_val,
            "pos": pos_val,
            "ref_base": ref_base,
            "alt_base": "T",
            "sample_id": "sample1",
            'var_no': self.line_no
        }
        var_dicts.append(var_dict)
        print(var_dict)
        return var_dicts
