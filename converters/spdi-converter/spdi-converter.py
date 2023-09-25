import re
from oakvar import BaseConverter
from typing import Union, Dict

class Converter(BaseConverter):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Currently var_no is not used in the SPDI format, but it is required by the BaseConverter class
        # Therefore, we initialize it to 0 and increment it for each variant. Maybe another implementation later?
        self.index = 0
        
    def check_format(self, input_path) -> bool:
        """
        This function checks if the input file is in SPDI format and returns a boolean value.
        """
        
        pattern = re.compile(r'^\w+:\d+:[\w-]+:[\w-]+$')

        with open(input_path, 'r') as input_file:
            for line in input_file:
                line = line.strip()
                # Skip empty or comment lines
                if not line or line.startswith("#"):
                    continue
                
                # Return True if the line matches the SPDI pattern, otherwise False
                return bool(pattern.match(line))
            
    def convert_line(self, line: str) -> Union[Dict[str, Union[str, int]], None]:
        """
        This function converts a single line of SPDI format to a dictionary.
        """
        parts = line.strip().split(':')
        if len(parts) < 4:
            return None  
        
        chrom, pos, ref_base, alt_base = parts
        
        var_no = self.index
        
        var_dict = {
            'var_no': var_no,  
            'chrom': chrom,
            'pos': int(pos) + 1,
            'ref_base': ref_base,
            'alt_base': alt_base
        }
        
        self.index += 1  
        
        return [var_dict]
