from typing import Tuple 
from typing import List
from typing import Dict
from typing import Iterator
import duckdb
from oakvar import BaseConverter


class Converter(BaseConverter):
    def check_format(self, f) -> bool:
        """
        Detect the format of an input file.
        Arguments:
            f: a file handle to an input file
        Returns:
            bool: True if the input file is for this converter,
                  False if not.
        The example below checks if the input file's first line indicates
        VCF file format.
        """
        if f.name.endswith('.parquet'):
            return True 
        if f.name.endswith('.parquet.gz'):
            return True 
        else: return False
        

    # If your converter module needs something else than
    # the standard way of opening a text input file,
    # read line by line, and coverting each line into
    # a list of dictionaries of variants,
    # you may want to start with modifying
    # convert_file method. In that case, uncomment
    # the below convert_file method and add your implementation.
    #
    def convert_file(self, file, *__args__, exc_handler=None, **__kwargs__) -> Iterator[Tuple[int, List[dict]]]:
        conn = duckdb.connect()
        row_q = 'row_group_num_rows'
        rows = conn.execute(f'SELECT row_group_num_rows FROM parquet_metadata({file})').df().to_dict()
        num_rows = rows[row_q][0]


        start = 0
        chunk = 10000
        end = chunk

        
        while start < num_rows:
            print(start)
            QUERY = f'SELECT * FROM parquet_scan("{file}*") LIMIT {chunk} OFFSET {start}'
            count = conn.execute(QUERY).df().to_dict('index')
            for variant in count:
                count.update({variant+start:count[variant]})
            for line_no in list(count.keys()):
                try:
                    yield line_no, self.convert_line(count[line_no])
                except Exception as e:
                    if exc_handler:
                        exc_handler(line_no,e)
            start+= chunk
            end+= chunk
        return None

    def convert_line(self, line) -> List[Dict]:
        """
        Converts a line from an input file to OakVar's variant dict.
        Arguments:
            l: a string of a line from an input file
        Returns:
            dict: a list of dicts, each dict for a variant collected
                  from the input line. Each dict should have
                  the following required fields:
                  chrom: chromosome name [str]
                  pos: chromosomal position [int]
                  ref_base: reference bases [str]
                  alt_base: altername bases [str]
                  Optional fields for each dict are:
                  sample_id: the ID or name of a sample having the variant [list[str]]
                  tags: a custom tag given to the variant [list[str]]
        """
        current = line
        var_dicts = []
        var_dict = {
        "chrom": "",
        "pos": 0,
        "ref_base": "",
        "alt_base": ""
        }
        var_dict['chrom'] = current['# [1]CHROM']
        var_dict['pos'] = current['[2]POS']
        var_dict["ref_base"] = current['[4]REF']
        var_dict["alt_base"] = current['[5]ALT']
        var_dicts.append(var_dict)
        var_dicts.append(var_dict)
        return var_dicts