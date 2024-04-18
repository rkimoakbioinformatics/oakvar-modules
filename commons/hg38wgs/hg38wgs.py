from typing import Optional
import twobitreader
import os
from oakvar import BaseCommonModule

class CommonModule (BaseCommonModule):
    def setup (self):
        self.wgs_reader = twobitreader.TwoBitFile(os.path.join(os.path.dirname(__file__), 'data', 'hg38.2bit'))
        self.revbases = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N', 'a':'t', 't':'a', 'g':'c', 'c':'g', 'n':'n'}

    def __getitem__ (self, chrom):
        return self.wgs_reader[chrom]

    def get_bases (self, chrom: str, start: int, end: Optional[int]=None, strand: Optional[str]=None, to_upper: bool=False):
        if end is None:
            end = start
        if chrom not in self.wgs_reader:
            return None
        if strand is None or strand == 1 or strand == '+':
            if start <= end:
                bases = self.wgs_reader[chrom][start - 1:end]
            else:
                bases = ''
                for pos in range(start - 1, end - 2, -1):
                    bases += self.wgs_reader[chrom][pos]
            if to_upper:
                bases = bases.upper()
            return bases
        elif strand == -1 or strand == '-':
            if start <= end:
                bases = ''
                for pos in range(end - 1, start - 2, -1):
                    bases += self.wgs_reader[chrom][pos]
                bases = ''.join([self.revbases[b] for b in bases])
            else:
                bases = ''.join([self.revbases[b] for b in self.wgs_reader[chrom][end - 1:start]])
            if to_upper:
                bases = bases.upper()
            return bases
        else:
            return None

    def slice (self, chrom, start, end=None):
        if end is None or end<start:
            end = start
        elif end == 0:
            raise IndexError(end)
        elif end > 0:
            end = end - 1
        if start <= 0:
            raise IndexError(start)
        else:
            start = start - 1
        return self.wgs_reader[chrom][start:end]
