# OakVar Dual License
# 
# Copyright (c) 2023 Oak Bioinformatics, LLC
# 
# This program is dual licensed under the Affero GPL-3.0 or later for 
# non-commercial and open source use, and under a commercial license, 
# which is available for purchase, for closed-source or commercial use.
# 
# For the commercial use, please contact Oak Bioinformatics, LLC 
# for obtaining such a license. OakVar commercial license does not impose 
# the Affero GPL open-source licensing terms, conditions, and limitations. 
# To obtain a commercial-use license of OakVar, please visit our website at
# https://oakbioinformatics.com or contact us at info@oakbioinformatics.com 
# for more information.
# 
# ================
# OpenCRAVAT
# 
# MIT License
# 
# Copyright (c) 2021 KarchinLab
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


from typing import Optional
from oakvar import BaseAnnotator
import sqlite3

class Annotator(BaseAnnotator):

    def annotate(self, input_data: dict, secondary_data: Optional[dict] = None):

        assert input_data is not None
        _ = secondary_data

        chrom = input_data["chrom"].upper()
        pos = input_data["pos"]
        ref = input_data["ref_base"]
        alt = input_data["alt_base"]

        self.cursor.execute(
            f"SELECT AF_ALL, AF_AFR, AF_AMR, AF_ASJ, AF_EAS, AF_SAS, AF_FIN, AF_MID, AF_NFE FROM {chrom} WHERE POS={pos} AND REF='{ref}' AND ALT='{alt}'") 
        qr = self.cursor.fetchone()

        if qr is not None:
            return {
            "AF_ALL": qr[0],
            "AF_AFR": qr[1],
            "AF_AMR": qr[2],
            "AF_ASJ": qr[3],
            "AF_EAS": qr[4],
            "AF_SAS": qr[5],
            "AF_FIN": qr[6],
            "AF_MID": qr[7],
            "AF_NFE": qr[8]
            }

if __name__ == '__main__':
    annotator = Annotator(sys.argv)
    annotator.run()

