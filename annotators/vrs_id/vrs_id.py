# Copyright (c) 2024 Oak Bioinformatics, LLC
# 
# All rights reserved.
# 
# Do not distribute or use this software without obtaining 
# a license from Oak Bioinformatics, LLC.
# 
# For personal use of non-commercial nature, you may use this software
# after registering with `ov store account create`.
# 
# For research use of non-commercial nature, you may use this software
# after registering with `ov store account create`.
# 
# For use by commercial entities, you must obtain a commercial license
# from Oak Bioinformatics, LLC. Please write to info@oakbioinformatics.com
# to obtain the commercial license.

from typing import Optional
import oakvar as ov
import warnings
warnings.simplefilter("ignore") # prevents warning message about JSON schema
from ga4gh.vrs.dataproxy import SeqRepoRESTDataProxy
from ga4gh.vrs.extras.translator import Translator

class Annotator(ov.BaseAnnotator):

    def annotate(self, input_data: dict, secondary_data: Optional[dict] = None):
        _ = secondary_data
        assert input_data is not None
        out = {}
        dp = SeqRepoRESTDataProxy(base_url="https://services.genomicmedlab.org/seqrepo")
        tlr = Translator(data_proxy=dp,
                         translate_sequence_identifiers=True,  
                         normalize=True,                       
                         identify=True)                        
        chrom = input_data['chrom']
        pos = input_data['pos']
        ref = input_data['ref_base']
        alt = input_data['alt_base']
        reader = ov.get_wgs_reader("hg38")
        if reader is None:
            raise ValueError("Could not get WGS reader")
        prev_base = reader.get_bases(chrom, pos - 1).upper() # get previous base for gnomad format
        if ref == '-':
            ref = prev_base
            alt = prev_base + alt
            pos = pos - 1
        if alt == '-':
            alt = prev_base
            ref = prev_base + ref
            pos = pos - 1
        gnomad_expr = f'{chrom}-{pos}-{ref}-{alt}'
        vrs_allele = tlr.translate_from(gnomad_expr, 'gnomad') # translate gnomad expression to VRS allele
        vrs_id = vrs_allele.as_dict()['_id']
        out['id'] = vrs_id
        return out
