from typing import Optional
from oakvar import BaseAnnotator

## Additional Packages
import warnings
warnings.simplefilter("ignore") # prevents warning message about JSON schema
from ga4gh.vrs import models
from ga4gh.vrs.dataproxy import SeqRepoRESTDataProxy
from ga4gh.vrs.extras.translator import Translator
from biocommons.seqrepo import SeqRepo

class Annotator(BaseAnnotator):

    def annotate(self, input_data: dict, secondary_data: Optional[dict] = None):
        
        assert input_data is not None
        out = {}
        dp = SeqRepoRESTDataProxy(base_url="https://services.genomicmedlab.org/seqrepo")
        ## Translator
        tlr = Translator(data_proxy=dp,
                         translate_sequence_identifiers=True,  # default
                         normalize=True,                       # default
                         identify=True)                        # default

        chrom = input_data['chrom'].strip('chr')
        pos = input_data['pos'] - 1 # becomes end position in VRS format
        ref = input_data['ref_base']
        alt = input_data['alt_base']
        allele = f'{chrom} : {pos} {ref} > {alt}' # get in format of e.g. "19 : 44908822 C > T"
        vrs_allele = tlr.translate_from(allele, 'beacon') # translate to vrs format
        computed_id = vrs_allele.as_dict()['_id'] # save allele ID

        out['id'] = computed_id
        return out


if __name__ == '__main__':
    annotator = Annotator(sys.argv)
    annotator.run()

