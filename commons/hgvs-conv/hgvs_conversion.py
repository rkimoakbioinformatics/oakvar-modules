from Bio import Entrez, SeqIO
from oakvar import BaseCommonModule
import hgvs.parser

class CommonModule(BaseCommonModule):
    var_example = "NM_001366781.1:c.90T>C"

    def obtain_chromosome(self,transcript_id):
        Entrez.email = "akashjrampersad@gmail.com"

        handle = Entrez.efetch(db='nucleotide', id=transcript_id, rettype="gb", retmode='text')
        record = SeqIO.read(handle,'genbank')

        for feature in record.features:
            if "chromosome" in feature.qualifiers:
                chromosome = feature.qualifiers["chromosome"][0]
                return chromosome

        return None

    def hgvs_parser(self, hgvs_code):

        hp = hgvs.parser.Parser()
        var_dict = {"transcript": None,
                    "pos": None,
                    "ref" : None,
                    "alt" : None,
                    'chrom' : None
                    }

        var_g = hp.parse_hgvs_variant(hgvs_code)
        var_dict["transcript"] = var_g.ac
        var_dict["pos"] = var_g.posedit.pos.start.base
        var_dict["ref"] = var_g.posedit.edit.ref
        var_dict["alt"] = var_g.posedit.edit.alt
        var_dict["chrom"] = self.obtain_chromosome(var_g.ac)
        return var_dict
