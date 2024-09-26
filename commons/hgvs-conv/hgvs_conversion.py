from Bio import Entrez, SeqIO
import hgvs
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.variantmapper
import hgvs.assemblymapper



def variant_converter(hgvs_variant):

    hp = hgvs.parser.Parser()
    components = hgvs_variant.split(":")
    variant_info = components[1]
    if variant_info[0] == "c":
        hdp = hgvs.dataproviders.uta.connect()
        var_c = hp.parse_hgvs_variant(hgvs_variant)
        am37 = hgvs.assemblymapper.AssemblyMapper(hdp,assembly_name='GRCh37',alt_aln_method = 'blat')
        var_g = am37.c_to_g(var_c)
        return var_g
    if variant_info[0] == "g":
        var_g = hp.parse_hgvs_variant(hgvs_variant)
        return var_g


def obtain_chromosome(transcript_id):
    Entrez.email = "akashjrampersad@gmail.com"

    handle = Entrez.efetch(db='nucleotide', id=transcript_id, rettype="gb", retmode='text')
    record = SeqIO.read(handle,'genbank')

    for feature in record.features:
        if "chromosome" in feature.qualifiers:
            chromosome = feature.qualifiers["chromosome"][0]
            return chromosome

    return None

def hgvs_parser(hgvs_code):
    var_dict = {"transcript": None,
                'chrom' : None,
                "pos": None,
                "ref" : None,
                "alt" : None
                }
    variant = variant_converter(hgvs_code)

    var_dict["transcript"] = variant.ac
    var_dict["pos"] = variant.posedit.pos.start.base
    var_dict["ref"] = variant.posedit.edit.ref
    var_dict["alt"] = variant.posedit.edit.alt
    var_dict["chrom"] = obtain_chromosome(variant.ac)
    return var_dict









