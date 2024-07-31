from fhir.resources.observation import Observation
from fhir.resources.bundle import Bundle
from fhir.resources.observation import Observation
from fhir.resources.core.fhirabstractmodel import *
from Bio import Entrez, SeqIO
import pyhgvs as hgvs
from pyensembl import EnsemblRelease
import ijson
import pathlib
import os

filename = pathlib.Path("C:\Parameters-FindSubjectMolecConseqOutput.json")
f2 = pathlib.Path("C:\Windows\System32\exampleinput__s0.fhir.json")

## Parse the JSON data into a FHIR Observation resource
#try:
#    Bundle.parse_file(filename)
#    print("True")
#except ValueError:
#    print("False")
#
#counter = 0
#
#bundle = Bundle.parse_file(filename)
#for entry in bundle.entry:
#    resource = entry.resource
#
#    if isinstance(resource,Observation):
#        counter += 1
#        print(resource.component)
#
#print(counter)
## Access attributes of the Observation resource
def check_format(file):
    start_terms = ['a','b','c']
    validate = False
    counter = 0
    dict = {}
    with open(file ,'r') as fn: 
        parser = ijson.parse(fn)
        state = None 
        key = None
        for prefix, event, value in parser:
            if event == 'map_key':
                component = value
                if value in start_terms:
                    validate = True
                    return validate
    return validate
check_format(filename)

def convert_file(file):
    start_terms = ['resourceType','name','resource','variant']
    var_dicts = []
    var_dict = {
        "chrom": None ,
        "pos":  None,
        "ref_base": None ,
        "alt_base": None,
        "sample_id": None,
    }
    ref_or_alt = None 
    fhir_dict = {}
    allele_codes = ["69547-8","69551-0"]
    HGVS_codes = ['48004-6','51958-7',]
    counter = 0 
    current_id = None
    with open(file, 'r') as json_file:
        parser = ijson.parse(json_file)
        for prefix,event,value in parser:

            #obtain component name
            if event == 'map_key':
                current_key = value
                if current_key in start_terms:
                        var_dict = {
                            "chrom": None ,
                            "pos":  None,
                            "ref_base": None ,
                            "alt_base": None,
                            "sample_id": None,
                        }

            #look at component value
            elif event =='string' or event == 'number':

                #create dictionary for new id
                #if current_key == 'id' or current_key == 'fullUrl':
                #    current_id = value
                #    var_dict['sample_id'] = current_id
                #    fhir_dict[current_id] = []

                #if current_key in start_terms:
                #    if current_key in empty_dict: 
                #        empty_dict[current_id].append({current_key:value})
                #    else:
                #        empty_dict[current_key] = [value]
                if current_key == 'code':
                    if value == "69547-8":
                        ref_or_alt = 'Ref'
                    if value == "69551-0":
                        ref_or_alt = 'Alt'
                if current_key == 'valueString' and ref_or_alt == "Ref":
                    var_dict['ref_base'] = value
                    #fhir_dict[current_id].append({'ref_allele': value})
                if current_key == 'valueString' and ref_or_alt == 'Alt':
                    var_dict['alt_base'] = value
                    #fhir_dict[current_id].append({'alt_allele' : value})
                    var_dicts.append(var_dict)

            elif event == 'end_map':
                pass
        print(var_dicts)
        return fhir_dict
convert_file(filename)
convert_file(f2)

#Entrez.email = 'akashjrampersad@gmail.com'
#def fetch_genomic_sequence(reference_genome):
#    try:
#        search_term = f'hg38[Human]'
#
#        search_handle = Entrez.esearch(db='nucleotide', term=search_term)
#        search_result = Entrez.read(search_handle)
#
#        genome_id = search_result['IdList'][0]
#
#        fetch_handle = Entrez.efetch(db='nucleotide', id=genome_id, rettype='fasta',retmode='text')
#        genome_sequence = fetch_handle.read()
#        return genome_sequence
#    except Exception as e:
#        print("An error has occured:", e)
#        return None
#
#
##chrom, offset, ref, alt = hgvs.parse_hgvs_name('NM_001395525.1:c.215A>G', fetch_genomic_sequence('bob'),get_transcript=get_transcript)
#
#ensembl = EnsemblRelease(97)
#print(ensembl)
#
#gene_symbol = "MYH7"
#chgvs_notation = "c.2323G>T"
#
#transcript = ensembl.transcripts_by_name(gene_symbol)
#print(transcript)
##transcript = ensembl.transcript_by_gene_symbol(gene_symbol)
##hgvs_parser = pyhgvs.HGVSName(transcript)
##genomic_coordinates = hgvs_parser.c_to_g(chgvs_notation)
#
## Extract information
##chromosome_number = genomic_coordinates.contig
##chromosome_position = genomic_coordinates.pos
##strand = genomic_coordinates.transcript.strand
##
##print("Chromosome Number:", chromosome_number)
##print("Chromosome Position:", chromosome_position)
##print("Strand (+/-):", strand)
##