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

from typing import List
from fhir_consts import *


def     get_codeable_concept(system: str, code: str, display: str):
    from fhir.resources.codeableconcept import CodeableConcept

    codeable_concept = CodeableConcept() # type: ignore
    coding = get_coding(system, code, display)
    codeable_concept.coding = [coding] # type: ignore
    return codeable_concept

def get_coding_generic(coding_json):
    coding = Coding()
    coding.code= coding_json['code']
    coding.system = coding_json['system']
    if coding_json['display']:
        coding.display = coding_json['display']
    return coding

def get_coding(system, code, display):
    from fhir.resources.fhirtypes import Uri
    from fhir.resources.coding import Coding

    coding = Coding() # type: ignore
    if system:
        coding.system = Uri(system)
    if code:
        coding.code = code
    if display:
        coding.display = display
    return coding
    
def get_codeable_concept_ensembl_transcript(ensembl: str):
    return get_codeable_concept("http://www.ensembl.org", ensembl, ensembl)

def get_codeable_concept_refseq_transcript(refseq: str):
    return get_codeable_concept("http://www.ncbi.nlm.nih.gov/refseq", refseq, refseq)

def get_codeable_concept_so(sos: List[str]):
    from fhir.resources.codeableconcept import CodeableConcept

    codings = []
    for so in sos:
        so_id = so_term_to_id_dict.get(so)
        if so_id:
            coding = get_coding("http://sequenceontology.org", so, so_id)
            codings.append(coding)
    return CodeableConcept( # type: ignore
        coding=codings
    )

def get_codeable_concept_hgvs(entity: str, change: str):
    c = f"{entity}:{change}"
    return get_codeable_concept("http://varnomen.hgvs.org", c, c)

def get_observation_code():
    return get_codeable_concept("http://loinc.org", "69548-6", "")

def get_codeable_concept_chrom(chrom: str, chrom_la: str):
    return get_codeable_concept("http://loinc.org", chrom_la, f"Chromosome {chrom}")

def get_observation_category_code():
    return get_codeable_concept(
        "http://terminology.hl7.org/CodeSystem/observation-category",
        "laboratory",
        ""
    )

