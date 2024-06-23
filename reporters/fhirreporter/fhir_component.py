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

from fhir.resources.fhirtypes import Uri
from fhir.resources.coding import Coding
from fhir.resources.codeableconcept import CodeableConcept
from fhir.resources.observation import ObservationComponent

def get_base_component_gene():
    code_coding = Coding() # type: ignore
    code_coding.system = Uri("http://loinc.org")
    code_coding.code = "48018-6" # type: ignore
    code_coding.display = "Gene studied [ID]" # type: ignore
    code = CodeableConcept(coding=[code_coding]) # type: ignore
    return ObservationComponent(code=code) # type: ignore

def get_base_component_achange_hgvs():
    code_coding = Coding() # type: ignore
    code_coding.system = Uri("http://loinc.org")
    code_coding.code = "48005-3" # type: ignore
    code_coding.display = "Amino acid change (pHGVS)" # type: ignore
    code = CodeableConcept(coding=[code_coding]) # type: ignore
    return ObservationComponent(code=code) # type: ignore

def get_base_component_cchange_hgvs():
    code_coding = Coding() # type: ignore
    code_coding.system = Uri("http://loinc.org")
    code_coding.code = "48004-6" # type: ignore
    code_coding.display = "DNA change (c.HGVS)" # type: ignore
    code = CodeableConcept(coding=[code_coding]) # type: ignore
    return ObservationComponent(code=code) # type: ignore

def get_base_component_transcript():
    code_coding = Coding() # type: ignore
    code_coding.system = "http://loinc.org" # type: ignore
    code_coding.code = "51958-7" # type: ignore
    code_coding.display = "Transcript reference sequence [ID]" # type: ignore
    code = CodeableConcept(coding=[code_coding]) # type: ignore
    return ObservationComponent(code=code) # type: ignore

def get_base_component_start_stop():
    code_coding = Coding() # type: ignore
    code_coding.system = Uri("http://loinc.org")
    code_coding.code = "81254-5" # type: ignore
    code_coding.display = "Genomic allele start-end" # type: ignore
    code = CodeableConcept(coding=[code_coding]) # type: ignore
    return ObservationComponent(code=code) # type: ignore

def get_base_component_ref_base():
    code_coding = Coding() # type: ignore
    code_coding.system = Uri("http://loinc.org")
    code_coding.code = "69547-8"  # type: ignore # always code for reference allele
    code_coding.display = "Ref nucleotide" # type: ignore
    code = CodeableConcept() # type: ignore
    code.coding = [code_coding] # type: ignore
    return ObservationComponent(code=code) # type: ignore

def get_base_component_alt_base():
    code_coding = Coding() # type: ignore
    code_coding.system = Uri("http://loinc.org")
    code_coding.code = "69551-0" # type: ignore
    code_coding.display = "Alt allele" # type: ignore
    code = CodeableConcept() # type: ignore
    code.coding = [code_coding] # type: ignore
    return ObservationComponent(code=code) # type: ignore

def get_base_component_chrom():
    code_coding = Coding() # type: ignore
    code_coding.system = Uri("http://loinc.org")
    code_coding.code = "48000-4" # type: ignore
    code = CodeableConcept() # type: ignore
    code.coding = [code_coding] # type: ignore
    return ObservationComponent(code=code) # type: ignore

def get_base_component_coord_system():
    code_coding = Coding() # type: ignore
    code_coding.system = Uri("http://loinc.org")
    code_coding.code = "92822-6" # type: ignore
    code = CodeableConcept() # type: ignore
    code.coding = [code_coding] # type: ignore
    return ObservationComponent(code=code) # type: ignore

def get_base_component_feature_consequence():
    code_coding = Coding() # type: ignore
    code_coding.system = "http://hl7.org/fhir/uv/genomics-reporting/STU2/CodeSystem-tbd-codes-cs" # type: ignore
    code_coding.code = "feature-consequence" # type: ignore
    code = CodeableConcept(coding=[code_coding]) # type: ignore
    return ObservationComponent(code=code) # type: ignore

