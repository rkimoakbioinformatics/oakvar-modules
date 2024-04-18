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
from fhir.resources import range, ratio , quantity

def get_base_component_gene():
    code_coding = Coding() # type: ignore
    code_coding.system = Uri("http://loinc.org")
    code_coding.code = "48018-6"   # type: ignore
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

def get_gnomad_allele_frequency():
    code_coding = Coding # type : ignore
    code_coding.system = Uri("http://loinc.org") # type: ignore
    code_coding.code = "92821-8"# type : ignorore
    code_coding.display = "gnomAD v3 Allele Frequency"# type: ignore
    code = CodeableConcept()# type: ignore
    code.coding = [code_coding]# type: ignore
    return ObservationComponent(code=code)# type:ignore

def get_range_comp(range_component):
    comp = range()
    if range_component['high']['value']:
        high_value = range_component['high']['value']
        high_quantity = quantity(value=high_value)
        comp.high= high_quantity
    if range_component['low']['value']:
        low_value = range_component['low']['value']
        low_quantity = quantity(value=low_value)
        comp.low = low_quantity
    return comp

def get_ratio_comp(ratio_component):
    comp = ratio()
    comp.numerator = Quantity(value=ratio_component['numerator']['value'])
    comp.denominator = Quantity(value=ratio_component['denominator']['value'])
    return comp

def get_string_comp(string_component):
    value = string_component['valueString']
    comp = String(value)
    return comp

def get_quantity_comp(quantity_component):
    from fhir.resources.fhirtypes import Code, Decimal
    comp = quantity()
    if quantity_component['code']:
        code = quantity_component['code']
        comp.code= Code(code)
    if quantity_component['system']:
        system = quantity_compnent['system']
        comp.system = Uri(system)
    if quantity_component['value']:
        value = quantity_component['value']
        comp.value = Decimal(value) 
    if quantity_component['comparator']:
        comparator = quantity_component['comparator']
        comp.comparator = Code(comparator)
    if quantity_component['unit']:
        unit = String(value= quantity_component['unit'])
        comp.unit = String(unit)
    return comp

def get_reference_comp(reference_component):
    from fhir.resources import reference
    from fhir.resources.fhirtypes import String, IdentifierType

    comp = reference()
    if reference_component['reference']:
        comp.reference = String(reference_component['reference'])
    if reference_component['type']:
        comp.type = Uri(reference_component['type'])
    if reference_component['identifier']:
        comp.identifier = IdentifierType(reference_component['identifier'])
    if reference_component['display']:
        comp.display = String(reference_component['display'])
    return comp 






def get_codeable_concept_generic(codeable_concept_component):
    comp = None
    system = codeable_concept_component['coding']['system']
    code = codeable_concept_component['coding']['code']
    if codeable_concept_component['coding']['display']:
        display = codeable_concept_component['coding']['display']
        comp = fhir_coding.get_codeable_concept(system, code, display)
    else:
        comp = fhir_coding.get_codeable_concept(system, code)
    return comp

def get_boolean_value_comp(boolean_component):
    from fhir.resources.fhirtypes import Boolean
    value = boolean_component['valueBoolean']
    comp = Boolean(value)
    return comp 

def get_effDateTime_comp(time_component):
    from fhir.resources.fhirtypes import datetime
    comp = datetime(time_component['effectiveDateTime'])
    return comp
def get_integer_comp(integer_comp):
    from fhir.resources.fhirtypes import Integer
    comp = Integer(integer_comp['valueInteger'])
    return comp

def get_period_comp(period_comp):
    from fhir.resources import period
    from fhir.resources.fhirtypes import DateTime
    comp = period()
    if period_comp['start']:
        comp.start = DateTime(period_comp['start'])
    if period_comp['end']:
        comp.end = DateTime(period_comp['end'])
    return comp

def get_attachment_comp(attachment_comp):
    from fhir.resources import attachment
    comp = attachment()
    if attachment_comp['contentType']:
        from fhir.resources.fhirtypes import Code
        comp.contentType = Code(attachment_comp['contentType'])
    if attachment_comp['language']:
        from fhir.resources.fhirtypes import Code
        comp.language = Code(attachment_comp['language'])
    if attachment_comp['data']:
        from fhir.resources.fhirtypes import Base64Binary 
        comp.data = Base64Binary(attachment_comp['data'])
    if attachment_comp['url']:
        from fhir.resources.fhirtypes import Url
        comp.url = Url(attachment_comp['url'])
    if attachment_comp['size']:
        from fhir.resources.fhirtypes import Integer64
        comp.size = Integer64(attachment_comp['size'])
    if attachment_comp['hash']:
        from fhir.resources.fhirtypes import Base64Binary
        comp.hash = Base64Binary(attachment_comp['hash'])
    if attachment_comp['title']:
        comp.title = String(attachment_comp['title'])
    if attachment_comp['creation']:
        from fhir.resources.fhirtypes import datetime
        comp.creation = datetime(attachment_comp['creation'])
    if attachment_comp['height']:
        from fhir.resources.fhirtypes import PositiveInt
        comp.height = PositiveInt(attachment_comp['height'])
    if attachment_comp['width']:
        from fhir.resources.fhirtypes import PositiveInt
        comp.width = PositiveInt(attachment_comp['width'])
    if attachment_comp['frames']:
        from fhir.resources.fhirtypes import PositiveInt
        comp.frames = PositiveInt(attachment_comp['frames'])
    if attachment_comp['duration']:
        from fhir.resources.fhirtypes import Decimal
        comp.duration = Decimal(attachment_comp['duration'])
    if attachment_comp['pages']:
        from fhir.resources.fhirtypes import PositiveInt
        comp.pages = PositiveInt(attachment_comp['pages'])
    return comp

def get_SampledData_comp(sampled_data_component):
    from fhir.resources import sampleddata
    from fhir.resources.fhirtypes import Code
    from fhir.resources.fhirtypes import PositiveInt

    comp = sampleddata.SampledData()
    comp.origin = quantity(sampled_data_component['origin'])
    comp.intervalUnit = Code(sampled_data_component['intervalUnit'])
    comp.dimensions = PositiveInt(sampled_data_component['dimensions'])

    if sampled_data_component['interval']:
        from fhir.resources.fhirtypes import Decimal
        comp.interval = Decimal(sampled_data_component['interval'])

    if sampled_data_component['factor']:
        from fhir.resources.fhirtypes import Decimal 
        comp.factor = Decimal(sampled_data_component['decimal'])

    if sampled_data_component['lowerLimit']:
        from fhir.resources.fhirtypes import Decimal
        comp.lowerLimit = Decimal(sampled_data_component['lowerLimit'])
    if sampled_data_component['upperLimit']:
        from fhir.resources.fhirtypes import Decimal
        comp.upperLimit = Decimal(sampled_data_component['upperLimit'])
    if sampled_data_component['codeMap']:
        from fhir.resources.fhirtypes import Canonical
        comp.codeMap = Canonical(sampled_data_component['codeMap'])
    if sampled_data_component['offsets']:
        from fhir.resources.fhirtypes import String
        comp.offsets = String(sampled_data_component['offsets'])
    if sampled_data_component['data']:
        from fhir.resources.fhirtypes import String
        comp.data = String(sampled_data_component['data'])
    
    return comp
    










    
