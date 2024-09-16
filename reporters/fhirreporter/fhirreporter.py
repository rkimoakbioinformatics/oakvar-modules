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
from typing import Dict
from typing import Tuple
import uuid
import hashlib
import yaml
import sqlite3
from oakvar import BaseReporter
from fhir.resources.patient import Patient
from fhir.resources.humanname import HumanName
from fhir.resources.codeableconcept import CodeableConcept
from fhir.resources.coding import Coding
from fhir.resources.reference import Reference
from fhir.resources.bundle import Bundle, BundleEntry
from fhir.resources import attachment, sampleddata
from fhir.resources.fhirtypes import Uri
from fhir.resources.fhirtypes import Boolean, String , Integer
from fhir.resources import range, ratio, period
from fhir.resources.fhirtypes import DateTime, Time 
from fhir.resources.fhirtypes import MetaType
from fhir.resources.observation import Observation
from fhir.resources.identifier import Identifier
from fhir.resources.quantity import Quantity
import os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from fhir_component import *
from fhir_coding import *
from fhir_consts import *





class Reporter(BaseReporter):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.base_id: str = ""
        self.levels_to_write: List[str] = self.conf.get("pages", ["variant"])
        self.filenames: List[str] = []
        self.counter: int = 0
        self.samples: List[str] = []
        self.num_rows: int = 0
        self.dict_entries: Dict[str, List[BundleEntry]] = {}
        self.dict_bundles: Dict[str, Bundle] = {}
        self.dict_patient: Dict[str, Reference] = {}
        self.ensembl_to_refseq_d: Dict[str, List[str]] = {}
        self.allele_frequency = None
        self.yaml_terms = {
            'code': fhir_coding.get_coding_generic,
            "valueAttachment": fhir_component.get_attachment_comp,
            "valueBoolean" : fhir_component.get_boolean_value,
            "valueCodeableConcept": fhir_component.get_codeable_concept_generic,
            "effectiveDateTime": fhir_component.get_effDateTime_comp,
            "valueInteger": fhir_component.get_integer_comp,
            "valuePeriod": fhir_component.get_period_comp,
            "valueQuantity": fhir_component.get_quantity_comp,
            "valueRange": fhir_component.get_range_comp,
            "valueRatio": fhir_component.get_ratio_comp,
            "valueReference": fhir_component.get_reference_comp ,
            "valueSampledData": fhir_component.get_SampledData_comp,
            "valueString": fhir_component.get_string_comp
            }
        self.yaml_components = []

    def setup(self):
        from oakvar import get_mapper
        from oakvar.lib.system import get_user_conf

        self.make_samples()
        # establish filename with fhir suffix
        self.prefix = self.savepath
        for sample in self.samples:
            # create sample id
            sample_id = f"{self.base_id} + {sample}"
            hex_sample = hashlib.md5(sample_id.encode("utf-8")).hexdigest()
            # create patient
            sample_2_patient = Patient() # type: ignore
            name = HumanName() # type: ignore
            name.use = "official" # type: ignore
            name.given = [sample] # type: ignore
            sample_2_patient.name = [name] # type: ignore
            patient_id = str(uuid.UUID(hex=hex_sample))
            # add to patient dict for reference in row observations
            subject = Reference(type="Patient") # type: ignore
            subject.resource_type = "Reference"
            subject.reference = f"urn:uuid:{patient_id}" # type: ignore
            self.dict_patient[sample] = subject
            patient_entry = BundleEntry( # type: ignore
                resource=sample_2_patient, fullUrl=f"urn:uuid:{patient_id}"
            )
            self.dict_entries[sample] = [patient_entry]
            # create sample bundle and and it to dictionary
            name = f"bundle + {sample} + self.str_id"
            hex_sample = hashlib.md5(name.encode("utf-8")).hexdigest()
            bundle = Bundle(
                type="collection",
                identifier=Identifier(value=str(uuid.UUID(hex=hex_sample))), # type: ignore
            )
            self.dict_bundles[sample] = bundle
        user_conf = get_user_conf()
        mapper_name = user_conf.get("genemapper", "gencodoe")
        mapper = get_mapper(mapper_name)
        if mapper:
            self.ensembl_to_refseq_d = mapper.get_ensembl_to_refseq_dict()

    def write_table_row(self, row):
        variant_id = self.get_variant_id(row)
        variant_uri = Uri(f"urn:uuid:{variant_id}")
        molecular_consequence_id = self.get_mcid(row)
        molecular_consequence_uri = Uri(f"urn:uuid:{molecular_consequence_id}")
        # get samples that have variant(row)
        samples = row["tagsampler__samples"].split(",")
        for sample in samples:
            entry = self.dict_entries[sample]
            observation = self.get_observation(sample, row)
            molecular_consequence = self.get_molecular_consequence(variant_id, sample, row)
            variant_bundle = BundleEntry(resource=observation, fullUrl=variant_uri) # type: ignore
            entry.append(variant_bundle)
            molecular_consequence_bundle = BundleEntry(resource=molecular_consequence, fullUrl=molecular_consequence_uri) # type: ignore
            entry.append(molecular_consequence_bundle)
            if self.module_options.get("all_transcripts") == "true":
                self.add_all_mappings(row, variant_id, sample, entry)

    def end(self):
        if self.module_options.get("all_transcripts") == "true":
            for sample in self.samples:
                self.dict_bundles[sample].entry = self.dict_entries[sample] # type: ignore
                obs = self.dict_bundles[sample]
                json_str = obs.json(indent=2)
                filename = str(self.prefix) + f"__{sample}.all.fhir.json"
                with open(filename, "w", encoding="utf-8") as wf:
                    wf.write(str(json_str))
                self.filenames.append(filename)
        else:
            for sample in self.samples:
                self.dict_bundles[sample].entry = self.dict_entries[sample] # type: ignore
                obs = self.dict_bundles[sample]
                json_str = obs.json(indent=2)
                filename = str(self.prefix) + f"__{sample}.primary.fhir.json"
                with open(filename, "w", encoding="utf-8") as wf:
                    wf.write(str(json_str))
                self.filenames.append(filename)
        return self.filenames

    def make_samples(self):
        conn = sqlite3.connect(self.dbpath)
        curs = conn.cursor()
        curs.execute("SELECT DISTINCT base__sample_id FROM sample")
        for sample in curs.fetchall():
            if sample[0].count(",") < 1:
                self.samples.append(sample[0])
        conn.close()

    def make_base_id(self):
        # get sample names
        conn = sqlite3.connect(self.dbpath)
        curs = conn.cursor()
        # get number of rows
        curs = conn.cursor()
        curs.execute("SELECT COUNT(*) from variant")
        self.num_rows = curs.fetchone()[0]
        # get str for id generation
        curs = conn.cursor()
        curs.execute('select colval from info where colkey="input_paths"')
        # get input_path and split is so that only path is part of id.
        self.base_id = curs.fetchone()[0].split(" ", 1)[-1]
        curs = conn.cursor()
        curs.execute('select colval from info where colkey="annotators"')
        self.base_id += curs.fetchone()[0][1:-1]
        self.base_id = self.base_id[1:-1]
        self.base_id = self.base_id[-32:]
        curs = conn.cursor()
        curs.execute('select colval from info where colkey="mapper"')
        self.base_id += curs.fetchone()[0]
        # make this bundle UUID

    def create_dict(self, keys):
        unique_dict = {}
        for key in keys:
            if key not in unique_dict:
                unique_dict[key] = []
        return unique_dict

    def uuid_maker(self, val: str):
        hex_str = hashlib.md5(val.encode("utf-8")).hexdigest()
        return uuid.UUID(hex=hex_str)

    def should_write_level(self, level):
        if self.levels_to_write is None:
            return True
        elif level in self.levels_to_write:
            return True
        else:
            return False

    def write_preface(self, level: str):
        if level not in self.levels_to_write:
            return

    def add_fhir_ref_base(self, row: dict, components: list):
        component = get_base_component_ref_base()
        component.valueString = row["base__ref_base"]
        components.append(component)

    def add_fhir_alt_base(self, row: dict, components: list):
        component = get_base_component_alt_base()
        component.valueString = row["base__alt_base"]
        components.append(component)

    def add_fhir_chrom(self, row: dict, components: list):
        component = get_base_component_chrom()
        chrom = row["base__chrom"]
        chrom_la = chrom_dict.get(chrom)
        chrom_short = chrom[3:] if chrom.startswith("chr") else chrom
        if not chrom:
            return
        component.valueCodeableConcept = \
            get_codeable_concept_chrom(chrom_la, chrom_short) # type: ignore
        components.append(component)

    def add_fhir_coord_system(self, row: dict, components: list):
        _ = row
        _ = components
        component = get_base_component_coord_system()
        component_code = CodeableConcept(
            coding=[
                (
                    Coding( # type: ignore
                        system=Uri("http://loinc.org"),
                        code="LA30102-0",
                        display=String("1-based character counting"),
                    )
                )
            ]
        )
        component.valueCodeableConcept = component_code # type: ignore
        components.append(component)

    def add_fhir_pos(self, row: dict, components: list):
        pos = row["base__pos"]
        pos_end = row["base__pos_end"]
        component = get_base_component_start_stop()
        st_sne_value_low = Quantity(value=pos) # type: ignore
        st_sne_value_high = Quantity(value=pos_end) # type: ignore
        component.valueRange = RangeType(
            low=st_sne_value_low, high=st_sne_value_high
        )
        components.append(component)

    def add_fhir_mc_ensembl_transcript(self, ensembl: str, mapping_comps: list):
        component = get_base_component_transcript()
        component.valueCodeableConcept = \
            get_codeable_concept_ensembl_transcript(ensembl) # type: ignore
        mapping_comps.append(component)

    def add_fhir_components(self, components: list):
        list_components = components['components']
        for component in list_components:
            code = component['code']
            coding = code['coding']
            component_code = coding['code']
            
    def add_gnomad_AF(self, frequency, components: list):
        component = get_gnomad_allele_frequency()
        frequency = valueQuantity(value=frequency)
        frequency.system = Uri('http://unitsofmeasure.org')
        frequency.code = String("1")
        component.valueQuantity = frequency
        components.append(component)


    def write_generic_comp_code(comp_yaml: dict):
        data = yaml.safe_load(comp_yaml)
        components = data['module_options']['fhirreporter']['components']
        list_components = [component['component'] for component in components]
        for component in list_components:
            keys = list(component.keys())
            counter = 0
            while counter <= len(keys-1):
                value_key = keys[counter]
                value_json = component[value_key]
                component_to_add = self.value_terms[value_key](value_json)
                self.yaml_components.append(component_to_add)

    
        










    def get_observation(self, sample, row):
        # Creates CC for Laboratory.
        components = []
        self.add_fhir_chrom(row, components)
        self.add_fhir_coord_system(row, components)
        self.add_fhir_pos(row, components)
        self.add_fhir_ref_base(row, components)
        self.add_fhir_alt_base(row, components)
        observation = Observation( # type: ignore
            status="final",
            code=get_observation_code(),
            subject=self.dict_patient[sample],
            category=[get_observation_category_code()],
        )
        observation.component = components
        observation.meta = MetaType()
        observation.meta.profile = [ # type: ignore
            "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/variant"
        ]
        return observation

    def get_variant_id(self, row):
        id_maker_v = (
            (self.base_id)
            + f"""{str(row['base__chrom']) 
                 + str(row['base__pos']) 
                 + str(row['base__ref_base']) 
                 + str(row['base__alt_base'])}"""
        )
        variant_id = self.uuid_maker(id_maker_v + self.base_id)
        return variant_id

    def get_mcid(self, row):
        id_maker_mc = (
            (self.base_id)
            + f"""{str(row['base__chrom']) 
                 + str(row['base__pos']) 
                 + str(row['base__ref_base']) 
                 + str(row['base__alt_base'])
                 + "molecular consequence 1"}"""
        )
        mc_id = self.uuid_maker(id_maker_mc + self.base_id)
        return mc_id

    def get_mc_row_init(self, variant_id, sample):
        mc_row = Observation(
            status="final",
            category=[
                CodeableConcept(
                    coding=[
                        Coding( # type: ignore
                            system=Uri(
                                "http://terminology.hl7.org/CodeSystem/observation-category"
                            ),
                            code="laboratory",
                        ),
                        Coding( # type: ignore
                            system=Uri(
                                "http://terminology.hl7.org/CodeSystem/v2-0074"
                            ),
                            code="GE",
                        ),
                    ]
                )
            ],
            subject=self.dict_patient[sample],
            code=CodeableConcept(
                coding=[
                    Coding( # type: ignore
                        system=Uri(
                            "http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs"
                        ),
                        code="molecular-consequence",
                        display="Molecular Consequence",
                    )
                ]
            ),
        )
        mc_row.derivedFrom = [
            Reference( # type: ignore
                type=Uri("variant"), 
                identifier=Identifier(value=str(variant_id)) # type: ignore
            )
        ]
        return mc_row

    def add_hugo_to_mc_row(self, row, mc_row):
        gene_id = row["base__hugo"]
        if gene_id is not None:
            component = get_base_component_gene()
            component.valueCodeableConcept = CodeableConcept(
                coding=[
                    Coding( # type: ignore
                        system=Uri("https://www.genenames.org/geneId"), code=gene_id
                    )
                ]
            )
            mc_row.component = [component] # type: ignore
        else:
            mc_row.component = []

    def add_ensembl_to_mc_row(self, row, mc_row):
        component = get_base_component_transcript()
        component.valueCodeableConcept =\
                get_codeable_concept_ensembl_transcript(row["base__transcript"]) # type: ignore
        mc_row.component.append(component) # type: ignore

    def add_refseq_to_mc_row(self, row, mc_row):
        component = get_base_component_transcript()
        component.valueCodeableConcept =\
                get_codeable_concept_refseq_transcript(row["base__refseq"]) # type: ignore
        mc_row.component.append(component) # type: ignore

    def add_so_to_mc_row(self, row, mc_row):
        SO = row["base__so"]
        if SO == " " or SO is None or SO == "":
            return
        component = get_base_component_feature_consequence()
        component.valueCodeableConcept = CodeableConcept(
            coding=[
                Coding( # type: ignore
                    system="http://sequenceontology.org",
                    code=SO,
                    display=so_term_to_id_dict[SO],
                )
            ]
        )
        mc_row.component.append(component) # type: ignore

    def add_achange_to_mc_row(self, row, mc_row):
        achange = row["base__achange"]
        if achange == "" or achange == " ":
            return
        ensembl = row["base__transcript"]
        component = get_base_component_achange_hgvs()
        component.valueCodeableConcept =\
            get_codeable_concept_hgvs(ensembl, achange) # type: ignore
        mc_row.component.append(component) # type: ignore
        refseq = row["base__refseq"]
        component = get_base_component_achange_hgvs()
        component.valueCodeableConcept =\
                get_codeable_concept_hgvs(refseq, achange) # type: ignore
        mc_row.component.append(component) # type: ignore

    def add_cchange_to_mc_row(self, row, mc_row):
        cchange = row["base__cchange"]
        if cchange == "" or cchange == " ":
            return
        ensembl = row["base__transcript"]
        component = get_base_component_cchange_hgvs()
        component.valueCodeableConcept =\
                get_codeable_concept_hgvs(ensembl, cchange) # type: ignore
        mc_row.component.append(component)
        refseq = row["base__refseq"]
        component = get_base_component_cchange_hgvs()
        component.valueCodeableConcept =\
            get_codeable_concept_hgvs(refseq, cchange) # type: ignore
        mc_row.component.append(component) # type: ignore

    def get_molecular_consequence(self, variant_id, sample, row):
        mc_row = self.get_mc_row_init(variant_id, sample)
        self.add_hugo_to_mc_row(row, mc_row)
        self.add_ensembl_to_mc_row(row, mc_row)
        self.add_refseq_to_mc_row(row, mc_row)
        self.add_so_to_mc_row(row, mc_row)
        self.add_achange_to_mc_row(row, mc_row)
        self.add_cchange_to_mc_row(row, mc_row)
        return mc_row

    def parse_mapping(self, mapping_str: str) -> Tuple[str, str, List[str], str, str]:
        transcript: str = ""
        uniprot_id: str = ""
        list_so: List[str] = []
        achange: str = ""
        cchange: str = ""
        mapping_list = mapping_str.split(":")
        if len(mapping_list) < 6:
            return transcript, uniprot_id, list_so, achange, cchange
        transcript = mapping_list[0].strip()
        uniprot_id = mapping_list[1].strip()
        so = mapping_list[3].strip()
        list_so = so.split(",")
        achange = mapping_list[4].strip()
        cchange = mapping_list[5].strip()
        return transcript, uniprot_id, list_so, achange, cchange

    def add_fhir_mc_refseq_transcript(self, ensembl: str, mapping_comps: list):
        refseqs = self.ensembl_to_refseq_d.get(ensembl)
        if not refseqs:
            return
        for refseq in refseqs:
            component = get_base_component_transcript()
            component.valueCodeableConcept =\
                get_codeable_concept_refseq_transcript(refseq) # type: ignore
            mapping_comps.append(component)

    def add_fhir_mc_uniprot(self, uniprot_id: str, components: list):
        if not uniprot_id:
            return
        component = get_base_component_gene()
        component.valueCodeableConcept = CodeableConcept( # type: ignore
            text=f"{uniprot_id}" # type: ignore
        )
        components.append(component)

    def add_fhir_mc_sos(self, sos: List[str], components: list):
        component = get_base_component_feature_consequence()
        component.valueCodeableConcept = get_codeable_concept_so(sos) # type: ignore
        components.append(component)

    def add_fhir_mc_cchange_ensembl(self, cchange: str, ensembl: str, mapping_comps: list):
        if not cchange or not ensembl:
            return
        component = get_base_component_cchange_hgvs()
        component.valueCodeableConcept =\
            get_codeable_concept_hgvs(ensembl, cchange) # type: ignore
        mapping_comps.append(component)

    def add_fhir_mc_cchange_refseq(self, cchange: str, ensembl: str, mapping_comps: list):
        if not cchange or not ensembl:
            return
        refseqs = self.ensembl_to_refseq_d.get(ensembl)
        if not refseqs:
            return
        for refseq in refseqs:
            component = get_base_component_cchange_hgvs()
            component.valueCodeableConcept =\
                get_codeable_concept_hgvs(refseq, cchange) # type: ignore
            mapping_comps.append(component)

    def add_fhir_mc_achange_ensembl(self, achange: str, ensembl: str, mapping_comps: list):
        if not achange or not ensembl:
            return
        component = get_base_component_achange_hgvs()
        component.valueCodeableConcept =\
            get_codeable_concept_hgvs(ensembl, achange) # type: ignore
        mapping_comps.append(component)

    def add_fhir_mc_achange_refseqs(self, achange: str, ensembl: str, mapping_comps: list):
        if not achange or not ensembl:
            return
        refseqs = self.ensembl_to_refseq_d.get(ensembl)
        if not refseqs:
            return
        for refseq in refseqs:
            component = get_base_component_achange_hgvs()
            component.valueCodeableConcept =\
                get_codeable_concept_hgvs(refseq, achange) # type: ignore
            mapping_comps.append(component)

    def add_all_mappings(self, row, variant_id, sample, entry: list):
        mc_num = 1
        all_mappings = row["base__all_mappings"].split(";")
        for mapping in all_mappings:
            mc_num += 1
            transcript, uniprot_id, sos, achange, cchange = self.parse_mapping(mapping)
            if not transcript:
                continue
            if not sos or sos[0].lower() == "unknown":
                continue
            molecular_consequence = self.get_mc_row_init(variant_id, sample)
            molecular_consequence_id = self.uuid_maker(f"{self.base_id}{variant_id}{transcript}{mc_num}")
            molecular_consequence_uri = Uri(f"urn:uuid:{molecular_consequence_id}")
            components = []
            self.add_fhir_mc_ensembl_transcript(transcript, components)
            self.add_fhir_mc_refseq_transcript(transcript, components)
            self.add_fhir_mc_uniprot(uniprot_id, components)
            self.add_fhir_mc_sos(sos, components)
            self.add_fhir_mc_achange_ensembl(achange, transcript, components)
            self.add_fhir_mc_achange_refseqs(achange, transcript, components)
            self.add_fhir_mc_cchange_ensembl(cchange, transcript, components)
            self.add_fhir_mc_cchange_refseq(cchange, transcript, components)
            molecular_consequence.component = components
            molecular_consequence_bundle = BundleEntry( # type: ignore
                resource=molecular_consequence, fullUrl=molecular_consequence_uri
            )
            entry.append(molecular_consequence_bundle)
