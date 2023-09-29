import uuid
import json
import requests
import hashlib
import sqlite3
from pathlib import Path
from oakvar import BaseReporter
from fhir.resources.patient import Patient
from fhir.resources.observation import Observation
from fhir.resources.observation import ObservationComponent
from fhir.resources.humanname import HumanName
from fhir.resources.codeableconcept import CodeableConcept
from fhir.resources.coding import Coding
from fhir.resources.reference import Reference
from fhir.resources.bundle import Bundle, BundleEntry
from fhir.resources.fhirtypes import (
    Uri,
    MetaType,
    IdentifierType,
    String,
    RangeType,
    Code,
    BackboneElementType,
)
from fhir.resources.quantity import Quantity
from fhir.resources.identifier import Identifier
from fhir.resources.list import List, ListEntry


class Reporter(BaseReporter):
    SO_dict = {
        "start_lost": "SO:0002012",
        "intron_variant": "SO:0001627",
        "frameshift_elongation": "SO:0001909",
        "missense_variant": "SO:0001583",
        "stop_lost": "SO:0001578",
        "frameshift_truncation": "SO:0001910",
        "lnc_RNA": "SO:0002127",
        "inframe_insertion": "SO:0001821",
        "5_prime_UTR_variant": "SO:0001623",
        "synonymous_variant": "SO:0001819",
        "splice_site_variant": "SO:0001629",
        "3_prime_UTR_variant": "SO:0001624",
        "2kb_upstream_variant": "SO:0001636",
        "2kb_downstream_variant": "SO:0002083",
        "stop_gained": "SO:0001587",
        "retained_intron": "SO:0002113",
        "inframe_deletion": "SO:0001822",
        "misc_RNA": "SO:0000673",
        "complex_substitution": "SO:1000005",
        "processed_transcript": "SO:0001503",
        "transcribed_unprocessed_pseudogene": "SO:0002107",
        "unprocessed_pseudogene": "SO:0001760",
        "miRNA": "SO:0000276",
        "processed_pseudogene": "SO:0000043",
        "snRNA": "SO:0000274",
        "transcribed_processed_pseudogene": "SO:0002109",
        "NMD_transcript_variant": "SO:0001621",
        "unconfirmed_transcript": "SO:0002139",
        "pseudogene": "SO:0000336",
        "transcribed_unitary_pseudogene": "SO:0002108",
        "NSD_transcript": "SO:0002130",
        "snoRNA": "SO:0000275",
        "scaRNA": "SO:0002095",
        "unitary_pseudogene": "O:0001759",
        "polymorphic_pseudogene": "SO:0001841",
        "rRNA": "SO:0000252",
        "IG_V_pseudogene": "SO:0002102",
        "ribozyme": "SO:0000374",
        "TR_V_gene": "SO:0002137",
        "TR_V_pseudogene": "SO:0002103",
        "TR_D_gene": "SO:0002135",
        "TR_J_gene": "SO:0002136",
        "TR_C_gene": "SO:0002134",
        "TR_J_pseudogene": "SO:0002104",
        "IG_C_gene": "SO:0002123",
        "IG_C_pseudogene": "SO:0002100",
        "IG_J_gene": "SO:0002125",
        "IG_J_pseudogene": "SO:0002101",
        "IG_D_gene": "SO:0002124",
        "IG_V_gene": "SO:0002126",
        "translated_processed_pseudogene": "SO:0002105",
        "scRNA": "SO:0000013",
        "vault_RNA": "SO:0000404",
        "translated_unprocessed_pseudogene": "SO:0002106",
        "Mt_tRNA": "SO:0002129",
        "Mt_rRNA": "SO:0002128",
        "start_retained_variant": "SO:0002019",
        "stop_retained_variant": "SO:0001567",
        "exon_loss_variant": "SO:0001572",
        "transcript_ablation": "SO:0001893",
        "pseudogene_rRNA": "SO:0002111",
        "sRNA": "SO:0002352",
    }

    # def __init__(self):
    #    #super().__init__()
    #    self.dict_entries = None
    #    self.dict_bundles = None
    #    self.dict_patient = None
    #    self.dict_nums = None
    #
    #    self.prefix = None
    #    self.wf = None

    def setup(self):
        self.levels_to_write = self.conf.get("pages", "variant")
        self.filenames = []
        self.counter = 0

        self.samples = []
        self.bundles = []

        # establish filename with fhir suffix
        self.prefix = self.savepath

        # get sample names
        conn = sqlite3.connect(self.dbpath)
        curs = conn.cursor()
        curs.execute("SELECT DISTINCT base__sample_id FROM 'sample' ")
        for sample in curs.fetchall():
            if sample[0].count(",") < 1:
                self.samples.append(sample[0])

        def create_dict(keys):
            unique_dict = {}
            for key in keys:
                if key not in unique_dict:
                    unique_dict[key] = []
            return unique_dict

        self.dict_entries = create_dict(self.samples)
        self.dict_bundles = create_dict(self.samples)
        self.dict_patient = create_dict(self.samples)
        self.dict_nums = create_dict(self.samples)

        # get number of rows
        curs = conn.cursor()
        curs.execute("SELECT COUNT(*) from variant")
        self.num_rows = curs.fetchone()[0]

        # get str for id generation
        curs = conn.cursor()
        curs.execute('select colval from info where colkey="input_paths"')
        # get input_path and split is so that only path is part of id.
        self.str_id = curs.fetchone()[0].split(" ", 1)[-1]

        curs = conn.cursor()
        curs.execute('select colval from info where colkey="annotators"')
        self.str_id += curs.fetchone()[0][1:-1]
        self.str_id = self.str_id[1:-1]
        self.str_id = self.str_id[-32:]

        curs = conn.cursor()
        curs.execute('select colval from info where colkey="mapper"')
        self.str_id += curs.fetchone()[0]
        # make this bundle UUID

        for sample in self.samples:
            # create sample id
            sample_name = f"{self.str_id} + {sample}"
            hex_sample = hashlib.md5(sample_name.encode("utf-8")).hexdigest()

            # create patient
            sample_2_patient = Patient()
            name = HumanName()
            name.use = "official"
            name.given = [sample]
            sample_2_patient.name = [name]
            patient_id = str(uuid.UUID(hex=hex_sample))
            # add to patient dict for reference in row observations
            # if for future cases a reference is needed it is here otherwise for the MVP it is not needed
            subject = Reference(type="Patient")
            subject.resource_type = "Reference"
            subject.reference = f"urn:uuid:{patient_id}"
            self.dict_patient[sample] = subject

            patient_entry = BundleEntry(
                resource=sample_2_patient, fullUrl=f"urn:uuid:{patient_id}"
            )
            self.dict_entries[sample].append(patient_entry)

            # create sample bundle and and it to dictionary
            name = f"bundle + {sample} + self.str_id"
            hex_sample = hashlib.md5(name.encode("utf-8")).hexdigest()
            bundle = Bundle(
                type="collection",
                identifier=Identifier(value=str(uuid.UUID(hex=hex_sample))),
            )
            self.dict_bundles[sample] = bundle
            # self.dict_entries[sample].append(BundleEntry(resource=subject))

        # create CodingResource for row ObservationResources to Use
        self.fhir_system = Uri("http://loinc.org")

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

    def end(self):
        if self.module_options.get("all_transcripts") == "true":
            for sample in self.samples:
                # fill in bundles with entries
                self.dict_bundles[sample].entry = self.dict_entries[sample]
                obs = self.dict_bundles[sample]
                filename = str(self.prefix) + f"__{sample}.all.fhir.json"
                self.wf = open(filename, "w", encoding="utf-8")
                json_str = obs.json(indent=2)
                self.wf.write(json_str)
                self.filenames.append(filename)
            self.wf.close()
        else:
            for sample in self.samples:
                filename = str(self.prefix) + f"__{sample}.primary.fhir.json"
                self.wf = open(filename, "w", encoding="utf-8")
                self.dict_bundles[sample].entry = self.dict_entries[sample]
                obs = self.dict_bundles[sample]
                json_str = obs.json(indent=2)
                self.wf.write(json_str)
                self.filenames.append(filename)
            self.wf.close()

        return self.filenames

        # assign entries to correct bundles
        # self.bundle.entry = self.entries

        # create a json_str from FHIR BundleResource
        # for x in self.bundles:
        #    json_str = x.json(indent=2)
        #    self.wf.write(json_str)
        # return[str(self.savepath)]
        # json_str = self.bundles.json(indent=2)

        # write json_file
        # self.wf.write(json_str)

        # self.wf.close()
        # return [str(self.savepath)]

    def write_preface(self, level: str):
        if level not in self.levels_to_write:
            return

    def write_table_row(self, row):
        # get samples that have variant(row)
        sample_with_variants = row["tagsampler__samples"].split(",")

        for sample in sample_with_variants:
            # create codingType for row  Variant Observation
            coding = Coding()
            coding.system = Uri("http://loinc.org")
            coding.code = "69548-6"
            code = CodeableConcept()
            code.coding = [coding]
            # create CC for Laboratory
            cat_coding = Coding()
            cat_coding.system = Uri(
                "http://terminology.hl7.org/CodeSystem/observation-category"
            )
            cat_coding.code = "laboratory"
            cat_cc = CodeableConcept()
            cat_cc.coding = [cat_coding]

            #get primary transcript ENSEMBL
            primary_transcript = row['base__transcript']
            coding_primary_transcript = Coding()
            coding_primary_transcript.system = "http://loinc.org"
            coding_primary_transcript.code = "51958-7"
            coding_primary_transcript.display = "Transcript reference sequence [ID]"
            code_transcript = CodeableConcept(coding=[coding_primary_transcript])
            comp_primary_transcript = ObservationComponent(code=code_transcript)
            coding_pt_comp = Coding()
            coding_pt_comp.system = "http://www.ensembl.org"
            coding_pt_comp.code = primary_transcript
            comp_primary_transcript.valueCodeableConcept = CodeableConcept(
                coding=[coding_pt_comp]
            )
            # Get Alleles from sqlite file
            ref = row["base__ref_base"]
            alt = row["base__alt_base"]

            # Make Component for reference allele (UNIVERSAL COMPONENT)
            coding_ref = Coding()
            coding_ref.system = Uri("http://loinc.org")
            coding_ref.code = "69547-8"  # always code for reference allele
            coding_ref.display = "Ref nucleotide"
            code_ref = CodeableConcept()
            code_ref.coding = [coding_ref]
            comp_ref = ObservationComponent(code=code_ref)
            comp_ref.valueString = ref


            # Make Component for (alt)ernate allele (UNIVERSAL COMPONENT)
            coding_alt = Coding()
            coding_alt.system = Uri("http://loinc.org")
            coding_alt.code = "69551-0"
            coding_alt.display = "Alt allele"
            code_alt = CodeableConcept()
            code_alt.coding = [coding_alt]
            comp_alt = ObservationComponent(code=code_alt)
            comp_alt.valueString = alt


            # get chrom and pos information
            chrom_location = row["base__chrom"]
            pos = row["base__pos"]


            #Make Component for Chrom Location
            coding_chrom = Coding()
            coding_chrom.system = Uri("http://loinc.org")
            coding_chrom.code = "48001-2"
            coding_chrom.display = "Cytogenetic (chromosome) location"
            code_chrom = CodeableConcept()
            code_chrom.coding = [coding_chrom]
            comp_chrom = ObservationComponent(code=code_chrom)
            comp_chrom.valueCodeableConcept = CodeableConcept(
                text=String(f"{chrom_location}")
            )
            #
            ##Make Component for character based counting 
            coding_counting = Coding()
            coding_counting.system = Uri("http://loinc.org")
            coding_counting.code = "92822-6"
            code_counting = CodeableConcept()
            code_counting.coding = [coding_counting]
            comp_counting = ObservationComponent(code=code_counting)
            cc_counting = CodeableConcept(
                coding=[
                    (
                        Coding(
                            system=Uri("http://loinc.org"),
                            code="LA30102-0",
                            display=String("1-based character counting"),
                        )
                    )
                ]
            )
            comp_counting.valueCodeableConcept = cc_counting
    
            # Make Component for Start and End (sne) ##contains position 
            coding_st_sne = Coding()
            coding_st_sne.system = Uri("http://loinc.org")
            coding_st_sne.code = "81254-5"
            coding_st_sne.display = "Genomic allele start-end"
            st_sne_value_low = Quantity(value=pos)
            st_sne_value_high = Quantity(value=row["base__pos_end"])
            code_sne = CodeableConcept(coding=[coding_st_sne])
            comp_sne = ObservationComponent(code=code_sne)
            comp_sne.valueRange = RangeType(
                low=st_sne_value_low, high=st_sne_value_high
            )
            # create Observation Resource for row


            #create single observation for row
            obs_row = Observation(
                status="final",
                code=code,
                subject=self.dict_patient[sample],
                category=[cat_cc],
            )
            
            #make ID for single observation
            conn = sqlite3.connect(self.dbpath)
            curs = conn.cursor()
            curs.execute("SELECT COUNT(*) from variant")
            self.num_rows = curs.fetchone()[0]
            id_maker = (
                (self.str_id)
                + f"{str(row['base__chrom']) + str(row['base__pos']) + str(row['base__ref_base']) + str(row['base__alt_base'])}"
            )
            id = self.uuid_maker(id_maker + self.str_id)
            uri_maker = Uri(f"urn:uuid:{id}")


            # obs_row.meta = MetaType()
            # obs_row.meta.profile = ["http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/variant"]

            # Make Component for Gene ID
            gene_id = row["base__hugo"]
            if gene_id is not None:
                coding_id = Coding()
                coding_id.system = Uri("http://loinc.org")
                coding_id.code = "48018-6"
                coding_id.display = "Gene studied [ID]"
                code_gene_id = CodeableConcept(coding=[coding_id])
                comp_gene_id = ObservationComponent(code=code_gene_id)
                comp_gene_id.valueCodeableConcept = CodeableConcept(
                    text=f"{gene_id}"
                )
                obs_row.component = [comp_gene_id,comp_primary_transcript, comp_ref, comp_alt, comp_chrom, comp_counting,comp_sne]
            else:
                obs_row.component = [comp_primary_transcript,comp_ref, comp_alt, comp_chrom, comp_counting, comp_sne]


    
            SO = row["base__so"]
            if SO != ' ' and SO!= None or SO != '' and SO != None:
                SO_coding = Coding()
                SO_coding.system = "http://hl7.org/fhir/uv/genomics-reporting/STU2/CodeSystem-tbd-codes-cs"
                SO_coding.code = "feature-consequence"
                code_SO = CodeableConcept(coding = [SO_coding])
                comp_SO = ObservationComponent(code=code_SO)
                
                comp_SO.valueCodeableConcept = CodeableConcept(
                    coding=[Coding(system="http://sequenceontology.org",
                                  code=SO,
                                  display=self.SO_dict[SO])
                    ])
                obs_row.component.append(comp_SO)




            aa_change = row["base__achange"]
            c_change = row["base__cchange"]


            # Make Component for achange (change)
            if aa_change != "" and aa_change != " ":
                coding_change = Coding()
                coding_change.system = Uri("http://loinc.org")
                coding_change.code = "48006-1"
                coding_change.display = "Amino Acid Change [type]"
                code_achange = CodeableConcept(coding=[coding_change])
                comp_achange = ObservationComponent(code=code_achange)
                comp_achange.valueCodeableConcept = CodeableConcept(
                    text=f"{primary_transcript}:{aa_change}"
                )
                obs_row.component.append(comp_achange)
            #Make Component for cchange (change)
            if c_change != '' and c_change != ' ':
                coding_c_change = Coding()
                coding_c_change.system = Uri("http://loinc.org")
                coding_c_change.code = "48004-6"
                coding_c_change.display = "DNA change (c.HGVS)"
                code_c_change = CodeableConcept(coding=[coding_c_change])
                comp_c_change = ObservationComponent(code=code_c_change)
                comp_c_change.valueCodeableConcept = CodeableConcept(
                    text=f"{primary_transcript}:{c_change}"
                )
                obs_row.component.append(comp_c_change)



            # add primary observation to bundle
            converted_ent = BundleEntry(resource=obs_row, fullUrl=uri_maker)
            self.dict_entries[sample].append(converted_ent)

            #begin all_transcript module optionloop
            if self.module_options.get("all_transcripts") == "true":
                all_mappings = row["base__all_mappings"].split(";")


                
                for mapping in all_mappings:
                    mapping_comps = []
                    obs_mapping = Observation(
                        status="final",
                        code=code,
                        subject=self.dict_patient[sample],
                        category=[cat_cc],
                    )
                    mapping_list = mapping.split(":")
                    transcript = mapping_list[0].strip()
                    id_maker = (
                        (self.str_id)
                        + f"{str(row['base__chrom']) + str(row['base__pos']) + str(row['base__ref_base']) + str(row['base__alt_base'] + transcript)}"
                    )
                    mapping_id = self.uuid_maker(id_maker + self.str_id)
                    uri_maker = Uri(f"urn:uuid:{mapping_id}")

                    # Ensembl transcript component
                    if len(transcript) > 0:
                        coding_transcript = Coding()
                        coding_transcript.system = "http://loinc.org"
                        coding_transcript.code = "51958-7"
                        coding_transcript.display = "Transcript reference sequence [ID]"
                        code_transcript = CodeableConcept(coding=[coding_transcript])
                        comp_transcript = ObservationComponent(code=code_transcript)
                        coding_comp = Coding()
                        coding_comp.system = "http://www.ensembl.org"
                        coding_comp.code = transcript
                        comp_transcript.valueCodeableConcept = CodeableConcept(
                            coding=[coding_comp]
                        )
                        mapping_comps.append(comp_transcript)
                        if len(mapping_list) > 1:
                            uniprot_id = mapping_list[1].strip()
                            sequence_ontology = mapping_list[3].strip()
                            list_so = sequence_ontology.split(",")
                            amino_acid_change = mapping_list[4].strip()
                            chromosome_change = mapping_list[5].strip()
    
                        #create SO system for SO list
                        so_system_coding = Coding()
                        so_system_coding.system = "http://hl7.org/fhir/uv/genomics-reporting/STU2/CodeSystem-tbd-codes-cs"
                        so_system_coding.code = "feature-consequence"
                        code_system = CodeableConcept(coding = [so_system_coding])
                        comp_SO_all = ObservationComponent(code=code_system)
    
                        coding_list = []
                            
                        for so in list_so:
                            if so != "unknown" and so != "" and so != " ":
                                so_coding = Coding()
                                so_coding.system = "http://sequenceontology.org"
                                so_coding.code = self.SO_dict[so]
                                so_coding.display = so
                                coding_list.append(so_coding)
    
    
    
                        comp_SO_all.valueCodeableConcept = CodeableConcept(
                                coding=coding_list
                        )
    
                        mapping_comps.append(comp_SO_all)
                        
                        # make a_change component
                        if  amino_acid_change != "" and amino_acid_change != " ":
                        
                            coding_change = Coding()
                            coding_change.system = Uri("http://loinc.org")
                            coding_change.code = "48006-1"
                            coding_change.display = "Amino Acid Change [type]"
                            code_achange = CodeableConcept(coding=[coding_change])
                            comp_achange = ObservationComponent(code=code_achange)
                            comp_achange.valueCodeableConcept = CodeableConcept(
                                text=f"{transcript}:{amino_acid_change}"
                            )
                            mapping_comps.append(comp_achange)
                        #
                        # make c_change component (ENSEMBL)
                        if chromosome_change != '' and chromosome_change != ' ':
                            coding_c_change = Coding()
                            coding_c_change.system = Uri("http://loinc.org")
                            coding_c_change.code = "48004-6"
                            coding_c_change.display = "DNA change (c.HGVS)"
                            code_c_change = CodeableConcept(coding=[coding_c_change])
                            comp_c_change = ObservationComponent(code=code_c_change)
                            comp_c_change.valueCodeableConcept = CodeableConcept(
                                text=f"{transcript}:{chromosome_change}"
                            )
                            mapping_comps.append(comp_c_change)
                            #
                            ##make rc_change component (RefSeq)
                            coding_rc_change = Coding()
                            coding_rc_change.system = Uri("http://loinc.org")
                            coding_rc_change.code = "48004-6"
                            coding_rc_change.display = "DNA change (c.HGVS)"
                            code_rc_change = CodeableConcept(coding=[coding_rc_change])
                            comp_rc_change = ObservationComponent(code=code_rc_change)
                            comp_rc_change.valueCodeableConcept = CodeableConcept(
                                text=f"{row['base__refseq']}:{chromosome_change}"
                            )
                            mapping_comps.append(comp_rc_change)
                        obs_mapping.component = mapping_comps
                        mapping_ent = BundleEntry(resource=obs_mapping, fullUrl=uri_maker)
                        self.dict_entries[sample].append(mapping_ent)
