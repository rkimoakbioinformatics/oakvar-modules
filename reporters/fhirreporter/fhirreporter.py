import uuid
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
from fhir.resources.fhirtypes import Uri, MetaType, IdentifierType, String, RangeType
from fhir.resources.quantity import Quantity
from fhir.resources.identifier import Identifier


class Reporter(BaseReporter):
    
    
    def setup(self):
        # establish filename with fhir suffix
        self.prefix = self.savepath
        self.wf = None
        self.filenames = []
        self.counter = 0

        self.levels_to_write = self.confs.get("pages", "variant")
        self.samples = []
        self.bundles = []
        self.SO_dict = {
            "start_lost" :"SO:0002012",
            "intron_variant" : "SO:0001627",
            "frameshift_elongation" : "SO:0001909",
            "missense_variant" : "SO:0001583",
            "stop_lost" : "SO:0001578",
            "frameshift_truncation" : "SO:0001910",
            'lnc_RNA': "SO:0002127",
            "inframe_insertion" : "SO:0001821",
            "5_prime_UTR_variant" : "SO:0001623",
            "synonymous_variant" : "SO:0001819",
            "splice_site_variant" : "SO:0001629",
            "3_prime_UTR_variant" : "SO:0001624",
            "2kb_upstream_variant" : "SO:0001636",
            "2kb_downstream_variant":"SO:0002083",
            "stop_gained" : "SO:0001587",
            "retained_intron":"SO:0002113",
            "inframe_deletion" : "SO:0001822",
            "misc_RNA" : "SO:0000673",
            "complex_substitution" : "SO:1000005",
            "processed_transcript" : "SO:0001503",
            "transcribed_unprocessed_pseudogene" : "SO:0002107",
            "unprocessed_pseudogene": "SO:0001760",
            "miRNA" : "SO:0000276",
            "processed_pseudogene" : "SO:0000043",
            "snRNA":"SO:0000274",
            "transcribed_processed_pseudogene":"SO:0002109",
            "NMD_transcript_variant" : "SO:0001621",
            "unconfirmed_transcript": "SO:0002139",
            "pseudogene" :"SO:0000336",
            "transcribed_unitary_pseudogene" : "SO:0002108",
            "NSD_transcript" : "SO:0002130",
            "snoRNA":"SO:0000275",
            "scaRNA" : "SO:0002095",
            "unitary_pseudogene" : "O:0001759",
            "polymorphic_pseudogene" : "SO:0001841",
            "rRNA" : "SO:0000252",
            "IG_V_pseudogene" : "SO:0002102",
            "ribozyme" : "SO:0000374",
            "TR_V_gene":"SO:0002137",
            "TR_V_pseudogene" : "SO:0002103",
            "TR_D_gene" : "SO:0002135",
            "TR_J_gene" : "SO:0002136",
            "TR_C_gene" : "SO:0002134",
            "TR_J_pseudogene" : "SO:0002104",
            "IG_C_gene" : "SO:0002123" , 
            "IG_C_pseudogene" : "SO:0002100",
            "IG_J_gene" : "SO:0002125",
            "IG_J_pseudogene" : "SO:0002101",
            "IG_D_gene" : "SO:0002124",
            "IG_V_gene" : "SO:0002126",
            "translated_processed_pseudogene" : "SO:0002105",
            "scRNA" : "SO:0000013",
            "vault_RNA":"SO:0000404",
            "translated_unprocessed_pseudogene" : "SO:0002106",
            "Mt_tRNA" : "SO:0002129",
            "Mt_rRNA" : "SO:0002128",
            "start_retained_variant" : "SO:0002019",
            "stop_retained_variant" : "SO:0001567",
            "exon_loss_variant" : "SO:0001572",
            "transcript_ablation" : "SO:0001893",
            "pseudogene_rRNA" : "SO:0002111",
            "sRNA" : "SO:0002352"
        }
        # get sample names
        conn = sqlite3.connect(self.dbpath)
        curs = conn.cursor()
        curs.execute("SELECT DISTINCT base__sample_id FROM 'sample' ")
        for sample in curs.fetchall():
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
            # print(type(sample))
            # id_er = uuid.UUID(sample)
            # print(id_er)

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
        for sample in self.samples:
            filename = str(self.prefix) + f"__{sample}.fhir.json"
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
        sample_with_variants = row["tagsampler__samples"].split(":")

        for sample in sample_with_variants:
            self.dict_nums[sample].append(self.counter)

            #get RefSeq from sqlite file 
            refseq = row["base__refseq"]


            # Get Alleles from sqlite file
            ref = row["base__ref_base"]
            alt = row["base__alt_base"]

            #get chrom and pos information
            chrom_location = row["base__chrom"]
            pos = row["base__pos"]




            # create codingType for row  Variant Observation
            coding = Coding()
            coding.system = Uri("http://loinc.org")
            coding.code = "69548-6"
            code = CodeableConcept()
            code.coding = [coding]

            #create CC for Laboratory 
            cat_coding = Coding()
            cat_coding.system = Uri("http://terminology.hl7.org/CodeSystem/observation-category")
            cat_coding.code = "laboratory"
            cat_cc = CodeableConcept()
            cat_cc.coding=[cat_coding]


            # create Observation Resource for row
            obs_row = Observation(
                status="final", code=code, subject=self.dict_patient[sample], category=[cat_cc]
            )
            #obs_row.meta = MetaType()
            #obs_row.meta.profile = ["http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/variant"]

            # Make Component for reference allele
            coding_ref = Coding()
            coding_ref.system = Uri("http://loinc.org")
            coding_ref.code = "69547-8"  # always code for reference allele
            coding_ref.display = "Ref nucleotide"
            code_ref = CodeableConcept()
            code_ref.coding = [coding_ref]
            comp_ref = ObservationComponent(code=code_ref)
            comp_ref.valueString = ref

            # Make Component for (alt)ernate allele
            coding_alt = Coding()
            coding_alt.system = Uri("http://loinc.org")
            coding_alt.code = "69551-0"
            code_alt = CodeableConcept()
            code_alt.coding = [coding_alt]
            comp_alt = ObservationComponent(code=code_alt)
            comp_alt.valueString = alt

            ##Make Component for Chrom 
            coding_chrom = Coding()
            coding_chrom.system = Uri("http://loinc.org")
            coding_chrom.code = "48001-2"
            coding_chrom.display = "Cytogenetic (chromosome) location"
            code_chrom = CodeableConcept()
            code_chrom.coding = [coding_chrom]
            comp_chrom = ObservationComponent(code=code_chrom)
            comp_chrom.valueCodeableConcept = CodeableConcept(text=String(f"{chrom_location}"))
#
            ##Make Component for Pos 
            coding_pos = Coding()
            coding_pos.system = Uri("http://loinc.org")
            coding_pos.code = "92822-6"
            code_pos = CodeableConcept()
            code_pos.coding = [coding_pos]
            comp_pos = ObservationComponent(code=code_pos)
            cc_pos = CodeableConcept(coding=
                                     [(Coding(system=Uri("http://loinc.org"),
                                             code="LA30102-0",
                                             display=String("1-based character counting")
            ))])
            comp_pos.valueCodeableConcept = cc_pos

            #Make Component for Sequence Ontology 

            SO = row["base__so"]
            comp_so = None
            if SO is not None:
                so_val_coding = Coding()
                so_val_coding.system = ("http://sequenceontology.org")
                so_val_coding.code = self.SO_dict[SO]
                so_val_coding.display = SO
                so_value = CodeableConcept(coding=[so_val_coding])
                comp_so = ObservationComponent(code=so_value)


            #Make Component for Start and End (sne)

            coding_st_sne = Coding()
            coding_st_sne.system = Uri("http://loinc.org")
            coding_st_sne.code = "81254-5"
            coding_st_sne.display = "Genomic allele start-end"
            st_sne_value_low = Quantity(value=row["base__pos"])
            st_sne_value_high = Quantity(value=row["base__pos_end"])
            code_sne = CodeableConcept(coding=[coding_st_sne])
            comp_sne = ObservationComponent(code=code_sne)
            comp_sne.valueRange= RangeType(low=st_sne_value_low, high = st_sne_value_high)

            #Make Component for cchange (change)
            aa_change = row["base__achange"]
            c_change = row["base__cchange"]

            coding_change = Coding()
            coding_change.system = Uri("http://loinc.org")
            coding_change.code = '48006-1'
            coding_change.display  = "Amino Acid Change [type]"
            code_change = CodeableConcept(coding=[coding_change])
            comp_change = ObservationComponent(code=code_change)
            comp_change.valueCodeableConcept = CodeableConcept(text=f"{aa_change}")

            coding_c_change = Coding()
            coding_c_change.system = Uri("http://loinc.org")
            coding_c_change.code = "48004-6"
            coding_c_change.display = "DNA change (c.HGVS)"
            code_c_change = CodeableConcept(coding=[coding_c_change])
            comp_c_change = ObservationComponent(code=code_c_change)
            comp_c_change.valueCodeableConcept = CodeableConcept(text=f"{refseq}:{c_change}")

            
        




            

            # add componenets to row observation
            if comp_so is not None:
                obs_row.component = [comp_ref, comp_alt,comp_chrom,comp_pos, comp_so, comp_sne,comp_change,comp_c_change]
            else:obs_row.component = [comp_ref, comp_alt,comp_chrom,comp_pos, comp_sne,comp_change,comp_c_change]

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
            converted_ent = BundleEntry(resource=obs_row, fullUrl=uri_maker)

            self.dict_entries[sample].append(converted_ent)
