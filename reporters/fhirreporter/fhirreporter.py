import uuid
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
from fhir.resources.fhirtypes import Uri, MetaType, IdentifierType, String
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


            # create codingType for row Observation
            coding = Coding()
            coding.system = Uri("http://loinc.org")
            coding.code = "8480-6"
            code = CodeableConcept()
            code.coding = [coding]

            # Get Alleles from sqlite file
            ref = row["base__ref_base"]
            alt = row["base__alt_base"]

            #get chrom and pos information
            chrom_location = row["base__chrom"]
            print(chrom_location)
            pos = row["base__pos"]
            print(pos)

            # create Observation Resource for row
            obs_row = Observation(
                status="final", code=code, subject=self.dict_patient[sample]
            )
            # obs_row.meta = MetaType()
            # obs_row.meta.profile = ["http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/variant"]

            # Make Component for reference allele
            coding_ref = Coding()
            coding_ref.system = Uri("http://loinc.org")
            coding_ref.code = "69547-8"  # always code for reference allele
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
            

            # add componenets to row observation
            obs_row.component = [comp_ref, comp_alt,comp_chrom,comp_pos]

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
