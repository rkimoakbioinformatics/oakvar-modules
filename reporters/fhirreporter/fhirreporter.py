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
from fhir.resources.fhirtypes import Uri


class Reporter(BaseReporter):
    def setup(self):
        # establish filename with fhir suffix
        self.prefix = self.savepath
        self.wf = None
        self.filenames = []
        print(self.savepath)
        if self.savepath == None:
            self.savepath = Path("oakvar_result.json")
        else:
            if self.savepath.suffix != ".json":
                self.savepath = Path(str(self.prefix) + ".fhir.json")

        self.wf = open(self.savepath, "w", encoding="utf-8")
        self.levels_to_write = self.confs.get("pages", "variant")

        # create FHIR bundle resource
        self.bundle = Bundle(type="collection")

        # get patient name
        conn = sqlite3.connect(self.dbpath)
        curs = conn.cursor()
        curs.execute("SELECT 'primary_transcript' FROM 'info' ")
        patient_name = curs.fetchone()[0]

        # get number of rows
        curs = conn.cursor()
        curs.execute("SELECT COUNT(*) from variant")
        self.num_rows = curs.fetchone()[0]

        # get str for id generation
        curs = conn.cursor()
        curs.execute('select colval from info where colkey="input_paths"')
        # get input_path and split is so that only path is part of id.
        self.str_id = curs.fetchone()[0].split(" ", 1)[-1]

        curs.execute('select colval from info where colkey="annotators"')
        self.str_id += curs.fetchone()[0][1:-1]
        self.str_id = self.str_id[1:-1]
        self.str_id = self.str_id[-32:]
        self.obs_list = []
        self.entries = []
        self.counter = 0

        # create and fill in PatientResource
        self.patient = Patient()
        name = HumanName()
        name.use = "official"
        name.given = [patient_name]
        self.patient.name = [name]
        id = hashlib.md5(self.str_id.encode("utf-8")).hexdigest()
        self.patient.id = f"{str(uuid.UUID(hex=id))}"
        patient_entry = BundleEntry(
            resource=self.patient, fullUrl=f"urn:uuid:{str(uuid.UUID(hex=id))}"
        )
        self.entries.append(patient_entry)

        # create Reference resource for PatientResource
        self.subject = Reference(type="Patient")
        self.subject.reference = f"urn:uuid:{self.patient.id}"

        # create CodingResource for row ObservationResources to Use
        coding = Coding()
        coding.system = Uri("http://loinc.org")
        self.code = CodeableConcept()
        self.code.coding = [coding]

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
        print(self.savepath)

        self.wf.close()
        return [str(self.savepath)]

    def write_preface(self, level: str):
        if level not in self.levels_to_write:
            return

    def write_table_row(self, row):
        # create codingType for row Observation
        coding = Coding()
        coding.system = Uri("http://loinc.org")
        coding.code = "8480-6"
        code = CodeableConcept()
        code.coding = [coding]

        # Get Alleles from sqlite file
        ref = row["base__ref_base"]
        alt = row["base__alt_base"]

        # create Observation Resource for row
        obs_row = Observation(status="final", code=code, subject=self.subject)

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

        # add componenets to row observation
        obs_row.component = [comp_ref, comp_alt]

        # add row observation to list of observations
        self.obs_list.append(obs_row)

        self.counter += 1
        id_maker = (self.str_id) + str(100 * self.counter)
        id = self.uuid_maker(id_maker + self.str_id)
        uri_maker = Uri(f"urn:uuid:{id}")
        converted_ent = BundleEntry(resource=obs_row, fullUrl=uri_maker)
        self.entries.append(converted_ent)

        # check if all rows are done
        if len(self.obs_list) == self.num_rows:
            self.bundle.entry = self.entries

            # create a json_str from FHIR BundleResource
            json_str = self.bundle.json(indent=2)

            # write json_file
            self.wf.write(json_str)
