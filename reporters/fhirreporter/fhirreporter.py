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
from fhir.resources.fhirtypes import Uri, MetaType, IdentifierType
from fhir.resources.identifier import Identifier


class Reporter(BaseReporter):
    def setup(self):
        # establish filename with fhir suffix
        self.prefix = self.savepath
        self.wf = None
        self.filenames = []



        self.levels_to_write = self.confs.get("pages", "variant")
        self.samples = []
        self.bundles = []

        # get sample names
        conn = sqlite3.connect(self.dbpath)
        curs = conn.cursor()
        curs.execute("SELECT DISTINCT base__sample_id FROM 'sample' ")
        for sample in curs.fetchall():
            self.samples.append(sample[0])
        print(self.samples)
        self.dict_bundles = dict.fromkeys(self.samples)
        self.dict_entries = dict.fromkeys(self.samples,[])
        self.dict_patient = dict.fromkeys(self.samples)


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
        #make this bundle UUID


        for sample in self.samples:
            #print(type(sample))
            #id_er = uuid.UUID(sample)
            #print(id_er)
            print(sample)
            
            #create sample id 
            sample_name = f"{self.str_id} + {sample}"
            hex_sample = hashlib.md5(sample_name.encode("utf-8")).hexdigest()


            #create patient 
            sample_2_patient = Patient()
            name = HumanName()
            name.use = "offical"
            name.given = [sample]
            sample_2_patient.name  = [name]
            patient_id = hex_sample
            #add to patient dict for reference in row observations 
            subject = Reference(type="Patient")
            subject.reference = f"urn:uuid:{patient_id}"
            self.dict_patient[sample] = subject

            patient_entry = BundleEntry(resource = sample_2_patient)
            self.dict_entries[sample].append(patient_entry)



            #create sample bundle and and it to dictionar
            name = f"bundle + {sample}" 
            bundle= Bundle(type="collection",identifier=Identifier(value=str(uuid.UUID(hex=hex_sample))))
            self.dict_bundles[sample] = bundle
            self.dict_entries[sample].append(BundleEntry(resource=subject))
        print(self.dict_patient["s0"])


            
        self.counter = 0


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
        print(len(self.dict_entries["s0"]))
        print(len(self.dict_entries["s1"]))
        print(self.dict_entries.keys())

        for sample in self.samples:
            filename = str(self.prefix) + f"__{sample}.fhir.json"
            self.wf = open(filename,"w",encoding="utf-8")
            self.dict_bundles[sample].entry = self.dict_entries[sample]
            obs = self.dict_bundles[sample]
            json_str = obs.json(indent=2)
            self.wf.write(json_str)
            print("works")
            self.filenames.append(filename)
        print(self.filenames)
        self.wf.close()
        return self.filenames

        #assign entries to correct bundles
        #self.bundle.entry = self.entries

        # create a json_str from FHIR BundleResource
        #for x in self.bundles:
        #    json_str = x.json(indent=2)
        #    self.wf.write(json_str)
        #return[str(self.savepath)]
        #json_str = self.bundles.json(indent=2)

        # write json_file
        #self.wf.write(json_str)

        #self.wf.close()
        #return [str(self.savepath)]

    def write_preface(self, level: str):
        if level not in self.levels_to_write:
            return

    def write_table_row(self, row):
        #get samples that have variant(row)
        sample_with_variants = row["tagsampler__samples"].split(":")


        for sample in sample_with_variants:
            
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
            obs_row = Observation(status="final", code=code, subject=self.dict_patient[sample])
            #obs_row.meta = MetaType()
            #obs_row.meta.profile = ["http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/variant"]

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

            conn = sqlite3.connect(self.dbpath)
            curs = conn.cursor()
            curs.execute("SELECT COUNT(*) from variant")
            self.num_rows = curs.fetchone()[0]



            self.counter += 1
            id_maker = (self.str_id) + f"{str(row['base__chrom']) + str(row['base__pos']) + str(row['base__ref_base']) + str(row['base__alt_base'])}"
            id = self.uuid_maker(id_maker + self.str_id)
            uri_maker = Uri(f"urn:uuid:{id}")
            converted_ent = BundleEntry(resource=obs_row, fullUrl=uri_maker)


            print(len(self.dict_entries[sample]))
            self.dict_entries[sample].append(converted_ent)


