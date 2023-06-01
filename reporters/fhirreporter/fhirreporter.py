import zipfile
import os
import sqlite3
from pathlib import Path
from oakvar import BaseReporter
from fhir.resources.patient import Patient
from fhir.resources.observation import Observation
from fhir.resources.observation import ObservationComponent
from fhir.resources.humanname import HumanName
from fhir.resources.codeableconcept import  CodeableConcept
from fhir.resources.coding import Coding 
from fhir.resources.reference import Reference
from fhir.resources.bundle import Bundle, BundleEntry
from fhir.resources.fhirtypes import Uri

class Reporter(BaseReporter):
    def setup(self):

        #establish filename with fhir suffix
        self.prefix = self.savepath 
        self.wf = None 
        self.filenames = []
        print(self.savepath)
        if self.savepath == None:
            self.savepath = Path("oakvar_result.json")
        else:
            if self.savepath.suffix != ".json":
                self.savepath = Path(str(self.prefix)+"-fhir.json")
     
        self.wf = open(self.savepath, "w", encoding="utf-8")
        self.levels_to_write = self.confs.get("pages","variant")

        #create FHIR bundle resource
        self.bundle = Bundle(type="batch")

        #get patient name 
        conn = sqlite3.connect(self.dbpath)
        curs = conn.cursor()
        curs.execute("SELECT 'primary_transcript' FROM 'info' ")
        patient_name = curs.fetchone()[0]

        #get number of rows
        curs = conn.cursor()
        curs.execute("SELECT COUNT(*) from variant")
        self.num_rows = curs.fetchone()[0]

        #create Reference resource for PatientResource
        self.subject = Reference()
        self.subject.reference = "Patient"

        #create and fill in PatientResource
        self.patient = Patient()
        name = HumanName()
        name.use = "official"
        name.family = " "
        name.given = [patient_name]
        self.patient.name = [name]

        #create CodingResource for row ObservationResources to Use
        coding = Coding()
        coding.system = "http://loinc.org"
        coding.code = "8480-6"
        self.code = CodeableConcept()
        self.code.coding = [coding]

        self.obs_list = []





        

    def should_write_level(self,level):
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
        
        #create codingType for row Observation
        coding = Coding()
        coding.system = "loinc"
        coding.code = "8480-6"
        code = CodeableConcept()
        code.coding = [coding]

        #Get Alleles from sqlite file  
        ref = row['base__ref_base']
        alt = row['base__alt_base']

        #create Observation Resource for row  
        obs_row = Observation(status="final", code=code, subject=self.subject)

        #Make Component for reference allele 
        coding_ref = Coding()
        coding_ref.system = "loinc"
        coding_ref.code = "69547-8" #always code for reference allele
        code_ref = CodeableConcept()
        code_ref.coding = [coding_ref]
        comp_ref = ObservationComponent(code=code_ref)
        comp_ref.valueString = ref

        #Make Component for (alt)ernate allele
        coding_alt = Coding()
        coding_alt.system = "loinc"
        coding_alt.code = "69551-0"
        code_alt = CodeableConcept()
        code_alt.coding = [coding_alt]
        comp_alt = ObservationComponent(code=code_alt)
        comp_alt.valueString = alt

        #add componenets to row observation
        obs_row.component = [comp_ref, comp_alt]

        #add row observation to list of observations 
        self.obs_list.append(obs_row)

        #check if all rows are done
        if len(self.obs_list) == self.num_rows:
            
            #Turn each observation into a Bundle Entry of type resource
            entries = []
            for ind_obs in self.obs_list:
                full_uri= Uri("urn:oakvar-test123")
                ind_resource = ind_obs
                converted_entry = BundleEntry(fullUrl=full_uri,resource=ind_resource)
                entries.append(converted_entry)
            self.bundle.entry = entries

            #create a json_str from FHIR BundleResource
            json_str = self.bundle.json()
            
            #write json_file
            self.wf.write(json_str)
    