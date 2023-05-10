from oakvar import BaseReporter
from typing import Union
from typing import Any
from typing import Dict
from typing import List
from fhir.resources.genetics import(
    GeneticVariantReport,
    Meta,
    Reference,
    CodeableConcept,
    Group,
    Sequence
)
import os
import sqlite3
import json 



class Reporter(BaseReporter):
    def __init__(self, *args, **kwargs):
        self.wf = None 
        self.filename_prefix = None 
        self.filename = None

    def set_output_filename(self):
        if self.savepath == None:
            self.filename_prefix = Path("oakvar_result")
        else:
            self.filename_prefix = self.savepath   
    def setup(self):
        self.filenames=[] 
    
    def end(self): 
        if self.wf is not None:
            self.wf.close()
        report.json = json.dumps(report.as_json())
        return self.filename


        
    """Reporter module template."""

    def write_preface(self, level: str):
        """Method to set up output files and write "preface".
        Preface means content before the headers for output columns.
        Opening an output file and saving the file handler should be done
        in this method.

        Args:
            level (str): Level of reporting. For example, "variant" or "gene".
        """
        self.level = level 
        if self.wf:
            self.wf.close()
        
        sql_conn = sqlite3.connect(self.dbpath)
        curs = sql_conn.cursor()
        curs.execute('SELECT ')
        report = GeneticVariantReport()
        report.subject = Reference(reference="test")

    def write_header(self, level: str):
        """Method to write column headers.
        This method is supposed to write column headers.
        The file handler saved in write_preface method is supposed to be
        used here to write column headers.

        Args:
            level (str): Level of reporting. For example, "variant" or "gene".
        """
        self.level  = level
        pass

    def write_table_row(self, row: Union[Dict[str, Any], List[Any]]):
        """Method to write annotation results.

        The file handler saved in write_preface method is supposed to be
        used here to write column headers.

        Variants are sent to this method one by one as `row`. It is either
        a dictionary or a list. If `self.dictrow` is `True`, a dict will be
        passed to this method. If not, a list will. The dict's keys will be
        output column names. If a list is passed, `self.columns` will have
        the column information in the same order as in the list.

        Args:
            row (Union[Dict[str, Any], List[Any]]): annotation result for a variant.
        """
        _ = row
        pass
    
