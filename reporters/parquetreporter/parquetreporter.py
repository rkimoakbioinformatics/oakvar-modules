from oakvar import BaseReporter
import sqlite3
import os
import pyarrow as pa
import pyarrow.parquet as pq
import duckdb as db
import zipfile

class Reporter(BaseReporter): 
    def setup(self):
        self.wf = None
        self.filenames = []
        #list of dict -> parquet file 
        self.filename = None
        self.filename_prefix = None
        self.chunkno = 0
        self.chunk_size = 100000
        self.batch = []
        self.rowno = 1
        self.list_tbl = []
        


        
        if self.savepath == None:
            self.filename_prefix = "cravat_result"
        else:
            self.filename_prefix = self.savepath
        self.levels_to_write = self.get_standardized_module_option(
            self.confs.get("pages", "variant")
        )

        #for now use sqlite3 package to get total variant row count and store it for chunking purposes
        sql_conn = sqlite3.connect(self.dbpath)
        curs = sql_conn.cursor()
        curs.execute('SELECT COUNT(*) from "variant"')
        self.results = curs.fetchone()[0]
    


        self.zip = (
            self.get_standardized_module_option(self.confs.get("zip", "false")) == True
        )
        self.module_col_sep = "."
        self.filename_postfix = ".parquet"


    def should_write_level(self,level):
        if self.levels_to_write is None:
            return True
        elif level in self.levels_to_write:
            return True
        else:
            return False
          
    def end(self): 
        if self.wf is not None:
            self.wf.close()
        if self.zip and self.filename_prefix:
            zipfile_path = self.filename_prefix + f"{self.filename_postfix}.zip"
            zf = zipfile.ZipFile(
                zipfile_path, mode="w", compression=zipfile.ZIP_DEFLATED
            )
            for filename in self.filenames:
                zf.write(
                    filename, os.path.relpath(filename, start=os.path.dirname(filename))
                )
            zf.close()
        else:
            zipfile_path = self.filenames
        
        return self.savepath
    


    def write_preface(self,level):

        #determine if current level needs to be written, if not return
        if level not in self.levels_to_write:
            return
        if self.wf is not None:
            self.wf.close()

        

    def write_table_row(self,row):
        self.batch.append(row)
        

        if self.rowno % self.chunk_size == 0 or self.rowno == self.results:
            
            
            #set the file name with correct chunk number
            self.filename = f"{self.filename_prefix}{self.chunkno}{self.filename_postfix}"
            self.filenames.append(self.filename)

            #connect to duckDB            
            conn = db.connect()

            #create a pyarrow table so that schema is set for duckDB table
            tbl = pa.Table.from_pylist(self.batch)
            tbl_name = f"my_tbl{self.chunkno}"

            

            #create a duckDB table from pyarrow table 
            conn.sql(f"CREATE TABLE {tbl_name} AS SELECT * from tbl")
            

            #create parquet file from duckDB table 
            result = conn.execute(f"COPY {tbl_name} to '{self.filename}' (FORMAT 'PARQUET')")

            self.list_tbl.append(tbl_name)

            #reset parameters for next chunk
            conn.close()
            self.batch = []

            self.rowno += 1
            self.chunkno += 1
        else:
            self.rowno += 1