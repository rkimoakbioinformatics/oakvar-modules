from oakvar import BaseReporter
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
        self.batch_size = 10000
        self.batch = []
        self.rowno = 0 

        conn = db.connect()
        conn.execute('INSTALL parquet')
        conn.execute('LOAD parquet')
        
        if self.savepath == None:
            self.filename_prefix = "cravat_result"
        else:
            self.filename_prefix = self.savepath
        self.levels_to_write = self.get_standardized_module_option(
            self.confs.get("pages", "variant")
        )

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
        conn = db.connect()
        columns = list(self.batch[0].keys())

        #create a list of data types for the duckDB table. 
        tbl = pa.Table.from_pylist(self.batch)
        print(tbl)
        conn.sql("CREATE TABLE my_tbl AS SELECT * from tbl")
        conn.sql("INSERT INTO my_tbl SELECT * FROM tbl")
        print(self.filename)
        #conn.execute(f"CREATE TABLE my_table ({', '.join(f'{col} {type(val).__name__.upper()}' for col, val in self.batch[0].items())})")
        #for row in self.batch:
        #   values = ','.join(f"'{val}" if isinstance(val,str) else str(val) for val in row.values())
        #    conn.execute(f"INSERT INTO my_table ({', '.join(row.keys())}) VALUES ({values})")
        result = conn.execute(f"COPY my_tbl to '{self.filename}' (FORMAT 'PARQUET')")

        return self.savepath
    


    def write_preface(self,level):

        #determine if current level needs to be written, if not return
        if level not in self.levels_to_write:
            return
        if self.wf is not None:
            self.wf.close()

        #write file name for specific level:
        
        self.filename = f"{self.filename_prefix}{self.chunkno}{self.filename_postfix}"
        self.filenames.append(self.filename)
        

    def write_table_row(self,row):
        if self.rowno == self.batch_size:
            print(self.batch)
            conn = db.connect()
            conn.execute(f"COPY tbl TO '{self.filename}' (FORMAT 'PARQUET')")
        else:
            self.batch.append(row)
            self.rowno += 1