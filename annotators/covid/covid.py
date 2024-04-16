import sqlite3
import pandas as pd
from oakvar import BaseAnnotator
import os

class Annotator(BaseAnnotator):
    
    def setup(self):
        # Convert CSV to SQLite
        csv_file_path = "./data/covid_data.csv"
        db_file_path = "./data/covid_variants.sqlite"
        table_name = "OakVar_Covid_Annotations"
        
        # If the SQLite database file already exists, delete it
        if os.path.exists(db_file_path):
            os.remove(db_file_path)
        
        conn = sqlite3.connect(db_file_path)
        df = pd.read_csv(csv_file_path)
        df.to_sql(table_name, conn, if_exists='replace', index=False)
        conn.close()

    def annotate(self, input_data: dict):
        # Path to SQLite file
        db_file_path = "./data/covid_variants.sqlite"
        table_name = "OakVar_Covid_Annotations"

        chrom = input_data["chrom"]
        pos = input_data["pos"]
        alt_base = input_data["alt_base"]

        conn = sqlite3.connect(db_file_path)
        cursor = conn.cursor()

        query = f"SELECT desc FROM {table_name} WHERE chrom=? AND pos=? AND alt_base=?"
        cursor.execute(query, (chrom, pos, alt_base))

        description = cursor.fetchone()
        conn.close()  # Close the connection

        return {"desc": description[0] if description else ""}
