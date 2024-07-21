import pandas as pd
import numpy as np
import sqlite3 
import re
import os

def main():
    dbgap_gen()


def dbgap_gen():
    db = sqlite3.connect('dbgap.sqlite')
    cur = db.cursor()
    # cur.execute('DROP TABLE "dbgap"')

    df = pd.read_csv("PheGenI_Association_full.tab", sep = "\t")
    df2 = df.loc[:,['Trait', 'SNP rs', 'P-Value','PubMed']]
    df2.fillna(0)
    df2.to_sql('temp', db, index= False)
    try:
        cur.execute("""
                CREATE TABLE 'dbgap' (
                    trait TEXT,
                    rsid INTEGER,
                    p_value REAL,
                    pubmed INTEGER

                )
            """)
        cur.execute("INSERT INTO dbgap (trait, rsid, p_value, pubmed) SELECT * FROM temp")
        db.commit()

        cur.execute("DROP TABLE temp")

        cur.execute("""CREATE TABLE 'trait' (
            
            uid INTEGER,
            trait TEXT, 
            PRIMARY KEY('uid')
        )""")

        cur.execute("""
                INSERT INTO trait (trait) SELECT DISTINCT(trait) FROM dbgap
                    """)
        db.commit()
    finally:
        db.close()
if __name__ == "__main__":
    os.system("/bin/bash -c \"" + 'wget -O PheGenI_Association_full.tab https://www.ncbi.nlm.nih.gov/projects/gap/eqtl/EpiViewBE.cgi?type=dl.tab' + "\"")
    main()