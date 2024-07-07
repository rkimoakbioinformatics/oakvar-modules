import sqlite3
import os
import requests
import argparse
from cyvcf2 import VCF

# Function to create new table in sqlite database
def create_table(chr_num):
    cursor.execute(f'''
        CREATE TABLE IF NOT EXISTS {table_name} (
            CHROM TEXT,
            POS INTEGER,
            REF TEXT,
            ALT TEXT,
            AF_ALL REAL,
            AF_AFR REAL,
            AF_AMR REAL,
            AF_ASJ REAL,
            AF_EAS REAL,
            AF_SAS REAL,
            AF_FIN REAL,
            AF_MID REAL,
            AF_NFE REAL
        )
    ''')
    return

# Function to obtain VCF file from gnomad4 database
def get_vcf(chr_num):
    url = f'https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.chr{chr_num}.vcf.bgz'

    with requests.get(url, stream=True) as response:
        response.raise_for_status()
        with open(filename, 'wb') as file:
            for chunk in response.iter_content(chunk_size=8192):
                file.write(chunk)

    vcf = VCF(filename)
    return vcf

# Function to add data to sqlite database
def add_data(chr_num, vcf):
    for var in vcf:
        CHROM = var.CHROM
        POS = var.POS
        REF = var.REF
        ALT = var.ALT[0]
        AF_ALL = var.INFO.get('AF', 0)
        AF_AFR = var.INFO.get('AF_afr', 0)
        AF_AMR = var.INFO.get('AF_amr', 0)
        AF_ASJ = var.INFO.get('AF_asj', 0)
        AF_EAS = var.INFO.get('AF_eas', 0)
        AF_SAS = var.INFO.get('AF_sas', 0)
        AF_FIN = var.INFO.get('AF_fin', 0)
        AF_MID = var.INFO.get('AF_mid', 0)
        AF_NFE = var.INFO.get('AF_nfe', 0)

        values = (CHROM, POS, REF, ALT, AF_ALL, AF_AFR, AF_AMR, AF_ASJ, AF_EAS, AF_SAS, AF_FIN, AF_MID, AF_NFE)    
        placeholders = ','.join(['?'] * len(values))
        table_name = f'CHR{chr_num}'
        filename = f'chr{chr_num}_exome.vcf.bgz'

        cursor.execute(f'''
                       INSERT INTO {table_name} (CHROM, POS, REF, ALT, AF_ALL, AF_AFR, AF_AMR, AF_ASJ, AF_EAS, AF_SAS, AF_FIN, AF_MID, AF_NFE)
                       VALUES ({placeholders})
                       ''', values)
        

parser = argparse.ArgumentParser(description='Get allele frequency data from gnomad4')
parser.add_argument('-n', '--chr_num', type=str, help='Chromosome to be added to sqlite database.')
args = parser.parse_args()
chr_num = args.chr_num

# Open connection to sqlite database
conn = sqlite3.connect('gnomad4.sqlite')
cursor = conn.cursor()

# Global variables
table_name = f'CHR{chr_num}'
filename = f'chr{chr_num}_exome.vcf.bgz'

# Drop table if it already exists
cursor.execute(f"SELECT count(*) FROM sqlite_master WHERE type='table' AND name='{table_name}'")
table_exists = cursor.fetchone()[0]
if table_exists:
    cursor.execute(f"DROP TABLE {table_name}")

# Main
create_table(chr_num)
vcf = get_vcf(chr_num)
add_data(chr_num, vcf)
print(f"Data added for Chromosome {chr_num}\n")

# Commit changes and close connection
conn.commit()
conn.close()

os.remove(filename)

