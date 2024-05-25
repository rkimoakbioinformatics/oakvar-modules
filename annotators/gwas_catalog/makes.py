import sys
import os
import duckdb

def read_tsv(file_path):
    filtered_tsv = []
    with open(file_path, 'r') as file:
        headers = file.readline().strip().split('\t')
        indices = [headers.index(col) for col in ["SNPS", "STRONGEST SNP-RISK ALLELE", "DISEASE/TRAIT", "P-VALUE", "P-VALUE (TEXT)", "OR or BETA", "PUBMEDID"]]
            
        for line in file:
            values = line.strip().split('\t')
            if values[indices[0]].startswith('rs'):
                rsid_values = [values[idx] for idx in indices]
                filtered_tsv.append(rsid_values)
    
    return filtered_tsv

def split_rows(filtered_tsv):
    split_data = []
    for row in filtered_tsv:
        snps = row[0].split(';') if ';' in row[0] else row[0].split('x') if 'x' in row[0] else row[0].split(',')
        alleles = row[1].split(';') if ';' in row[1] else row[1].split('x') if 'x' in row[1] else row[1].split(',')
        for snp, allele in zip(snps, alleles):
            split_row = [snp.strip(), allele.strip(), *row[2:]]
            split_data.append(split_row)
    return split_data

def split_sci_not(number):
    notation_str = str(number)
    coeff, exp = notation_str.split('E')
    coeff = int(coeff)
    exp = int(exp.lstrip('-'))
    exp = min(exp, 255)

    pval_coeff_bits = coeff & 0xFF
    pval_coeff_bits = pval_coeff_bits << 8
    pval_exp_bits = exp & 0xFF
    pval_bits = pval_coeff_bits | pval_exp_bits

    return pval_bits

def pval_column(split_data):
    pval_rows = []
    seen = set()
    for row in split_data:
        pval_bits = split_sci_not(row[3])
        row_tuple = tuple(row)
        if row_tuple not in seen:
            row[3] = pval_bits
            pval_rows.append(row)
            seen.add(row_tuple)
    return pval_rows

def clean_data(pval_rows):
    cleaned_data = []
    for row in pval_rows:
        rsid = row[0].replace('rs', '').strip()
        allele = ''.join(filter(lambda x: x in 'ACGT', row[1].strip().upper()))
        if rsid.isdigit():
            or_beta = row[5]
            if or_beta:
                or_beta = float(or_beta)
            else:
                or_beta = None 
            cleaned_data.append([rsid, allele] + row[2:5] + [or_beta] + row[6:])
    return cleaned_data

def group_data(cleaned_data):
    grouped_data = {}
    for row in cleaned_data:
        snps = row[0]
        if snps in grouped_data:
            grouped_data[snps].append(row[1:])
        else:
            grouped_data[snps] = [row[1:]]
    return grouped_data

def create_tables(connection, bin_index):
    table_name = f"bin_{bin_index}"
    connection.execute(f"""
        CREATE TABLE IF NOT EXISTS {table_name} (
            rsid UBIGINT PRIMARY KEY,
            data STRUCT(allele VARCHAR, trait VARCHAR, pval_bits USMALLINT, pval_info VARCHAR, or_beta REAL, pubmed_id INTEGER)[]
        )
    """)
    return table_name

def find_bin_index(rsid):
    bin_index = rsid % 1000
    return bin_index

def insert_data(connection, grouped_data):
    for rsid, data in grouped_data.items():
        rsid = int(rsid)
        bin_index = find_bin_index(rsid)
        table_name = f"bin_{bin_index}"
        nested_data = []
        for allele, trait, pval_bits, pval_info, or_beta, pubmed_id in data:
            nested_data.append({
                "allele": allele,
                "trait": trait,
                "pval_bits": pval_bits,
                "pval_info": pval_info,
                "or_beta": or_beta,
                "pubmed_id": pubmed_id
            })
        connection.execute(f"INSERT INTO {table_name} VALUES (?, ?)", (rsid, nested_data))

def process_data(file_path):
    filtered_tsv = read_tsv(file_path)
    split_data = split_rows(filtered_tsv)
    pval_rows = pval_column(split_data)
    cleaned_data = clean_data(pval_rows)
    grouped_data = group_data(cleaned_data)
    return grouped_data

def process_and_create_db(file_path, db_path):
    grouped_data = process_data(file_path)
    conn = duckdb.connect(db_path)
    for i in range(1000):
        create_tables(conn, i)
    insert_data(conn, grouped_data)
    conn.close()

if __name__ == "__main__":
    import argparse
    script_dir = os.path.dirname(os.path.realpath(__file__))
    data_dir = os.path.join(script_dir, 'data')
    os.makedirs(data_dir, exist_ok=True)
    default_db_path = os.path.join(data_dir, 'gwas_catalog.duckdb')
    parser = argparse.ArgumentParser(description='Process GWAS catalog data and create a DuckDB database.')
    parser.add_argument('--file', type=str, required=True, help='Path to the GWAS catalog TSV file')

    args = parser.parse_args()

    process_and_create_db(args.file, default_db_path)
