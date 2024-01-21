import sqlite3 
import pandas as pd
import numpy as np
import re
import os        
def variant_conditionals(variant_pmkb,pos,ec_identifier,variant,d): 
    if ('exon' in pos or 'exon' in pos) and 'missense' in pos:
            pos = re.sub(r'(codon-|exon)...\s(\d*.*[0-9])*\s(missense)$',f'_{ec_identifier}:{",".join(d).replace(" ","")}:missense:_any:_any',pos)
            variant_pmkb.at[variant, 'achange'] = pos
    elif ('exon' in pos or 'exon' in pos) and 'nonsense' in pos:
        pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(nonsense)$',f'_{ec_identifier}:{",".join(d).replace(" ","")}:nonsense:_any:_any',pos)
        variant_pmkb.at[variant, 'achange'] = pos
        
    #Frameshift
    elif ('codon' in pos or 'exon' in pos) and 'frameshift' in pos:
        pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\sframeshift$',f'_{ec_identifier}:{",".join(d).replace(" ","")}:frameshift:_any:_any',pos)
        variant_pmkb.at[variant, 'achange'] = pos
        
    #insertion
    elif ('codon' in pos or 'exon' in pos) and 'insertion' in pos:
        print(pos)
        pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\sinsertion$',f'_{ec_identifier}:{",".join(d).replace(" ","")}:insertion:_any:_any',pos)
        variant_pmkb.at[variant, 'achange'] = pos
        print(pos)
    #deletion
    elif ('codon' in pos or 'exon' in pos) and 'deletion' in pos:
        pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\sdeletion$',f'_{ec_identifier}:{",".join(d).replace(" ","")}:deletion:_any:_any',pos)
        variant_pmkb.at[variant, 'achange'] = pos
    #any 
    elif ('exon' in pos or 'codon' in pos) and 'any' in pos:
        pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\sany$',f'_{ec_identifier}:{",".join(d).replace(" ","")}:_any:_any:_any',pos)
        variant_pmkb.at[variant, 'achange'] = pos
    #indel
    elif 'exon' in pos and 'indel' in pos:
        pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\sindel$',f'_exon:{",".join(d).replace(" ","")}:indel:_any:_any',pos)
        variant_pmkb.at[variant, 'achange'] = pos
def variant_initial_manipulation(csv_file, db_conn):
    variant_pmkb = pd.read_csv(csv_file)
    a_change = variant_pmkb.insert(
    loc = 2,
    column = 'achange',
    value = variant_pmkb['Description'].str.split(" ").str[1:]
    )


    variant_pmkb.drop('Description',inplace = True, axis = 1)
    variant_pmkb.drop(variant_pmkb.index[variant_pmkb['Variant'] == 'CNV'], inplace = True)
    variant_pmkb.drop(variant_pmkb.index[variant_pmkb['Variant'] == 'rearrangement'], inplace = True)
    variant_pmkb = variant_pmkb[variant_pmkb['Amino Acid Change'] != 'Unknown']
    variant_pmkb = variant_pmkb.reset_index(drop= True)
    variant_pmkb['achange'] = variant_pmkb['achange'].str.join(" ")

    #iterate through specific column 
    for variant in range(len(variant_pmkb.loc[:,"achange"])):
        pos = variant_pmkb.loc[variant,"achange"]
        get_variant  = re.search(r'(codon|exon)...\s(\d*.*[0-9])*\s(missense|nonsense|frameshift|insertion|deletion|any|indel)$',pos)
        any_category = re.search(r'^any\s*.*',pos)
        #Handling any category ex: any mutation, any frameshift
        if any_category:
            m = any_category.group().split(' ')
            m[1] = '_any' if m[1] == 'mutation' else m[1]
            pos = re.sub(r'any\s*.*', f'_codon:_any:{m[1]}:_any:_any',pos)
            variant_pmkb.at[variant,'achange'] = pos

        if get_variant:
            get_variant = get_variant.group()
            match = re.search(r'(\d.*)[^A-Za-z]',get_variant)
            if match:
                d = match.group().split(',')
                ec_identifier = pos.split(" ")[0].split("(")[0]
                variant_conditionals(variant_pmkb,pos,ec_identifier,variant,d)
            
    variant_pmkb.to_csv('test.csv')
                #misssense / nonsense
                # if ('exon' in pos or 'exon' in pos) and 'missense' in pos:
                #     pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(missense)$',f'_{ec_identifier}:{",".join(d).replace(" ","")}:missense:_any:_any',pos)
                #     variant_pmkb.at[variant, 'achange'] = pos
                # elif ('exon' in pos or 'exon' in pos) and 'nonsense' in pos:
                #     pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(nonsense)$',f'_{ec_identifier}:{",".join(d).replace(" ","")}:nonsense:_any:_any',pos)
                #     variant_pmkb.at[variant, 'achange'] = pos
                    
                # #Frameshift
                # elif ('codon' in pos or 'exon' in pos) and 'frameshift' in pos:
                #     pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\sframeshift$',f'_{ec_identifier}:{",".join(d).replace(" ","")}:frameshift:_any:_any',pos)
                #     variant_pmkb.at[variant, 'achange'] = pos
                    
                # #insertion
                # elif ('codon' in pos or 'exon' in pos) and 'insertion' in pos:
                #     print(pos)
                #     pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\sinsertion$',f'_{ec_identifier}:{",".join(d).replace(" ","")}:insertion:_any:_any',pos)
                #     variant_pmkb.at[variant, 'achange'] = pos
                #     print(pos)
                # #deletion
                # elif ('codon' in pos or 'exon' in pos) and 'deletion' in pos:
                #     pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\sdeletion$',f'_{ec_identifier}:{",".join(d).replace(" ","")}:deletion:_any:_any',pos)
                #     variant_pmkb.at[variant, 'achange'] = pos
                # #any 
                # elif ('exon' in pos or 'codon' in pos) and 'any' in pos:
                #     pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\sany$',f'_{ec_identifier}:{",".join(d).replace(" ","")}:_any:_any:_any',pos)
                #     variant_pmkb.at[variant, 'achange'] = pos
                # #indel
                # elif 'exon' in pos and 'indel' in pos:
                #     pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\sindel$',f'_exon:{",".join(d).replace(" ","")}:indel:_any:_any',pos)
                #     variant_pmkb.at[variant, 'achange'] = pos
    # variant_pmkb.to_sql(os.path.splitext(csv_file)[0], db_conn, if_exists = 'replace', index = False )

#Upload csv files into the database 
# os.system("/bin/bash -c \"" + 'wget -O pmkb_interpretations.csv https://pmkb.weill.cornell.edu/therapies/download.csv' + "\"")

# os.system("/bin/bash -c \"" + 'wget -O pmkb_variants.csv https://pmkb.weill.cornell.edu/variants/download.csv' + "\"")


# def upload_csv(db_conn, file_csv):
#     try:
#         df = pd.read_csv(file_csv)
#         df.to_sql(os.path.splitext(file_csv)[0], db_conn, if_exists = 'replace', index = False)
#         print(f'Uploaded {file_csv} in to the database')
#     except Exception as e:
#         print(f'Could not upload {file_csv} due to {str(e)}')
# #Connect database        
# db_conn = sqlite3.connect('pmkb.sqlite', timeout = 30)
# #list files in the directory
# files_csv = os.listdir('.')
# #upload pmkb files into the database
# for file_csv in files_csv:
#     if 'pmkb' in file_csv:
#         upload_csv(db_conn, file_csv)


db_conn = sqlite3.connect('pmkb.sqlite')
variant_cur = db_conn.cursor()
variants_csv = 'pmkb_variants.csv'


variant_initial_manipulation(variants_csv, db_conn)

    
