import sqlite3 
import pandas as pd
import re
import os
from oakvar.lib.module.local import get_module_dir, get_module_conf
import oakvar as ov
from refactor import variant_conditionals, variant_sequel_conditional

class File_manipulation():
    def __init__(self):
        #initialize the object as connection to database using API and read as csv and then do data manipulation
        pass
    #Variants dataset manipulation
    def variant_initial_manipulation(self, db_conn, csv_file):
        variant_pmkb = pd.read_csv(csv_file)
        a_change = variant_pmkb.insert(
            loc = 2,
            column = 'achange',
            value = variant_pmkb['Description'].str.split(" ").str[1:]
        )
        # print(variant_pmkb['achange'])
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
            # any_mutation = re.search(r'any\smutation', pos)
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
        variant_pmkb.to_sql(os.path.splitext(csv_file)[0], db_conn, if_exists = 'replace', index = False )

    #variants new parsing
    def variants_sequel_parsing(self, conn, csv_file):
        variants_ds = pd.read_csv(csv_file)
        
        variant_sequel_conditional(variants_ds)
        
        variants_ds.to_sql(os.path.splitext(csv_file)[0], conn, if_exists = 'replace', index= False)

    def final_variants_codon_exon_trans(self, db_conn, csv_file):
        #To be investigated
        gencode_dir = ov.lib.module.local.get_module_dir("gencode")
        gencode_data_dir = gencode_dir / "data"
        gencode_sqlite_path = gencode_data_dir / "gene_43_10000.sqlite"
        conn = sqlite3.connect(gencode_sqlite_path)
        # conn_pmkb = db_conn
        variants_ds = pd.read_csv(csv_file)
        c = conn.cursor()
        chrom_query = """
                    SELECT
                        chrom
                    FROM
                        chroms
                    WHERE 
                        chromid = (
                        SELECT 
                            chromid
                        FROM 
                            transcript
                        WHERE 
                            name LIKE :transcript)
                    """
        transcript_query = """
                        SELECT 
                            tid
                        FROM 
                            transcript
                        WHERE 
                            name LIKE :transcript
                        """
        for row in variants_ds.index:
            if '_exon' in variants_ds.loc[row , 'achange'] and "NOTCH1" not in variants_ds.loc[row,'Gene']:
                data = variants_ds.loc[row,'achange'].split(":")
                # print(data)
                if "," in variants_ds.loc[row,'achange']:
                    exon = variants_ds.loc[row,'achange'].split(":")[1].split(",")
                elif "-" in variants_ds.loc[row,'achange']:
                    exon = variants_ds.loc[row,'achange'].split(":")[1].split("-")
                else:
                    exon = [variants_ds.loc[row,'achange'].split(":")[1]]
                cstart = []
                cend = []
                for exon_no in exon:
                    exon_no = int(exon_no)
                # print(variants_ds.loc[row,'Transcript ID (GRCh37/hg19)'])
                    t = variants_ds.loc[row, 'Transcript ID (GRCh37/hg19)']
                    # print(t)
                    c.execute(chrom_query, {"transcript": t+"%"})
                    chromd = c.fetchmany()[0]
                    c.execute(transcript_query, {"transcript": t+"%"})
                    tid = c.fetchall()[0][0]
                    # print(tid)
                    t_name = f'transcript_frags_{chromd[0]}'
                    # print(t_name)
                    f = f'SELECT cstart FROM {t_name} WHERE tid = {tid} AND kind = 24 AND exonno = {exon_no - 1}'
                    c.execute(f)
                    f_f = c.fetchmany()
                    if f_f != []:
                        cstart.append(f_f[0][0]//3  + 1)
                    c_end = f'SELECT cstart FROM {t_name} WHERE tid = {tid} AND kind = 24 AND exonno = {exon_no}'
                    c.execute(c_end)
                    c_end_fetch = c.fetchmany()
                    if len(c_end_fetch) == 0:
                        c.execute("""
                        SELECT 
                            alen
                        FROM 
                            transcript
                        WHERE 
                            name LIKE :transcript
                        """, {"transcript": t+"%"})
                        alen = c.fetchmany()
                        cend.append(alen[0][0])
                    else:
                        cend.append(c_end_fetch[0][0]//3+1)
                    # print(f_f)
                codon_end = ",".join([str(i) for i in cend])
                codon_no = ",".join([str(i) for i in cstart])
                variants_ds.loc[row,'achange'] = f'_codon:{codon_no}:{codon_end}:{data[2]}:{data[3]}:{data[4]}'
        variants_ds.to_sql(os.path.splitext(csv_file)[0], db_conn, if_exists='replace', index = False)

        '''
        Interpretations dataset methods
        '''

    def interpretations_cleaning(self, conn, csv_file):
        c = conn.cursor()
        interpretations_pmkb = pd.read_csv(csv_file)
        #Filling null values
        variants_null = interpretations_pmkb[interpretations_pmkb['Variant(s)'].isnull()]['Gene']
        c.execute("SELECT Gene FROM pmkb_variants")
        rows = c.fetchall()
        all_genes = [list(row) for row in rows]
        # print(all_genes)
        for ind in interpretations_pmkb[interpretations_pmkb['Variant(s)'].isnull()].index:
            l = []
            l.append(variants_null.loc[ind])
            if l in all_genes:
                c.execute("SELECT achange FROM pmkb_variants WHERE Gene = ? AND (achange LIKE '_codon%' OR achange LIKE '_exon%')", (variants_null.loc[ind],))
                variant = c.fetchone()
                if l[0] == "WHSC1":
                    pass
                if variant:
                    interpretations_pmkb.loc[ind, 'Variant(s)'] = variant[0]
                    # print(interpretations_pmkb.loc[ind, 'Gene'])
        interpretations_pmkb = interpretations_pmkb.dropna(subset=["Variant(s)"])
        interpretations_pmkb = interpretations_pmkb.reset_index(drop = True)
        interpretations_pmkb.to_csv('pmkb_interpretations.csv')
        # interpretations_pmkb.to_sql(os.path.splitext(csv_file)[0], conn, if_exists = 'replace', index = False)

        
    def interpretations_manipulation(self, conn, csv_file):
        inter = pd.read_csv(csv_file, index_col = False)

        #loop through each row
        for ind in inter.index:
            #get the Varaint row
            variants = inter.loc[ind,'Variant(s)'].split("|")
            #put each manipulated variant in a new list
            l = []

            for j in variants:
                # handle any pos firstly
                if "_codon:_any" in j or "_exon:_any" in j:
                    description  = " ".join(j.split(" "))
                    inter.loc[ind,'Variant(s)'] = description
                #handle copy number variant
                elif "copy number" in j :
                    description  = " ".join(j.split(" "))
                    inter.loc[ind,"Variant(s)"] = description

                elif "_codon:" in j or "_exon:" in j:
                    pass
                else:
                    #remove the gene name part and get the variant
                    description = " ".join(j.split(" ")[1:])
                    if 'any' not in description and 'copy' not in description and 'rearrangement' not in description and 'exon' not in description and 'codon' not in description and "_codon" not in description and "_exon" not in description:
                        description = description
                        l.append(description)
                    else:
                        match = re.search(r'(codon|exon)...\s(\d*.*[0-9])*\s(missense|nonsense|frameshift|insertion|deletion|any|indel)$',description)
                        if match:
                            get_v = match.group()
                            get_v = re.search(r'(\d.*)[^A-Za-z]',get_v)
                            if get_v:
                            # print(get_v)
                        # print(match)
                                d = get_v.group().split(',')
                                if 'exon' in description and 'missense' in description:
                                    description = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(missense|nonsense)$',f'_exon:{",".join(d).strip().replace(" ","")}:missense:_any:_any',description)
                                    l.append(description)
                                elif 'codon' in description and 'missense' in description:
                                    description = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(missense|nonsense)$',f'_codon:{",".join(d).strip().replace(" ","")}:missense:_any:_any',description)
                                    l.append(description)
                                #nonsense
                                elif 'exon' in description and 'nonsense' in description:
                                    description = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(missense|nonsense)$',f'_exon:{",".join(d).strip().replace(" ","")}:nonsense:_any:_any',description)
                                    l.append(description)
                                elif 'codon' in description and 'nonsense' in description:
                                    description = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(missense|nonsense)$',f'_codon:{",".join(d).strip().replace(" ","")}:nonsense:_any:_any',description)
                                    l.append(description)
                                #Frameshift
                                elif 'exon' in description and 'frameshift' in description:
                                    description = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(frameshift)$',f'_exon:{",".join(d).strip().replace(" ","")}:frameshift:_any:_any',description)
                                    l.append(description)
                                elif 'codon' in description and 'frameshift' in description:
                                    description = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(frameshift)$',f'_codon:{",".join(d).strip().replace(" ","")}:frameshift:_any:_any',description)
                                    l.append(description)
                                #insertion
                                elif 'codon' in description and 'insertion' in description:
                                    description = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(insertion)$',f'_codon:{",".join(d).strip().replace(" ","")}:insertion:_any:_any',description)
                                    l.append(description)
                                elif 'exon' in description and 'insertion' in description:
                                    description = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(insertion)$',f'_exon:{",".join(d).strip().replace(" ","")}:insertion:_any:_any',description)
                                    l.append(description)

                                #deletion
                                elif 'codon' in description and 'deletion' in description:
                                    description = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(deletion)$',f'_codon:{",".join(d).strip().replace(" ","")}:deletion:_any:_any',description)
                                    l.append(description)
                                elif 'exon' in description and 'deletion' in description:
                                    description = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(deletion)$',f'_exon:{",".join(d).strip().replace(" ", "")}:deletion:_any:_any',description)
                                    l.append(description)
                                #any
                                elif 'codon' in description and 'any' in description:
                                    description = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(any)$',f'_codon:{",".join(d).strip().replace(" ","")}:any:_any:_any',description)
                                    l.append(description)
                                    # print(description)
                                elif 'exon' in description and 'any' in description:
                                    description = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(any)$',f'_exon:{",".join(d).strip().replace(" ","")}:any:_any:_any',description)        
                                    l.append(description)
                        else:
                            match_any = re.search(r'^any\s.*',description)
                            if match_any:
                                m = match_any.group().split(" ")
                                #if any:mutation convert it to _any
                                m[1] = '_any' if m[1] == 'mutation' else m[1]
                                description = "_codon:" + ":".join(m) +":_any:_any"
                                l.append(description)
                            else:
                                l.append(description)
                    #join the variants with | separator
            inter.loc[ind,'Variant(s)'] = ("|").join(l)
        #Drop rows containing copy number 
        inter = inter.drop(index=[row for row in inter.index if 'copy' in inter.loc[row,'Variant(s)']])
        inter = inter.reset_index(drop = True)
        #drop rows containing rearrangement
        inter = inter.drop(index=[row for row in inter.index if 'rearrangement' in inter.loc[row,'Variant(s)']])
        inter = inter.reset_index(drop = True)
        
        inter.to_csv('pmkb_interpretations.csv')
        inter.to_sql(os.path.splitext(csv_file)[0], conn, if_exists = 'replace', index = False)

    def interpretations_split_achange(self,conn):
        # connect_pmkb = sqlite3.connect('pmkb.sqlite')
        c = conn.cursor()
        c.execute('ALTER TABLE pmkb_interpretations RENAME COLUMN Gene TO gene_name')
        c.execute('ALTER TABLE pmkb_interpretations RENAME COLUMN "Variant(s)" TO variants')
        c.execute('ALTER TABLE pmkb_interpretations RENAME COLUMN "Tumor Type(s)" TO tumor_type')
        c.execute('ALTER TABLE pmkb_interpretations RENAME COLUMN "Tissue Type(s)" TO tissue_type')
        c.execute('ALTER TABLE pmkb_interpretations RENAME COLUMN Interpretations TO interpretations')
        c.execute('ALTER TABLE pmkb_interpretations RENAME COLUMN Citations TO citations')
        c.execute('ALTER TABLE pmkb_interpretations RENAME COLUMN "PMKB URL" TO pmkb_url_interpretation')
        
        c.execute("""
            CREATE TABLE temporary_interpretations AS
            WITH RECURSIVE split(gene_name, variants, rest, tumor_type, tissue_type,pmkb_url_interpretation, interpretations,citations) AS(
                SELECT gene_name, "",variants||"|", tumor_type, tissue_type,pmkb_url_interpretation, interpretations, citations FROM pmkb_interpretations
                UNION ALL SELECT
                gene_name,
                substr(rest,0,instr(rest,"|")),
                substr(rest, instr(rest,"|")+1),
                tumor_type, tissue_type, interpretations, citations, pmkb_url_interpretation
                FROM split WHERE rest !=""	)
                SELECT gene_name, variants, tumor_type, tissue_type, interpretations, citations, pmkb_url_interpretation FROM split WHERE variants !="" 
        """)

        c.execute('DROP TABLE pmkb_interpretations')
        c.execute('ALTER TABLE temporary_interpretations RENAME TO pmkb_interpretations')
        final_table = pd.read_sql_query('SELECT * FROM pmkb_interpretations', conn)
        
        final_table.to_csv('pmkb_interpretations.csv')

    def interpretations_final_manipulation(self,conn,csv_file):
        # interpretations_split = pd.read_sql_query('SELECT * FROM pmkb_interpretations', conn)
        interpretations_split = pd.read_csv(csv_file, index_col = False)
        
        for ind in interpretations_split.index:
            desc = interpretations_split.loc[ind, 'variants']
            if "_codon:" not in desc and  "_exon" not in desc:
                #Q61L
                #capturing missense variants only
                missense = re.search(r'([A-Z])(\d+)([A-Z])',desc)
                if missense:
                    #capture the ref_allele, pos, alt_allele
                    missense_group= missense.groups()
                    desc = f'_codon:{missense_group[1]}:{missense_group[1]}:missense:{missense_group[0]}:{missense_group[2]}'
                    interpretations_split.loc[ind,'variants'] = desc
                else:
                    nonsense = re.search(r'([A-Z])(\d+)(\*)',desc)
                    if nonsense:
                        nonsense_group = nonsense.groups()
                        # print(nonsense_group)
                        desc = f'_codon:{nonsense_group[1]}:{nonsense_group[1]}:nonsense:{nonsense_group[0]}:{nonsense_group[2]}'
                        interpretations_split.loc[ind,'variants'] = desc
                    #capturing deletion only
                    deletion = re.search(r'^(\w)(\d+)del$',desc)
                    if deletion:
                        deletion_group = deletion.groups()
                        desc = f'_codon:{deletion_group[1]}:{deletion_group[1]}:deletion:{deletion_group[0]}:-'
                        interpretations_split.loc[ind,'variants'] = desc
                    #capturing deletion from pos to pos
                    deletions = re.search(r'(\w)(\d*)_?(\w)(\d+)(delins|del)(\w*)',desc)
                    if deletions:
                        deletions_group = deletions.groups()
                        if deletions_group[5] == '':
                            desc = f'_codon:{deletions_group[1]}:{deletions_group[3]}:deletion:{deletions_group[0]}X{deletions_group[2]}:-'
                            interpretations_split.loc[ind,'variants'] = desc
                        else:
                            if deletions_group[0] == 'G':
                                pos = "".join(deletions_group[1:4])
                                desc = f'_codon:{pos}:{pos}:indel:{deletions_group[0]}:{deletions_group[5]}'
                                interpretations_split.loc[ind,'variants'] = desc
                            else:
                                # print(deletions_group)
                                desc = f'_codon:{deletions_group[1]}:{deletions_group[3]}:indel:{deletions_group[0]}X{deletions_group[2]}:{deletions_group[5]}'
                                interpretations_split.loc[ind,'variants']  = desc
                    frameshift = re.search(r'(\w)(\d+)(fs)$',desc)
                    if frameshift:
                        frameshift_group = frameshift.groups()
                        desc = f'_codon:{frameshift_group[1]}:{frameshift_group[1]}:frameshift:{frameshift_group[0]}:-'
                        interpretations_split.loc[ind,'variants'] = desc
                    insertion = re.search(r'(\w)(\d*)_?(\w)(\d+)ins(\w*)',desc)
                    if insertion:
                        insertion_group = insertion.groups()
                        desc = f'_codon:{insertion_group[1]}:{insertion_group[3]}:insertion:{insertion_group[0]}:{insertion_group[4]}'
                        interpretations_split.loc[ind,'variants'] = desc
        interpretations_split.to_csv('pmkb_interpretations.csv')
        interpretations_split.to_sql(os.path.splitext(csv_file)[0], conn, if_exists = 'replace', index = False)

#Upload csv files into the database 
os.system("/bin/bash -c \"" + 'wget -O pmkb_interpretations.csv https://pmkb.weill.cornell.edu/therapies/download.csv' + "\"")
os.system("/bin/bash -c \"" + 'wget -O pmkb_variants.csv https://pmkb.weill.cornell.edu/variants/download.csv' + "\"")


def upload_csv(db_conn, file_csv):
    try:
        df = pd.read_csv(file_csv)
        df.to_sql(os.path.splitext(file_csv)[0], db_conn, if_exists = 'replace', index = False)
        print(f'Uploaded {file_csv} in to the database')
    except Exception as e:
        print(f'Could not upload {file_csv} due to {str(e)}')
#Connect database        
db_conn = sqlite3.connect('pmkb.sqlite', timeout = 30)
#list files in the directory
files_csv = os.listdir('.')
#upload pmkb files into the database
for file_csv in files_csv:
    if 'pmkb_' in file_csv:
        upload_csv(db_conn, file_csv)



'''
Variants data set manipulation
'''
variant_cur = db_conn.cursor()
variants_csv = 'pmkb_variants.csv'
#Initialize pmkb variants object
pmkb_variants = File_manipulation()
#Perform initial manipulation
pmkb_variants.variant_initial_manipulation(db_conn, variants_csv)
#Save new manipulated csv file in the working directory
variants_df = pd.read_sql_query('SELECT * FROM pmkb_variants', db_conn)
variants_df.to_csv('pmkb_variants.csv',index = False)
#Initiate a file object
variants_df = File_manipulation()
variants_df.variants_sequel_parsing(db_conn, variants_csv)
variants_final_df = pd.read_sql_query('SELECT * FROM pmkb_variants', db_conn)
variants_final_df.to_csv('pmkb_variants.csv', index = False)
#Initial file object
variants_final_df = File_manipulation()
variants_final_df.final_variants_codon_exon_trans(db_conn, variants_csv)

###########################################################################
'''
Interpretations manipulation
# '''
interpretations_csv = 'pmkb_interpretations.csv'
pmkb_interpretations = File_manipulation()
pmkb_interpretations.interpretations_cleaning(db_conn, interpretations_csv)
pmkb_inter = pd.read_csv(interpretations_csv, index_col = False)
pmkb_inter.to_sql(os.path.splitext(interpretations_csv)[0], db_conn, if_exists = 'replace', index = False)

pmkb_interpretations.interpretations_manipulation(db_conn, interpretations_csv)
pmkb_interpretations.interpretations_split_achange(db_conn)
pmkb_interpretations.interpretations_final_manipulation(db_conn, interpretations_csv)

variant_cur.execute('ALTER TABLE pmkb_variants RENAME COLUMN "PMKB URL" TO pmkb_url_variants')
variant_cur.execute('ALTER TABLE pmkb_variants RENAME COLUMN Gene TO gene')
variant_cur.execute('ALTER TABLE pmkb_variants RENAME TO variant')
variant_cur.execute('ALTER TABLE pmkb_interpretations RENAME TO interpretations_final')