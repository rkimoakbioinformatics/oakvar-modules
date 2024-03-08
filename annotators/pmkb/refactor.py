import sqlite3 
import pandas as pd
import numpy as np
import re
import os        

def variant_conditionals(variant_pmkb,pos,ec_identifier,variant,d):
    if ('exon' in pos or 'codon' in pos) and 'missense' in pos:
            pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\smissense$',f'_{ec_identifier}:{",".join(d).replace(" ","")}:missense:_any:_any',pos)
            variant_pmkb.at[variant, 'achange'] = pos
    elif ('exon' in pos or 'codon' in pos) and 'nonsense' in pos:
        pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\snonsense$',f'_{ec_identifier}:{",".join(d).replace(" ","")}:nonsense:_any:_any',pos)
        variant_pmkb.at[variant, 'achange'] = pos
        
    #Frameshift
    elif ('codon' in pos or 'exon' in pos) and 'frameshift' in pos:
        pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\sframeshift$',f'_{ec_identifier}:{",".join(d).replace(" ","")}:frameshift:_any:_any',pos)
        variant_pmkb.at[variant, 'achange'] = pos
        
    #insertion
    elif ('codon' in pos or 'exon' in pos) and 'insertion' in pos:
        pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\sinsertion$',f'_{ec_identifier}:{",".join(d).replace(" ","")}:insertion:_any:_any',pos)
        variant_pmkb.at[variant, 'achange'] = pos
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

def variant_sequel_conditional(variants_ds):

    for variant in variants_ds.index:
        desc = variants_ds.loc[variant,'Variant']
        v = variants_ds.loc[variant,'achange']
        if '_codon' not in v and '_exon' not in v:
            if desc == 'missense':
                ex = re.search(r'([a-zA-Z])([0-9]*_[A-Za-z0-9]*|[0-9]*)([A-Za-z]*.)',v)
                groups = ex.groups()
                f = re.search(r'(\w)(\d+)_?(\w?)(\d+)?(delins|del)(\w)?',"".join(groups))
                if f:
                    group_f = f.groups()
                    if group_f[4] == 'delins':
                        v = f'_codon:{group_f[1]}:{group_f[3]}:indel:{group_f[0]}X{group_f[2]}:Y'
                        variants_ds.loc[variant,'achange'] = v
                    else:
                        #D419del
                        v = f'_codon:{group_f[1]}:{group_f[1]}:{group_f[4]}etion:{group_f[0]}:-'
                        variants_ds.loc[variant, 'achange'] = v
                else:
                    v = f'_codon:{groups[1]}:{groups[1]}:{desc}:{groups[0]}:{groups[2]}'
                    variants_ds.loc[variant,'achange'] = v

            elif desc == 'nonsense':
                a_change_nonsense = re.search(r'(\w)(\d+)(\*)',v)
                groups_nonsense = a_change_nonsense.groups()
                v = f'_codon:{groups_nonsense[1]}:{groups_nonsense[1]}:{desc}:{groups_nonsense[0]}:*'
                variants_ds.loc[variant,'achange'] = v

            elif desc == 'indel':
                if '_' in v and 'delins' in v:
                    a_change_indel = re.search(r'(\w)(\d+)_(\w)(\d+)[a-z]+([A-Z0-9\s]+)',v)
                    groups_indel = a_change_indel.groups()
                    v = f'_codon:{groups_indel[1]}:{groups_indel[3]}:indel:{groups_indel[0]}X{groups_indel[2]}:{groups_indel[4]}'
                    variants_ds.loc[variant,'achange'] = v

                else:
                    a_change_indel = re.search(r'(\w)(\d+|\d+_\w\d+)[a-z]+(\w+|\s)',v)
                    groups_indel = a_change_indel.groups()
                    f = re.search(r'(\w)(\d+)_(\w)(\d+)[a-z]',"".join(groups_indel))
                    if f:
                        group_f = f.groups()
                        v = f'_codon:{group_f[1]}:{group_f[3]}:indel:{group_f[0]}X{group_f[2]}:{group_f[2]}'
                        variants_ds.loc[variant,'achange'] = v
                    else:
                        v = f'_codon:{groups_indel[1]}:indel:{groups_indel[0]}:{groups_indel[2]}'
                        variants_ds.loc[variant,'achange'] = v
            elif desc == 'frameshift':
                a_change_fs = re.search(r'(\w)(\d+)(?:\w+|(\*)?)',v)
                fs_groups = a_change_fs.groups()
                if fs_groups[2]:
                    v = f'_codon:{fs_groups[1]}:{desc}:{fs_groups[0]}:*'
                    variants_ds.loc[variant,'achange'] = v
                else:
                    v = f'_codon:{fs_groups[1]}:{desc}:{fs_groups[0]}:-'
                    variants_ds.loc[variant, 'achange'] = v

            elif desc == 'insertion':
                achange_ins = re.search(r'(\w)(\d+)_(\w)(\d+)ins(\w+)',v)
                groups_ins = achange_ins.groups()
                v = f'_codon:{groups_ins[3]}:{groups_ins[3]}:{desc}:-:{groups_ins[4]}'
                variants_ds.loc[variant, 'achange'] = v
            elif desc == 'deletion':
                achange_del = re.search(r'(\w)(\d+)_?(\w?)(\d+)?del',v)
                groups_del = achange_del.groups()
                if groups_del[3] == None:
                    v = f'_codon:{groups_del[1]}:{desc}:{groups_del[0]}:-'
                    variants_ds.loc[variant,'achange'] = v
                else:
                    v = f'_codon:{groups_del[1]}:{groups_del[3]}:{desc}:{groups_del[0]}X{groups_del[2]}:-'
                    variants_ds.loc[variant,'achange'] = v
                
def variant_initial_manipulation(csv_file):
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


db_conn = sqlite3.connect('pmkb.sqlite')
variant_cur = db_conn.cursor()
variants_csv = 'pmkb_variants.csv'

test_csv = 'test.csv'

# variant_initial_manipulation(variants_csv, db_conn)
# variant_sequel_manipulation(test_csv)
    
