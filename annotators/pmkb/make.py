import sqlite3
import pandas as pd
import numpy as np
import re

conn = sqlite3.connect('pmkb.sqlite')
variant_pmkb = pd.read_csv('pmkb_variants.csv')
interpretations_pmkb = pd.read_csv('pmkb_interpretations.csv')
def main():
    # data_manipulation()
    variant_parsing_new()
    # interpretations =interpretations_cleaning(interpretations_pmkb)
    # interpretations_manipulation(interpretations)
def interpretations_cleaning(interpretations_pmkb):
    c = conn.cursor()
    # s = "CSF3R T618I|CSF3R any nonsense|CSF3R any frameshift".split("|")
    # l = []
    # for i in range(len(s)):
    #     l.append(" ".join(s[i].split(" ")[1:]))
    # print("|".join(l))
    #Filling null values 

    variants_null = interpretations_pmkb[interpretations_pmkb['Variant(s)'].isnull()]['Gene']
    # print(variants_null.loc[0])
    # c.execute("SELECT achange FROM variants WHERE gene = 'CSF1R'")
    # g = c.fetchone()[0]
    # print(g)
    c.execute("SELECT Gene FROM pmkb_variants")
    rows = c.fetchall()
    all_genes = [list(row) for row in rows]
    # print(all_genes)
    for ind in interpretations_pmkb[interpretations_pmkb['Variant(s)'].isnull()].index:
        l = []
        l.append(variants_null.loc[ind])
        if l in all_genes:
            c.execute("SELECT achange FROM variants WHERE Gene = ? AND (achange LIKE '_codon%' OR achange LIKE '_exon%')", (variants_null.loc[ind],))
            if variants_null.loc[ind] == "AKT1":
                print(l)
                print(c.fetchone()[0])
            variant = c.fetchone()
            if variant:
                interpretations_pmkb.at[ind, 'Variant(s)'] = variant[0]                    
                # print(interpretations_pmkb.loc[ind, 'Gene'])
    # print(interpretations_pmkb['Variant(s)'].isnull().sum())
    interpretations_pmkb = interpretations_pmkb.dropna(subset=["Variant(s)"])
    interpretations_pmkb = interpretations_pmkb.reset_index(drop = True)
    # interpretations_pmkb.to_csv('interpretations.csv')
    # print(interpretations_pmkb)
    return interpretations_pmkb
    
def interpretations_manipulation(inter):
    # print(inter.loc[8,'Variant(s)'].split("|"))
    # print(inter['Variant(s)'].isnull().sum())
    for ind in inter.index:
        variants = inter.loc[ind,'Variant(s)'].split("|")
        l = []
        for j in variants:
            if "_codon:_any" in j or "_exon:_any" in j:
                description  = " ".join(j.split(" "))
                # inter.loc[ind,'Variant(s)'] = description
                inter.loc[ind,'Variant(s)'] = description
                # l.append(description)
                # print(description)
                # print(ind)
            elif "copy number" in j:
                description  = " ".join(j.split(" "))
                inter.loc[ind,"Variant(s)"] = description
            else:
                description = " ".join(j.split(" ")[1:])
                if 'any' not in description and 'copy' not in description and 'rearrangement' not in description and 'exon' not in description and 'codon' not in description and "_codon" not in description and "_exon" not in description:
                    description = "p." + description
                    l.append(description)
                else:
                # print(description)
                    match = re.search(r'(codon|exon)...\s(\d*.*[0-9])*\s(missense|nonsense|frameshift|insertion|deletion|any)$',description)
                    if match:
                        get_v = match.group()
                        get_v = re.search(r'(\d.*)[^A-Za-z]',get_v)
                        if get_v:
                        # print(get_v)
                    # print(match)
                            d = get_v.group().split(',')
                            if 'exon' in description and 'missense' in description:
                                description = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(missense|nonsense)$',f'_exon:{",".join(d).strip()}:missense',description)
                                l.append(description)
                            elif 'codon' in description and 'missense' in description:
                                description = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(missense|nonsense)$',f'_codon:{",".join(d).strip()}:missense',description)
                                l.append(description)
                            #nonsense
                            elif 'exon' in description and 'nonsense' in description:
                                description = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(missense|nonsense)$',f'_exon:{",".join(d).strip()}:nonsense',description)
                                l.append(description)
                            elif 'codon' in description and 'nonsense' in description:
                                description = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(missense|nonsense)$',f'codon:{",".join(d).strip()}:nonsense',description)
                                l.append(description)
                            #Frameshift
                            elif 'exon' in description and 'frameshift' in description:
                                description = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(frameshift)$',f'_exon:{",".join(d).strip()}:frameshift',description)
                                l.append(description)
                            elif 'codon' in description and 'frameshift' in description:
                                description = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(frameshift)$',f'_codon:{",".join(d).strip()}:frameshift',description)
                                l.append(description)
                            #insertion
                            elif 'codon' in description and 'insertion' in description:
                                description = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(insertion)$',f'_codon:{",".join(d).strip()}:insertion',description)
                                l.append(description)
                            elif 'exon' in description and 'insertion' in description:
                                description = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(insertion)$',f'_exon:{",".join(d).strip()}:insertion',description)
                                l.append(description)
                                    
                            #deletion
                            elif 'codon' in description and 'deletion' in description:
                                description = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(deletion)$',f'_codon:{",".join(d).strip()}:deletion',description)
                                l.append(description)
                            elif 'exon' in description and 'deletion' in description:
                                description = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(deletion)$',f'_exon:{",".join(d).strip()}:deletion',description)
                                l.append(description)
                            #any
                            elif 'codon' in description and 'any' in description:
                                description = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(any)$',f'_codon:{",".join(d).strip()}:any',description)
                                l.append(description)
                                # print(description)
                            elif 'exon' in description and 'any' in description:
                                description = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(any)$',f'_exon:{",".join(d).strip()}:any',description)        
                                l.append(description)
                    else:
                        match_any = re.search(r'^any\s.*',description)
                        if match_any:
                            description = "_codon:" + ":".join(match_any.group().split(" "))
                            l.append(description)
                            
                        else:
                            l.append(description)
                inter.loc[ind,'Variant(s)'] = ("|").join(l)
    
    # print(inter)
    # inter['Variant(s)'].replace('p.', np.nan,inplace = True)
    
    # inter = inter[inter['Variant(s)'].notna()]
    # # print(inter)
    # inter.dropna(subset = ['Variant(s)'],inplace = True)
    inter = inter.reset_index(drop = True)
    # print(inter)
    # inter.to_csv("interpretations.csv")

# def interpretations_split():
# CREATE TABLE temporary_interpretations AS
# 	WITH RECURSIVE split(gene_name, variants, rest, tumor_type, tissue_type, interpretations,citations, pmkb_url) AS(
# 		SELECT gene_name, "",variants||"|", tumor_type, tissue_type,interpretations, citations, pmkb_url FROM interpretations
# 		UNION ALL SELECT
# 		gene_name,
# 		substr(rest,0,instr(rest,"|")),
# 		substr(rest, instr(rest,"|")+1),
# 		tumor_type, tissue_type, interpretations, citations, pmkb_url
# 		FROM split WHERE rest !=""	)
# SELECT gene_name, variants, tumor_type, tissue_type, interpretations, citations, pmkb_url FROM split WHERE variants !="" 



def data_manipulation():
    '''
    Manipulates pmkb_variant dataset as converting codon(s) pos so to _codon:pos:so:ref_aa:alt_aa

        '''

    a_change = variant_pmkb.insert(
        loc = 2,
        column = 'achange',
        value = variant_pmkb['Description'].str.split(" ").str[1:]
    )
    # print(variant_pmkb['achange'])
    variant_pmkb.drop('Description',inplace = True, axis = 1)
    # print(variant_pmkb['achange'].info())
    variant_pmkb['achange'] = variant_pmkb['achange'].str.join(" ")
    variant_pmkb_missense = variant_pmkb.query('Variant == "missense"')
    variant_pmkb_missense = variant_pmkb_missense.reset_index(drop = True)
    aa = variant_pmkb_missense['achange'][~variant_pmkb_missense['achange'].str.contains('codon|exon|anymutation')]

    #iterate through specific column 
    for variant in range(len(variant_pmkb.loc[:,"achange"])):
        pos = variant_pmkb.loc[variant,"achange"]
        get_missense  = re.search(r'(codon|exon)...\s(\d*.*[0-9])*\s(missense|nonsense|frameshift|insertion|deletion|any)$',pos)
        # any_mutation = re.search(r'any\smutation', pos)
        any_category = re.search(r'^any\s*.*',pos)
        #Handling any category ex: any mutation, any frameshift
        if any_category:
            m = any_category.group().split(' ')
            pos = re.sub(r'any\s*.*', f'_codon:_any:{m[1]}:_any:_any',pos)
            variant_pmkb.at[variant,'achange'] = pos
        # if any_mutation:
        #     m = any_mutation.group().split(",")
        #     pos = re.sub(r'any\smutation', f'_codon:_any:any',pos)
        #     variant_pmkb.at[variant,'achange'] = pos
        if get_missense:
            get_missense = get_missense.group()
            match = re.search(r'(\d.*)[^A-Za-z]',get_missense)
            if match:
                d = match.group().split(',')
                #misssense / nonsense
                if 'exon' in pos and 'missense' in pos:
                    pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(missense|nonsense)$',f'_exon:{",".join(d).strip()}:missense:_any:_any',pos)
                    variant_pmkb.at[variant, 'achange'] = pos
                elif 'exon' in pos and 'nonsense' in pos:
                    pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(missense|nonsense)$',f'_exon:{",".join(d).strip()}:nonsense:_any:_any',pos)
                    variant_pmkb.at[variant, 'achange'] = pos
                
                elif 'codon' in pos and'missense' in pos:
                    pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(missense|nonsense)$',f'_codon:{",".join(d).strip()}:missense:_any:_any',pos)
                    variant_pmkb.at[variant, 'achange'] = pos

                elif 'codon' in pos and 'nonsense' in pos:
                    pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(missense|nonsense)$',f'_codon:{",".join(d).strip()}:nonsense:_any:_any',pos)
                    variant_pmkb.at[variant, 'achange'] = pos
                #Frameshift
                elif 'codon' in pos and 'frameshift' in pos:
                    pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\sframeshift$',f'_codon:{",".join(d).strip()}:frameshift:_any:_any',pos)
                    variant_pmkb.at[variant, 'achange'] = pos
                elif 'exon' in pos and 'frameshift' in pos:
                    pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\sframeshift$',f'_exon:{",".join(d).strip()}:frameshift:_any:_any',pos)
                    variant_pmkb.at[variant, 'achange'] = pos
                #insertion
                elif 'codon' in pos and 'insertion' in pos:
                    pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\sinsertion$',f'_codon:{",".join(d).strip()}:insertion:_any:_any',pos)
                    variant_pmkb.at[variant, 'achange'] = pos
                elif 'exon' in pos and 'insertion' in pos:
                    pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\sinsertion$',f'_exon:{",".join(d).strip()}:insertion:_any:_any',pos)
                    variant_pmkb.at[variant, 'achange'] = pos
                #deletion
                elif 'codon' in pos and 'deletion' in pos:
                    pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\sdeletion$',f'_codon:{",".join(d).strip()}:deletion:_any:_any',pos)
                    variant_pmkb.at[variant, 'achange'] = pos
                elif 'exon' in pos and 'deletion' in pos:
                    pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\sdeletion$',f'_exon:{",".join(d).strip()}:deletion:_any:_any',pos)
                    variant_pmkb.at[variant, 'achange'] = pos
                #any 
                elif 'exon' in pos and 'any' in pos:
                    pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\sany$',f'_exon:{",".join(d).strip()}:any:_any:_any',pos)
                    variant_pmkb.at[variant, 'achange'] = pos
                elif 'codon' in pos and 'any' in pos:
                    pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\sany$',f'_codon:{",".join(d).strip()}:any:_any:_any',pos)
                    variant_pmkb.at[variant, 'achange'] = pos
                    
        # elif 'any' not in pos and 'copy' not in pos and 'rearrangement' not in pos and 'exon' not in pos and 'codon' not in pos:
        #     pos = "p." + pos
        #     variant_pmkb.at[variant,'achange'] = pos
    variant_pmkb.to_csv('variant.csv')

def variant_parsing_new():
    '''
    change achange to locationkind:pos:so:ref_aa:alt_aa
    '''
    variants_ds = pd.read_csv('variant.csv')
    for variant in variants_ds.index:
        desc = variants_ds.loc[variant,'Variant']
        v = variants_ds.loc[variant,'achange']
        if desc == 'missense':
            if '_codon' not in v and '_exon' not in v:
                ex = re.search(r'([a-zA-Z])([0-9]*_[A-Za-z0-9]*|[0-9]*)([A-Za-z]*.)',v)
                groups = ex.groups()
                f = re.search(r'(\w)(\d+)_(\w)(\d+).*',"".join(groups))
                if f:
                    group_f = f.groups()
                    v = f'_codon:{group_f[1]}:{group_f[3]}:indel:{group_f[0]}X{group_f[2]}:Y'
                    variants_ds.loc[variant,'achange'] = v
                else:
                    v = f'_codon:{groups[1]}:{desc}:{groups[0]}:{groups[2]}'
                    variants_ds.loc[variant,'achange'] = v
        elif desc == 'nonsense':
            if "_codon" not in v and '_exon' not in v:
                a_change_nonsense = re.search(r'(\w)(\d+)(\*)',v)
                groups_nonsense = a_change_nonsense.groups()
                v = f'_codon:{groups_nonsense[1]}:{desc}:{groups_nonsense[0]}:stop'
                variants_ds.loc[variant,'achange'] = v
                # print(variants_ds.loc[variant,'achange'])
        elif desc == 'indel':
            pass
            # print(variants_ds.loc[variant,'achange'])
            
    variants_ds.to_csv('variant2.csv')
                
if __name__ == "__main__":
    main()

    # print(variant_pmkb[:,"Variant"])
    # codons = variant_pmkb['achange'][variant_pmkb['achange'].str.contains('codon')].values.tolist()
    # variant_pmkb_variant = variant_pmkb[variant_pmkb['Variant']]
    # codons = variant_pmkb_missense['achange'][variant_pmkb_missense['achange'].str.contains('codon')]
    # exons = variant_pmkb_missense['achange'][variant_pmkb_missense['achange'].str.contains('exon')]
    # codons = codons.reset_index(drop = True)
    # exons = exons.reset_index(drop = True)
    '''
    Excerpt from data manipulation
    
    codons = variant_pmkb_missense['achange'][variant_pmkb_missense['achange'].str.contains('codon')]
    exons = variant_pmkb_missense['achange'][variant_pmkb_missense['achange'].str.contains('exon')]
    codons = codons.reset_index(drop = True)
    exons = exons.reset_index(drop = True)

    
    
    match = re.search(r'(\d.*)[^A-Za-z]',codons[0])
    cod  = re.search(r'codon...\s\d*\smissense$',codons[0])
    # print(cod)
    # print(match)
    # print(codons)
    match = re.search(r'codon\(s\)\s(\d*.*[0-9])*\smissense$',codons[5]).groups()[0]
    # match = re.search(r'codon\(s\)\s(\d+([^\d]?+))+missense$',codons[5]).groups()[0]
    print(match.split(","))

    print(match)
    # print(re.sub(r'codon...\s[1-9]\s([,1-9])\smissense$',f'c.{s} missense',codons[20]))
    # for row in range(len(codons)):
    #     match = re.search(r'(\d.*)[^A-Za-z]',codons[row])
    #     d = match.group().split(',')
    #     # print(d)
    #     if len(d) == 1:
    #         codons[row] = re.sub(r'[codon|exon]\(s\)\s(\d*.*[0-9])*\smissense$',f'c.{d[0]} missense',codons[row])
    #     elif len(d) > 1:
    #         codons[row] = re.sub(r'[codon|exon]\(s\)\s(\d*.*[0-9])*\smissense$',f'_codon:{",".join(d)}:missense',codons[row])
    #     print(codons[row])
    #         # codons[row] = re.sub(r'codon...\s\d*\smissense$',f'_codon:{d}:missense',codons[row])
        # print(row)
        # print(codons.str.replace()
    # example = re.compile(r'^codon...\s\d*.missense$')
    # print(codons)
    # p = example.search(codons[0])
    # print(p)

            # print(d)
                # if len(d) == 1 and 'exon' in pos:
                #     pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(missense|nonsense)$',f'_exons.{d[0]}',pos)
                #     variant_pmkb.at[variant, 'achange'] = pos
                # elif len(d) == 1 and 'codon' in pos:
                #     pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(missense|nonsense)$',f'c.{d[0]}',pos)
                #     variant_pmkb.at[variant, 'achange'] = pos

    # working on exons
    for row in range(len(exons)):
        match = re.search(r'(\d.*)[^A-Za-z]',exons[row])
        d = match.group().split(',')
        # print(d)
        if len(d) == 1:
            exons[row] = re.sub(r'[codon|exon]\(s\)\s(\d*.*[0-9])*\smissense$',f'_exons:{d[0]}:missense',exons[row])
        elif len(d) > 1:
            exons[row] = re.sub(r'[codon|exon]\(s\)\s(\d*.*[0-9])*\smissense$',f'_exon:{",".join(d).strip()}:missense',exons[row])
        print(exons[row])
    '''