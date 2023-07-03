import sqlite3
import pandas as pd
import re

conn = sqlite3.connect('pmkb.sqlite')
variant_pmkb = pd.read_csv('pmkb_variants.csv')

def main():
    data_manipulation()

def data_manipulation():
    a_change = variant_pmkb.insert(
        loc = 2,
        column = 'achange',
        value = variant_pmkb['Description'].str.split(" ").str[1:]
    )
    # print(variant_pmkb['achange'])
    variant_pmkb.drop('Description',inplace = True, axis = 1)
    # print(variant_pmkb['achange'].info())
    variant_pmkb['achange'] = variant_pmkb['achange'].str.join(" ")
    #iterate through specific column 
    for variant in range(len(variant_pmkb.loc[:,"achange"])):
        pos = variant_pmkb.loc[variant,"achange"]
        get_missense  = re.search(r'(codon|exon)...\s(\d*.*[0-9])*\s(missense|nonsense)$',pos)
        if get_missense:
            get_missense = get_missense.group()
            match = re.search(r'(\d.*)[^A-Za-z]',get_missense)
            if match:
                d = match.group().split(',')
            # print(d)
                if len(d) == 1 and 'exon' in pos:
                    pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(missense|nonsense)$',f'_exons.{d[0]}',pos)
                    variant_pmkb.at[variant, 'achange'] = pos
                elif len(d) == 1 and 'codon' in pos:
                    pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(missense|nonsense)$',f'c.{d[0]}',pos)
                    variant_pmkb.at[variant, 'achange'] = pos
                if len(d) > 1 and 'exon' in pos:
                    pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(missense|nonsense)$',f'_exon:{",".join(d).strip()}:missense',pos)
                    variant_pmkb.at[variant, 'achange'] = pos
                    print(pos)
                elif len(d) > 1 and 'codon' in pos:
                    pos = re.sub(r'(codon|exon)...\s(\d*.*[0-9])*\s(missense|nonsense)$',f'_codon:{",".join(d).strip()}:missense',pos)
                    variant_pmkb.at[variant, 'achange'] = pos
                    print(variant)
                    print(pos)
    variant_pmkb.to_csv('variant.csv')
    # print(variant_pmkb[:,"Variant"])
    # codons = variant_pmkb['achange'][variant_pmkb['achange'].str.contains('codon')].values.tolist()
    # variant_pmkb_variant = variant_pmkb[variant_pmkb['Variant']]
    # print(variant_pmkb.query("Variant"))
    
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
if __name__ == "__main__":
    main()