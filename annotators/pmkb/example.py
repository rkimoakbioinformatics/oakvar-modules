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
    # codons = variant_pmkb['achange'][variant_pmkb['achange'].str.contains('codon')].values.tolist()
    variant_pmkb_missense = variant_pmkb[variant_pmkb['Variant'] == "missense"]
    # print(variant_pmkb_missense)
    codons = variant_pmkb_missense['achange'][variant_pmkb_missense['achange'].str.contains('codon')]
    codons = codons.reset_index(drop = True)
    match = re.search(r'(\d.*)[^A-Za-z]',codons[0])
    cod  = re.search(r'codon...\s\d*\smissense$',codons[0])
    # print(cod)
    # print(match)
    for row in range(len(codons)):
        match = re.search(r'(\d.*)[^A-Za-z]',codons[row])
        d = match.groups()[0]
        codons[row] = re.sub(r'codon...\s\d*\smissense$',f'c.{d} missense',codons[row])
        # print(row)
        # print(codons.str.replace()
    # example = re.compile(r'^codon...\s\d*.missense$')
    # print(codons)
    # p = example.search(codons[0])
    # print(p)
    print(codons)
if __name__ == "__main__":
    main()