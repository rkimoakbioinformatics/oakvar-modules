import csv

from oakvar import BaseReporter


class Reporter(BaseReporter):
    MAF_COLUMN_MAP = {
        "Hugo_Symbol": "base__hugo",
        "Entrez_Gene_Id": "ncbigene__entrez",
        "Center": "null",
        "NCBI_Build": "GRCh38",
        "Chromosome": "base__chrom",
        "Start_Position": "base__pos",
        "End_Position": "base__pos_end",
        "Strand": "+",
        "Variant_Classification": "base__so",
        "Variant_Type": "",
        "Reference_Allele": "base__ref_base",
        "Tumor_Seq_Allele1": "base__alt_base",
        "Tumor_Seq_Allele2": "base__alt_base",
        "dbSNP_RS": "",
        "dbSNP_Val_Status": "",
        "Tumor_Sample_Barcode": "tagsampler__samples",
        "Matched_Norm_Sample_Barcode": "null",
        "Match_Norm_Seq_Allele1": "null",
        "Match_Norm_Seq_Allele2": "null",
        "Tumor_Validation_Allele1": "null",
        "Tumor_Validation_Allele2": "null",
        "Match_Norm_Validation_Allele1": "null",
        "Match_Norm_Validation_Allele2": "null",
        "Verification_Status": "null",
        "Validation_Status": "null",
        "Mutation_Status": "null",
        "Sequencing_Phase": "null",
        "Sequence_Source": "null",
        "Validation_Method": "null",
        "Score": "null",
        "BAM_File": "null",
        "Sequencer": "null",
        "Tumor_Sample_UUID": "null",
        "Matched_Norm_Sample_UUID": "null",
        "HGVSc": "base__cchange",
        "HGVSp": "base__achange",
        "HGVSp_Short": "base__achange",
        "Transcript_ID": "",
        "Exon_Number": "base__exonno",
        "t_depth": "null",
        "t_ref_count": "null",
        "t_alt_count": "null",
        "n_depth": "null",
        "n_ref_count": "null",
        "n_alt_count": "null",
        "all_effects": "",
        "Allele": "base__alt_base",
        "Gene": "",
        "Feature": "",
        "Feature_type": "null",
        "One_Consequence": "base__so",
        "Consequence": "base__so",
        "cDNA_position": "null",
        "CDS_position": "base__cchange",
        "Protein_position": "base__achange",
        "Amino_acids": "",
        "Codons": "null",
        "Existing_variation": "",
        "ALLELE_NUM": "null",
        "DISTANCE": "null",
        "TRANSCRIPT_STRAND": "+",
        "SYMBOL": "base__hugo",
        "SYMBOL_SOURCE": "",
        "HGNC_ID": "base__hugo",
        "BIOTYPE": "",
        "CANONICAL": "null",
        "CCDS": "null",
        "ENSP": "null",
        "SWISSPROT": "base__all_mappings",
        "TREMBL": "null",
        "UNIPARC": "null",
        "RefSeq": "base__refseq",
        "SIFT": "",
        "PolyPhen": "",
        "EXON": "base__exonno",
        "INTRON": "null",
        "DOMAINS": "null",
        "GMAF": "thousandgenomes__af",
        "AFR_MAF": "thousandgenomes__afr_af",
        "AMR_MAF": "thousandgenomes__amr_af",
        "ASN_MAF": "",
        "EAS_MAF": "thousandgenomes__eas_af",
        "EUR_MAF": "thousandgenomes__eur_af",
        "SAS_MAF": "thousandgenomes__sas_af",
        "AA_MAF": "null",
        "EA_MAF": "null",
        "CLIN_SIG": "clinvar__sig",
        "SOMATIC": "null",
        "PUBMED": "null",
        "MOTIF_NAME": "null",
        "MOTIF_POS": "null",
        "HIGH_INF_POS": "null",
        "MOTIF_SCORE_CHANGE": "null",
        "IMPACT": "null",
        "PICK": "null",
        "VARIANT_CLASS": "base__so",
        "TSL": "null",
        "HGVS_OFFSET": "null",
        "PHENO": "null",
        "MINIMISED": "1",
        "ExAC_AF": "gnomad__af",
        "ExAC_AF_Adj": "",
        "ExAC_AF_AFR": "gnomad__af_afr",
        "ExAC_AF_AMR": "gnomad__af_amr",
        "ExAC_AF_EAS": "gnomad__af_eas",
        "ExAC_AF_FIN": "gnomad__af_fin",
        "ExAC_AF_NFE": "gnomad__af_nfe",
        "ExAC_AF_OTH": "gnomad__af_oth",
        "ExAC_AF_SAS": "gnomad__af_sas",
        "GENE_PHENO": "null",
        "FILTER": "",
        "CONTEXT": "null",
        "src_vcf_id": "null",
        "tumor_bam_uuid": "null",
        "normal_bam_uuid": "null",
        "case_id": "null",
        "GDC_FILTER": "null",
        "COSMIC": "cosmic__cosmic_id",
        "MC3_Overlap": "null",
        "GDC_Validation_Status": "null",
        "GDC_Valid_Somatic": "null",
        "vcf_region": "",
        "vcf_info": "",
        "vcf_format": "",
        "vcf_tumor_gt": "",
        "vcf_normal_gt": ""
    }

    # the columns that need to be skipped/deleted from the final file, if the 'somatic' type is selected.
    PROTECTED_COLS_TO_DELETE = ('vcf_region', 'vcf_info', 'vcf_format', 'vcf_tumor_gt',
                                'vcf_normal_gt', 'GDC_Valid_Somatic')

    # the columns that need to be cleared (set as empty) in the final file, if the 'somatic' type is selected.
    PROTECTED_COLS_TO_EMPTY = ('Match_Norm_Seq_Allele1', 'Match_Norm_Validation_Allele1',
                               'Match_Norm_Validation_Allele2', 'n_ref_count', 'n_alt_count')

    SPECIAL_CASE_COLS = ('Entrez_Gene_Id', 'Variant_Classification', 'Variant_Type', 'Tumor_Seq_Allele1',
                         'Tumor_Seq_Allele2', 'dbSNP_RS', 'dbSNP_Val_Status', 'Tumor_Sample_Barcode', 'Amino_acids',
                         'Existing_variation', 'SYMBOL_SOURCE', 'SWISSPROT', 'PolyPhen', 'VARIANT_CLASS',
                         'ExAC_AF_Adj', 'FILTER', 'COSMIC', 'vcf_region', 'vcf_info', 'vcf_format',
                         'vcf_tumor_gt', 'vcf_normal_gt', 'NCBI_Build', 'TRANSCRIPT_STRAND', 'MINIMISED', 'Strand',
                         'SIFT', 'ASN_MAF', 'all_effects', 'SIFT')

    # TODO: find all base__so equivalents if possible
    SO_TO_MAF_VARIANT_CLASSIFICATION = {
        '3_prime_UTR_variant': "3'UTR",
        '2': "3'Flank",
        '5_prime_UTR_variant': "5'UTR",
        '4': "5'Flank",
        '5': 'IGR',
        'intron_variant': 'Intron',  # check if true?
        '7': 'Frame_Shift_Del',
        '8': 'Frame_Shift_Ins',
        'inframe_deletion': 'In_Frame_Del',
        'inframe_insertion': 'In_Frame_Ins',
        'missense_variant': 'Missense_Mutation',
        '12': 'Nonsense_Mutation',
        '13': 'Silent',
        'splice_site_variant': 'Splice_Site',
        '15': 'Translation_Start_Site',
        '16': 'Nonstop_Mutation',
        '17': 'RNA',
        '18': 'Targeted_Region',
        '19': 'De_novo_Start_InFrame',
        '20': 'De_novo_Start_OutOfFrame'
    }

    MAF_VARIANT_TYPE = ('SNP', 'DNP', 'TNP', 'ONP', 'INS', 'DEL', 'Consolidated')

    MAF_ALL_EFFECTS = ('SYMBOL', 'Consequence', 'HGVSp_Short', 'Transcript_ID', 'RefSeq', 'HGVSc', 'IMPACT',
                       'CANONICAL', 'SIFT', 'PolyPhen', 'Strand')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # defaults to somatic type if no module options are given
        self.dictrow = True
        self.maf_file = None
        self.file_writer = None
        self.maf_type = 'somatic'
        self.filename_prefix = None
        self.filename = None
        # for test purposes it is .tsv - to be removed in the future
        self.extension = '.maf.tsv'
        self.headers = None

    def setup(self):
        self.set_maf_type()
        if self.savepath:
            self.filename_prefix = f'{self.savepath}.{self.maf_type}'
        else:
            self.filename_prefix = self.maf_type

        self.levels_to_write = ['variant']

    def set_maf_type(self):
        if 'type' in self.module_options:
            maf_type = self.module_options['type']
            if maf_type not in ['protected', 'somatic']:
                # if module type is given wrongly, defaults to somatic
                self.maf_type = 'somatic'
            else:
                self.maf_type = maf_type

    def should_write_level(self, level):
        if level != 'variant':
            return False

    def write_preface(self, level):
        self.filename = f'{self.filename_prefix}{self.extension}'

        if self.maf_file:
            self.maf_file.close()
        else:
            self.maf_file = open(self.filename, "w", encoding="utf-8", newline="")
            self.file_writer = csv.writer(self.maf_file, delimiter='\t')

    def write_header(self, level):
        self.headers = list(self.MAF_COLUMN_MAP.keys())

        if self.maf_type == 'protected':
            self.file_writer.writerow(self.headers)
        else:
            # first step to creating Somatic MAF files
            # TODO: after protected maf is done, somatic will require some filtration still
            self.file_writer.writerow(self.headers[:-6])

    def write_table_row(self, row):
        new_row = []
        for col, val in self.MAF_COLUMN_MAP.items():
            if col in self.PROTECTED_COLS_TO_EMPTY:
                new_row.append('')
            # checks if column is a special case and handles accordingly
            elif col in self.SPECIAL_CASE_COLS:
                proc_value = self.handle_special_cases(row, col)
                new_row.append(proc_value)
            elif val in ['', 'null']:
                new_row.append('')
            else:
                new_row.append(row[val])

        self.file_writer.writerow(new_row)

    def handle_special_cases(self, row, col):
        """
        Main function to handle rows which do not contain straight-forward data
        :param row: database row to handle
        :param col: database column
        :return:
        """
        if col == 'Entrez_Gene_Id':
            return row['ncbigene__entrez']

        if col == 'NCBI_Build':
            return 'GRCh38'

        if col in ['Strand', 'TRANSCRIPT_STRAND']:
            return '+'

        if col == 'MINIMISED':
            return '1'

        # TODO: after all is finished, check the notes for
        #   'base__so values without a MAF mapping can be converted by capitalizing'
        if col == 'Variant_Classification':
            return row['base__so']

        if col == 'Variant_Type':
            variant_ref_base = row['base__ref_base']
            variant_alt_base = row['base__alt_base']
            return self.write_variant_type(variant_ref_base, variant_alt_base)

        if col in ['Tumor_Seq_Allele1', 'Tumor_Seq_Allele2']:
            return row['base__alt_base']

        if col in ['dbSNP_RS', 'dbSNP_Val_Status']:
            return self.write_dbsnp(row, col)

        if col == 'Tumor_Sample_Barcode':
            return row['tagsampler__samples']

        if col == 'SYMBOL_SOURCE':
            return 'HUGO'

        # this is different from what is shown in the variant table because of documentation rules
        # TODO: check if this should be like this
        if col == 'SIFT':
            return self.write_sift(row)

        # used the polyphen2 HDIV algorithm data. Multiple sources seem to have the same rules for both algorithms
        # TODO: research more which should be used.
        if col == 'PolyPhen':
            return self.write_polyphen2(row)

        # TODO: seems to depend on values that are handled AFTER this is run - somehow this has to grab
        #   values and join them in their final state
        if col == 'all_effects':
            return self.write_all_effects(row)

    @staticmethod
    def write_variant_type(ref_base, alt_base):
        len_ref_base = len(ref_base)
        len_alt_base = len(alt_base)
        variant_type = ''

        if len_ref_base == len_alt_base == 1 and ref_base != alt_base and ref_base != alt_base != '-':
            variant_type = 'SNP'
        elif len_ref_base == len_alt_base == 2 and ref_base != alt_base and ref_base != alt_base != '-':
            variant_type = 'DNP'
        elif len_ref_base == len_alt_base == 3 and ref_base != alt_base and ref_base != alt_base != '-':
            variant_type = 'TNP'
        elif len_ref_base == len_alt_base > 3 and ref_base != alt_base and ref_base != alt_base != '-':
            variant_type = 'ONP'
        elif (len_ref_base > len_alt_base) or (ref_base != '-' and alt_base == '-'):
            variant_type = 'DEL'
        elif (len_ref_base < len_alt_base) or (alt_base != '-' and ref_base == '-'):
            variant_type = 'INS'

        return variant_type

    # TODO check the documentation rule here: how would I find if it is in different databases to add 'novel' or 'null'?
    @staticmethod
    def write_dbsnp(row, col):
        if col == 'dbSNP_RS':
            val = row['dbsnp__rsid']
            # reverting to null at the moment
            if val is None:
                return 'null'
            else:
                return val

        if col == 'dbSNP_Val_Status':
            return ''

    def write_all_effects(self, row):
        effects = []

        for effect in self.MAF_ALL_EFFECTS:
            if self.MAF_COLUMN_MAP[effect] in row:
                result_row = row[self.MAF_COLUMN_MAP[effect]]

                if result_row is None:
                    effects.append('')
                else:
                    effects.append(result_row)
            else:
                effects.append('')

        all_effects = f'[{str(";".join(effects))}]'

        return all_effects

    @staticmethod
    def write_sift(row):
        prediction = row['sift__prediction']
        score = row['sift__score']
        confidence = row['sift__confidence']
        sift_values_map = {
            'Tolerated': 'tolerated',
            'Damaging': 'deleterious',
            'Low': '_low_confidence',
            'High': ''
        }

        if prediction is None:
            return ''

        prediction_score = f'{sift_values_map[prediction]}{sift_values_map[confidence]}'

        if score is not None:
            prediction_score += f'({score})'

        return prediction_score

    @staticmethod
    def write_polyphen2(row):
        prediction = row['polyphen2__hdiv_pred']

        if row['polyphen2__hdiv_rank'] is not None:
            rank = float(row['polyphen2__hdiv_rank'])
        else:
            rank = ''

        polyphen_values_map = {
            'D': 'probably_damaging',
            'P': 'possibly_damaging',
            'B': 'benign'
        }
        prediction_rank = ''

        if prediction is None:
            return ''

        # TODO: the check seems to be in variant anyway, the conditionals are redundant
        if rank >= 0.957:
            prediction_rank = f'{polyphen_values_map[prediction]}({rank})'
        elif 0.453 <= rank <= 0.956:
            prediction_rank = f'{polyphen_values_map[prediction]}({rank})'
        elif rank < 0.453:
            prediction_rank = f'{polyphen_values_map[prediction]}({rank})'

        return prediction_rank

    def end(self):
        pass
