import csv
import re

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
        "dbSNP_Val_Status": "null",
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
        "HGVSp_Short": "",
        "Transcript_ID": "ensembl_regulatory_build__ensr",
        "Exon_Number": "base__exonno",
        "t_depth": "null",
        "t_ref_count": "null",
        "t_alt_count": "null",
        "n_depth": "null",
        "n_ref_count": "null",
        "n_alt_count": "null",
        "all_effects": "",
        "Allele": "base__alt_base",
        "Gene": "ensembl_regulatory_build__ensr",
        "Feature": "ensembl_regulatory_build__ensr",
        "Feature_type": "null",
        "One_Consequence": "base__so",
        "Consequence": "",
        "cDNA_position": "null",
        "CDS_position": "base__cchange",
        "Protein_position": "base__achange",
        "Amino_acids": "",
        "Codons": "null",
        "Existing_variation": "dbsnp__rsid",
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

    # columns which are considered special case since their values need to be rendered separately
    SPECIAL_CASE_COLS = ('Entrez_Gene_Id', 'Variant_Classification', 'Variant_Type', 'Tumor_Seq_Allele1',
                         'Tumor_Seq_Allele2', 'dbSNP_RS', 'Tumor_Sample_Barcode', 'Amino_acids',
                         'SYMBOL_SOURCE', 'SWISSPROT', 'PolyPhen', 'VARIANT_CLASS',
                         'ExAC_AF_Adj', 'FILTER', 'vcf_region', 'vcf_info', 'vcf_format',
                         'vcf_tumor_gt', 'vcf_normal_gt', 'NCBI_Build', 'TRANSCRIPT_STRAND', 'MINIMISED', 'Strand',
                         'SIFT', 'ASN_MAF', 'SIFT', 'Consequence', 'HGVSp_Short', 'BIOTYPE')

    # TODO: find all base__so equivalents if possible
    # MAF Variant Classification as detailed by the MAF Format Specification.
    # Used in column 9. Different from Column 96 'VARIANT_CLASS'
    SO_TO_MAF_VARIANT_CLASSIFICATION = {
        '3_prime_UTR_variant': "3'UTR",
        '2': "3'Flank",
        '5_prime_UTR_variant': "5'UTR",
        '4': "5'Flank",
        '5': 'IGR',
        'intron_variant': 'Intron',  # check if true?
        'frameshift_truncation': 'Frame_Shift_Del',
        'frameshift_elongation': 'Frame_Shift_Ins',
        'inframe_deletion': 'In_Frame_Del',
        'inframe_insertion': 'In_Frame_Ins',
        'missense_variant': 'Missense_Mutation',
        '12': 'Nonsense_Mutation',
        '13': 'Silent',
        'splice_site_variant': 'Splice_Site',  # check if true?
        '15': 'Translation_Start_Site',
        '16': 'Nonstop_Mutation',
        '17': 'RNA',  # which RNA, there are multiple possibilities
        '18': 'Targeted_Region',
        '19': 'De_novo_Start_InFrame',
        '20': 'De_novo_Start_OutOfFrame'
    }

    # MAF Variant Type as detailed by the MAF Format Specification
    MAF_VARIANT_TYPE = ('SNP', 'DNP', 'TNP', 'ONP', 'INS', 'DEL', 'Consolidated')

    # MAF Effects as detailed by the MAF Format Specification
    MAF_ALL_EFFECTS = ('SYMBOL', 'Consequence', 'HGVSp_Short', 'Transcript_ID', 'RefSeq', 'HGVSc', 'IMPACT',
                       'CANONICAL', 'SIFT', 'PolyPhen', 'Strand')

    # base__so to biotype map from gencode
    SO_TO_BIOTYPE = {
        'processed_transcript': 'processed_transcript',
        'transcribed_unprocessed_pseudogene': 'transcribed_unprocessed_pseudogene',
        'unprocessed_pseudogene': 'unprocessed_pseudogene',
        'miRNA': 'miRNA',
        'lnc_RNA': 'lnc_RNA',  # lncRNA in gencode biotype
        'processed_pseudogene': 'processed_pseudogene',
        'snRNA': 'snRNA',
        'transcribed_processed_pseudogene': 'transcribed_processed_pseudogene',
        'retained_intron': 'retained_intron',
        'NMD_transcript_variant': 'NMD_transcript',
        'misc_RNA': 'misc_RNA',  # miscRNA in gencode biotype
        'unconfirmed_transcript': 'TEC',
        'pseudogene': 'pseudogene',
        'transcribed_unitary_pseudogene': 'transcribed_unitary_pseudogene',
        'NSD_transcript': 'NSD_transcript',
        'snoRNA': 'snoRNA',
        'scaRNA': 'scaRNA',
        'pseudogene_rRNA': 'pseudogene_rRNA',  # rRNA_pseudogene in gencode biotype
        'unitary_pseudogene': 'unitary_pseudogene',
        'polymorphic_pseudogene': 'polymorphic_pseudogene',
        'rRNA': 'rRNA',
        'IG_V_pseudogene': 'IG_V_pseudogene',
        'ribozyme': 'ribozyme',
        'sRNA': 'sRNA',
        'TR_V_gene': 'TR_V_gene',
        'TR_V_pseudogene': 'TR_V_pseudogene',
        'TR_D_gene': 'TR_D_gene',
        'TR_J_gene': 'TR_J_gene',
        'TR_C_gene': 'TR_C_gene',
        'TR_J_pseudogene': 'TR_J_pseudogene',
        'IG_C_gene': 'IG_C_gene',
        'IG_C_pseudogene': 'IG_C_pseudogene',
        'IG_J_gene': 'IG_J_gene',
        'IG_J_pseudogene': 'IG_J_pseudogene',
        'IG_D_gene': 'IG_D_gene',
        'IG_V_gene': 'IG_V_gene',
        'IG_pseudogene': 'IG_pseudogene',
        'translated_processed_pseudogene': 'translated_processed_pseudogene',
        'scRNA': 'scRNA',
        'vault_RNA': 'vault_RNA',
        'translated_unprocessed_pseudogene': 'translated_unprocessed_pseudogene',
        'Mt_tRNA': 'Mt_tRNA',
        'Mt_rRNA': 'Mt_rRNA',
        '2kb_downstream_variant': 'protein_coding',
        '2kb_upstream_variant': 'protein_coding',
        '3_prime_UTR_variant': 'protein_coding',
        '5_prime_UTR_variant': 'protein_coding',
        'intron_variant': 'protein_coding',
        'unknown': 'unknown',
        'synonymous_variant': 'protein_coding',
        'start_retained_variant': 'protein_coding',
        'stop_retained_variant': 'protein_coding',
        'missense_variant': 'protein_coding',
        'complex_substitution': 'protein_coding',
        'stop_lost': 'protein_coding',
        'splice_site_variant': 'protein_coding',
        'stop_gained': 'protein_coding',
        'frameshift_truncation': 'protein_coding',
        'frameshift_elongation': 'protein_coding',
        'inframe_insertion': 'protein_coding',
        'inframe_deletion': 'protein_coding',
        'start_lost': 'protein_coding',
        'exon_loss_variant': 'protein_coding',
        'transcript_ablation': 'protein_coding'
    }

    # used in AA minification
    AANUM_TO_AA1 = {
        'Ala': 'A',
        'Cys': 'C',
        'Asp': 'D',
        'Glu': 'E',
        'Phe': 'F',
        'Gly': 'G',
        'His': 'H',
        'Ile': 'I',
        'Lys': 'K',
        'Leu': 'L',
        'Met': 'M',
        'Asn': 'N',
        'Pro': 'P',
        'Gln': 'Q',
        'Arg': 'R',
        'Ser': 'S',
        'Thr': 'T',
        'Val': 'V',
        'Trp': 'W',
        'Tyr': 'Y',
        'Ter': '*',
        'NOA': '_',
        'XAA': '?',
        'NDA': ' ',
    }

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
        # sets the MAF file type property
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
            self.file_writer = csv.DictWriter(self.maf_file, delimiter='\t', fieldnames='')

    def write_header(self, level):

        if self.maf_type == 'protected':
            self.file_writer.fieldnames = self.MAF_COLUMN_MAP.keys()

            self.file_writer.writeheader()
        else:
            # first step to creating Somatic MAF files
            # TODO: after protected maf is done, somatic will require some filtration still
            column_map = self.MAF_COLUMN_MAP

            for col_to_del in self.PROTECTED_COLS_TO_DELETE:
                column_map.pop(col_to_del)

            self.file_writer.fieldnames = column_map.keys()

            self.file_writer.writeheader()

    def write_table_row(self, row):
        new_row = {}
        for col, val in self.MAF_COLUMN_MAP.items():
            # by default in our reporter these are empty
            # they are also empty in a Somatic MAF
            # according to rules these should be populated in Protected MAF
            if col in self.PROTECTED_COLS_TO_EMPTY:
                new_row[col] = ''
            # checks if column is a special case and handles accordingly
            elif col in self.SPECIAL_CASE_COLS:
                proc_value = self.handle_special_cases(row, col)
                new_row[col] = proc_value
            elif val in ['', 'null']:
                new_row[col] = ''
            else:
                new_row[col] = row[val]

        # all_effects is an exception from SPECIAL_CASE_COLS and from the empty string mapping
        # writing all_effects requires values from multiple converted columns.
        new_row['all_effects'] = self.write_all_effects(new_row)

        self.file_writer.writerow(new_row)

    def handle_special_cases(self, row, col):
        """
        Main function to handle rows which do not contain straight-forward data.
        Contains a series of calls to different functions that handle data transformation
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


        # Not found values are capitalized
        # Column 9
        if col == 'Variant_Classification':
            return self.write_variant_classification(row)

        # TODO: Undocumented case: 'Consolidated' option for column Variant_Type
        # Column 10
        if col == 'Variant_Type':
            variant_ref_base = row['base__ref_base']
            variant_alt_base = row['base__alt_base']
            return self.write_variant_type(variant_ref_base, variant_alt_base)

        # TODO: document the case covered by note 3
        # Columns 12, 13
        if col in ['Tumor_Seq_Allele1', 'Tumor_Seq_Allele2']:
            return self.write_tumor_seq_alleles(row, col)

        # TODO: here also used to be the former special case 'dbSNP_Val_Status'
        #   for the moment it is empty since no references for Validation Statuses have been found
        # Column 14
        if col == 'dbSNP_RS':
            return self.write_dbsnp(row)

        # Column 16
        if col == 'Tumor_Sample_Barcode':
            return row['tagsampler__samples']

        # Column 37
        if col == 'HGVSp_Short':
            return self.write_hgvsp_short(row['base__achange'])

        # Column 52
        if col == 'Consequence':
            return self.write_consequence(row)

        # Column 56
        if col == 'Amino_acids':
            return self.write_aa(row)

        # Column 63
        if col == 'SYMBOL_SOURCE':
            return 'HUGO'

        # Column 65
        if col == 'BIOTYPE':

            if row['base__so'] in self.SO_TO_BIOTYPE:
                return self.SO_TO_BIOTYPE[row['base__so']]
            else:
                return ''

        # Column 69
        if col == 'SWISSPROT':
            return self.write_swissprot(row)

        # this is different from what is shown in the variant table because of documentation rules
        # Column 73
        if col == 'SIFT':
            return self.write_sift(row)

        # used the polyphen2 HDIV algorithm data. Multiple sources seem to have the same rules for both algorithms
        # TODO: research more which should be used.
        # Column 74
        if col == 'PolyPhen':
            return self.write_polyphen2(row)

        # TODO: according to documentation this is a combination of 1000G EAS and SAS
        #   but example MAF files don't even contain this column. Leaving empty for the moment.
        # Column 81
        if col == 'ASN_MAF':
            return ''

        # Column 96
        if col == 'VARIANT_CLASS':
            return self.write_variant_class(row)

        # TODO: the same as 'ASN_MAF' - seems to be missing from examples.
        #   furthermore, as in 'ASN_MAF' the column name is different from the documentation.
        # Column 102
        if col == 'ExAC_AF_Adj':
            return ''

        # TODO: this relies on VCF file information. Postponed.
        # Column 111
        if col == 'FILTER':
            return ''

        # TODO: at the moment it just handles 'vcf_region' by mapping:
        #   original_input__pos, original_input__ref_base, original_input__alt_base
        #   all others are skipped
        # Column 122-126
        if col in ['vcf_region', 'vcf_info', 'vcf_format', 'vcf_tumor_gt', 'vcf_normal_gt']:
            return self.write_vcf_data(row, col)

    def write_variant_classification(self, row):
        base_so = row['base__so']
        variant_class = self.SO_TO_MAF_VARIANT_CLASSIFICATION

        if base_so is None:
            return ''

        if base_so in variant_class:
            return variant_class[base_so]
        else:
            if '_' not in base_so:
                return base_so
            else:
                new_base_so = '_'.join([word.capitalize() for word in base_so.split('_')])
                return new_base_so

    @staticmethod
    def write_variant_type(ref_base, alt_base):
        len_ref_base = len(ref_base)
        len_alt_base = len(alt_base)
        variant_type = ''

        # The rules can be consulted at the MAF File Specification Site
        if (len_ref_base == len_alt_base == 1 or len_ref_base == 1 < len_alt_base) and ref_base != alt_base != '-':
            variant_type = 'SNP'
        elif (len_ref_base == len_alt_base == 2 or len_ref_base == 2 < len_alt_base) and ref_base != alt_base != '-':
            variant_type = 'DNP'
        elif (len_ref_base == len_alt_base == 3 or len_ref_base == 3 < len_alt_base) and ref_base != alt_base != '-':
            variant_type = 'TNP'
        elif len_ref_base > 3 and ref_base != alt_base != '-':
            variant_type = 'ONP'
        elif ref_base != '-' and alt_base == '-':
            variant_type = 'DEL'
        elif (len_ref_base < len_alt_base) or (alt_base != '-' and ref_base == '-'):
            variant_type = 'INS'

        return variant_type

    def write_tumor_seq_alleles(self, row, col):
        ref_base = row['base__ref_base']
        alt_base = row['base__alt_base']
        variant_type = self.write_variant_type(ref_base, alt_base)
        zygosity = ''
        genotype = False

        # The rules can be consulted at the MAF File Specification Site
        if ref_base == alt_base and variant_type == 'SNP':
            zygosity = f'{ref_base}{alt_base}'
            genotype = True
        elif ref_base != alt_base and variant_type == 'SNP' and alt_base != '-':
            zygosity = f'{ref_base}{alt_base}'
            genotype = True

        if col == 'Tumor_Seq_Allele1' and genotype:
            return zygosity[1]
        elif col == 'Tumor_Seq_Allele2' and genotype:
            return zygosity[0]
        elif col == 'Tumor_Seq_Allele1' and not genotype and alt_base != '-':
            return alt_base
        elif col == 'Tumor_Seq_Allele2' and not genotype:
            return ''
        else:
            return ''

    @staticmethod
    def write_consequence(row):
        all_mappings = row['base__all_mappings']

        if not all_mappings:
            return ''

        split_collection = all_mappings.split(';')
        split_terms = [so.split(':')[3] for so in split_collection]
        all_split_terms = [so.split(',') if ',' in so else so for so in split_terms]
        all_so = []

        for term in all_split_terms:
            if isinstance(term, list):
                for _ in term:
                    all_so.append(_)
            else:
                all_so.append(term)

        all_so = set(all_so)

        return str(';'.join(all_so))

    # TODO check the documentation rule here: how would I find if it is in different databases to add 'novel' or 'null'?
    @staticmethod
    def write_dbsnp(row):
        val = row['dbsnp__rsid']

        # reverting to null at the moment
        if val is None:
            return 'null'
        else:
            return val

    def write_hgvsp_short(self, row):
        if row is None:
            return ''

        for aa in self.AANUM_TO_AA1:
            while aa in row:
                row = row.replace(aa, self.AANUM_TO_AA1[aa])

        return row

    def write_variant_class(self, row):
        ref_base = row['base__ref_base']
        alt_base = row['base__alt_base']

        if len(ref_base) == 1 and len(alt_base) == 1:
            return 'SNV'

        if len(ref_base) == len(alt_base) and len(ref_base) >= 2 and len(alt_base) >= 2:
            return 'substitution'

        if ref_base != '-' and alt_base == '-':
            return 'deletion'

        if ref_base != '-' and 1 < len(alt_base) and alt_base != '-' and len(ref_base) != len(alt_base):
            return 'indel'

        if ref_base == '-' and alt_base != '-':
            return 'insertion'

        pass

    def write_aa(self, row):
        base_achange = row['base__achange']

        short_aa = self.write_hgvsp_short(base_achange)
        min_aa_re = re.compile('[A-Z*](?![a-z]+)')

        short_aa_only = re.findall(min_aa_re, short_aa)

        return '/'.join(short_aa_only)

    @staticmethod
    def write_swissprot(row):
        all_mappings = row['base__all_mappings']
        base_transcript = row['base__transcript']
        swiss_prot = ''
        if all_mappings is None:
            return ''

        if base_transcript is None:
            return ''

        if base_transcript in all_mappings.split(';')[0]:
            swiss_prot = all_mappings.split(';')[0].split(':')[0]

        return swiss_prot

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

        # The rules can be consulted at the MAF File Specification Site
        if rank >= 0.957:
            prediction_rank = f'{polyphen_values_map[prediction]}({rank})'
        elif 0.453 <= rank <= 0.956:
            prediction_rank = f'{polyphen_values_map[prediction]}({rank})'
        elif rank < 0.453:
            prediction_rank = f'{polyphen_values_map[prediction]}({rank})'

        return prediction_rank

    @staticmethod
    def write_vcf_data(row, col):
        vcf_input_pos = row["original_input__pos"] or ''
        vcf_input_ref = row["original_input__ref_base"] or ''
        vcf_input_alt = row["original_input__alt_base"] or ''

        if col == 'vcf_region':
            return f'{vcf_input_pos}:{vcf_input_ref}:{vcf_input_alt}'
        else:
            return ''


    # Column 46
    def write_all_effects(self, row):
        effects = []
        for effect in self.MAF_ALL_EFFECTS:
            result_row = row[effect]

            if result_row:
                effects.append(result_row)

        all_effects = ";".join(effects)

        if not all_effects:
            return ''
        else:
            all_effects = f'[{str(all_effects)}]'
            return all_effects

    def end(self):
        pass
