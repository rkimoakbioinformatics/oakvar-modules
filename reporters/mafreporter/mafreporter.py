import csv

from oakvar import BaseReporter


class Reporter(BaseReporter):

    MAF_COLUMN_MAP = {
        "Hugo_Symbol": "base__hugo",
        "Entrez_Gene_Id": "",
        "Center": "null",
        "NCBI_Build": "GRCh38",
        "Chromosome": "base__chrom",
        "Start_Position": "base__pos",
        "End_Position": "base__pos_end",
        "Strand": "+",
        "Variant_Classification": "base__so",
        "Variant_Type": "special_case",
        "Reference_Allele": "ref_base",
        "Tumor_Seq_Allele1": "alt_base",
        "Tumor_Seq_Allele2": "alt_base",
        "dbSNP_RS": "dbsnp__rsid",
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
        "Allele": "alt_base",
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
        "SYMBOL_SOURCE": "HUGO",
        "HGNC_ID": "base__hugo",
        "BIOTYPE": "",
        "CANONICAL": "null",
        "CCDS": "null",
        "ENSP": "null",
        "SWISSPROT": "base__all_mappings",
        "TREMBL": "null",
        "UNIPARC": "null",
        "RefSeq": "base__refseq",
        "SIFT": "sift__prediction (sift__score)",
        "PolyPhen": "",
        "EXON": "base__exonno",
        "INTRON": "null",
        "DOMAINS": "null",
        "GMAF": "thousandgenomes__af",
        "AFR_MAF": "thousandgenomes__afr_af",
        "AMR_MAF": "thousandgenomes__amr_af",
        "ASN_MAF": "thousandgenomes__eas_af + thousandgenomes__sas_af",
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
    PROTECTED_COLS_TO_DELETE = ['vcf_region', 'vcf_info', 'vcf_format', 'vcf_tumor_gt', 'vcf_normal_gt']

    # the columns that need to be cleared (set as empty) in the final file, if the 'somatic' type is selected.
    PROTECTED_COLS_TO_EMPTY = ['Match_Norm_Seq_Allele1', 'Match_Norm_Validation_Allele1',
                               'Match_Norm_Validation_Allele2', 'n_ref_count', 'n_alt_count']

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # defaults to somatic type if no module options are given
        self.maf_file = None
        self.file_writer = None
        self.maf_type = 'somatic'
        self.filename_prefix = None
        self.filename = None
        self.extension = '.maf'
        self.headers = None

    def setup(self):
        self.set_maf_type()
        if self.savepath:
            self.filename_prefix = f'{self.savepath}.{self.maf_type}'
        else:
            self.filename_prefix = self.maf_type

        self.levels_to_write = ['variant']

    def should_write_level(self, level):
        if level != 'variant':
            return False

    def end(self):
        pass

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
            # somatic files are identical to protected files, except that they have the last 6 columns removed
            self.file_writer.writerow(self.headers[:-6])


    def write_table_row(self, row):
        pass

    def set_maf_type(self):
        if 'type' in self.module_options:
            maf_type = self.module_options['type']
            if maf_type not in ['protected', 'somatic']:
                # if module type is given wrongly, defaults to somatic
                self.maf_type = 'somatic'
            else:
                self.maf_type = maf_type
