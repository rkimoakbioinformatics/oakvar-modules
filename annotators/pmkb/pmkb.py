from typing import Optional
from typing import List
from oakvar import BaseAnnotator
import re


class Annotator(BaseAnnotator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.query = """ 
                    SELECT 
                        interpretations_final.gene_name, 
                        interpretations_final.tumor_type,
                        interpretations_final.tissue_type,
                        interpretations_final.pmkb_url_interpretations,
                        interpretations_final.interpretations, 
                        interpretations_final.citations,
                        variant.achange, 
                        variant.pmkb_url_variants 
                    FROM 
                        interpretations_final
                    LEFT JOIN
                        variant
                        ON variant.gene = interpretations_final.gene_name
                    WHERE 
                        achange = :achange_pmkb AND variant.gene = :gene
                """
        self.gene_query: str = """
            SELECT 
                achange
            FROM 
                variant
            WHERE 
                gene = :gene 
            """

    def add_to_qr(self, achange_pmkb: str, gene: str, qr):
        if not self.cursor:
            return
        self.cursor.execute(
            f"{self.query}", {"achange_pmkb": achange_pmkb, "gene": gene}
        )
        qr.append(self.cursor.fetchone())

    def annotate_mis(self, achange: str, gene: str, exonno: int, qr, pmkb_variants):
        if not self.cursor:
            return
        # search for the description of missense variant
        input_ref_alt_pos = re.search(r"(\w{3})(\d+)(\w{3}|\?)", achange)
        if not input_ref_alt_pos:
            raise Exception(f"{achange} cannot be parsed.")
        ref_alt_pos_catch = input_ref_alt_pos.groups()
        input_pos = int(ref_alt_pos_catch[1])
        input_ref_allele = seq1(ref_alt_pos_catch[0])
        input_alt_allele = seq1(ref_alt_pos_catch[2])
        # loop through results to get a match
        for row in pmkb_variants:
            achange_pmkb = row[0]
            pmkb_achange = achange_pmkb.split(":")
            # get pos:ref_allele_alt_allele
            pmkb_pos = codon_trans(pmkb_achange[1])
            pmkb_ref_allele = pmkb_achange[-2]
            pmkb_alt_allele = pmkb_achange[-1]
            pmkb_so = pmkb_achange[-3]
            pmkb_kind = pmkb_achange[0]
            # handle _any kind
            if pmkb_so in ["_any", "any"] and pmkb_pos[0] == -1:
                self.add_to_qr(achange_pmkb, gene, qr)
            # match results of pmkb achange with input data
            elif input_pos in pmkb_pos and (
                (
                    pmkb_ref_allele == input_ref_allele
                    and pmkb_alt_allele == input_alt_allele
                )
                or (
                    pmkb_ref_allele == "_any"
                    and pmkb_alt_allele == "_any"
                    and pmkb_kind != "_exon"
                )
                or (pmkb_kind == "_exon" and exonno in pmkb_pos)
            ):
                self.add_to_qr(achange_pmkb, gene, qr)

    def annotate_frameshift(
        self, achange: str, gene: str, exonno: int, qr, pmkb_variants
    ):
        if not self.cursor:
            return
        # get input data for querying the pmkb database
        input_ref_alt_pos = re.search(r"(\w{3})(\d+)", achange)
        if input_ref_alt_pos is None:
            raise Exception(f"{achange} cannot be parsed.")
        ref_alt_pos_catch = input_ref_alt_pos.groups()
        input_pos = ref_alt_pos_catch[1]
        input_ref_allele = seq1(ref_alt_pos_catch[0])
        # get pmkb pos:ref_allele:alt_allele
        for row in pmkb_variants:
            # get pos:ref_allele_alt_allele
            achange_pmkb = row[0]
            pmkb_achange = achange_pmkb.split(":")
            pmkb_pos = codon_trans(pmkb_achange[1])
            pmkb_ref_allele = pmkb_achange[-2]
            pmkb_alt_allele = pmkb_achange[-1]
            pmkb_so = pmkb_achange[-3]
            pmkb_kind = pmkb_achange[0]
            # handle _any kind
            if (pmkb_so == "_any" or pmkb_so == "any") and pmkb_pos[0] == -1:
                self.add_to_qr(achange_pmkb, gene, qr)
            elif input_pos in pmkb_pos and (
                pmkb_ref_allele == input_ref_allele
                or (pmkb_ref_allele == "_any" and pmkb_alt_allele == "_any")
                or (pmkb_kind == "_exon" and exonno in pmkb_pos)
            ):
                self.add_to_qr(achange_pmkb, gene, qr)

    def annotate_inframe_insertion(
        self, achange: str, gene: str, exonno: int, qr, pmkb_variants
    ):
        if not self.cursor:
            return
        input_ref_alt_pos = re.search(r"(\w{3})(\d+)_?(\w{3})(\d+)ins(\w+)", achange)
        if input_ref_alt_pos is None:
            raise Exception(f"{achange} cannot be parsed.")
        ref_alt_pos_catch = input_ref_alt_pos.groups()
        input_start_pos = ref_alt_pos_catch[1]
        input_alt_allele = ref_alt_pos_catch[-1]
        # Query pmkb database
        # loop through results to get a match
        for row in pmkb_variants:
            achange_pmkb = row[0]
            pmkb_achange = achange_pmkb.split(":")
            pmkb_pos = codon_trans(pmkb_achange[1])
            pmkb_ref_allele = pmkb_achange[-2]
            pmkb_alt_allele = pmkb_achange[-1]
            input_alt_allele_seq1 = ""
            pmkb_so = pmkb_achange[-3]
            pmkb_kind = pmkb_achange[0]
            # convert alt_allele to one letter amino acid notation
            if input_alt_allele != None:
                for i in range(0, len(input_alt_allele) + 1, 3):
                    if i < len(input_alt_allele):
                        input_alt_allele_seq1 += seq1(input_alt_allele[i : i + 3])
            # handle _any kind
            if (pmkb_so == "_any" or pmkb_so == "any") and pmkb_pos[0] == -1:
                self.add_to_qr(achange_pmkb, gene, qr)
            elif input_start_pos in pmkb_pos and (
                pmkb_alt_allele == input_alt_allele_seq1
                or (pmkb_ref_allele == "_any" and pmkb_alt_allele == "_any")
                or (pmkb_kind == "_exon" and exonno in pmkb_pos)
            ):
                self.add_to_qr(achange_pmkb, gene, qr)

    def annotate_inframe_deletion(
        self, achange: str, gene: str, exonno: int, qr, pmkb_variants
    ):
        input_ref_alt_pos = re.search(r"(\w{3})(\d+)_?(\w{3})?(\d+)?del", achange)
        if input_ref_alt_pos is None:
            raise Exception(f"{achange} cannot be parsed.")
        ref_alt_pos_catch = input_ref_alt_pos.groups()
        input_ref_allele = seq1(ref_alt_pos_catch[0])
        input_alt_allele = (
            "" if ref_alt_pos_catch[2] == None else seq1(ref_alt_pos_catch[2])
        )
        input_start_pos = ref_alt_pos_catch[1]
        input_end_pos = ref_alt_pos_catch[3]
        input_ref_alt = input_ref_allele + "X" + input_alt_allele
        for row in pmkb_variants:
            achange_pmkb = row[0]
            pmkb_achange = achange_pmkb.split(":")
            pmkb_start_pos = pmkb_achange[1]
            pmkb_ref_allele = pmkb_achange[-2]
            pmkb_alt_allele = pmkb_achange[-1]
            pmkb_end_pos = pmkb_achange[2] if type(pmkb_achange[2]) is int else ""
            pmkb_so = pmkb_achange[-3]
            pmkb_kind = pmkb_achange[0]
            # convert input ref_allele to one letter notation
            # handle _any so
            if (pmkb_so == "_any" or pmkb_so == "any") and pmkb_start_pos == "_any":
                self.add_to_qr(achange_pmkb, gene, qr)
            # if the deletion is in an interval
            elif (
                (
                    (
                        input_end_pos != ""
                        and pmkb_end_pos != ""
                        and input_alt_allele != ""
                    )
                    and (
                        pmkb_start_pos == input_start_pos
                        and pmkb_ref_allele == input_ref_alt
                        and pmkb_end_pos == input_end_pos
                    )
                )
                or (
                    pmkb_ref_allele == input_ref_allele
                    and pmkb_start_pos == input_start_pos
                )
                or (
                    pmkb_ref_allele == "_any"
                    and pmkb_alt_allele == "_any"
                    and pmkb_start_pos == input_start_pos
                )
                or (pmkb_kind == "_exon" and exonno in pmkb_start_pos)
            ):
                self.add_to_qr(achange_pmkb, gene, qr)

    def annotate_css(self, achange: str, gene: str, exonno: int, qr, pmkb_variants):
        if not self.cursor:
            return
        input_ref_alt_pos = re.search(
            r"(\w{3})(\d+)_?(\w{3})?(\d+)?delins(\w+)", achange
        )
        if input_ref_alt_pos is None:
            raise Exception(f"{achange} cannot be parsed.")
        ref_alt_pos_catch = input_ref_alt_pos.groups()
        input_start_pos = ref_alt_pos_catch[1]
        input_end_pos = ref_alt_pos_catch[3] if ref_alt_pos_catch[3] != None else ""
        input_alt_allele = ref_alt_pos_catch[-1]
        for row in pmkb_variants:
            achange_pmkb = row[0]
            pmkb_achange = achange_pmkb.split(":")
            # get pos:ref_allele_alt_allele
            pmkb_start_pos = pmkb_achange[1]
            # get pos:alt_allele
            pmkb_end_pos = pmkb_achange[2]
            pmkb_ref_allele = pmkb_achange[-2]
            pmkb_alt_allele = pmkb_achange[-1]
            input_alt_allele_seq1 = ""
            pmkb_so = pmkb_achange[-3]
            pmkb_kind = pmkb_achange[0]
            # convert alt_allele to one letter amino acid notation
            if input_alt_allele != None:
                for i in range(0, len(input_alt_allele) + 1, 3):
                    if i < len(input_alt_allele):
                        input_alt_allele_seq1 += seq1(input_alt_allele[i : i + 3])
            # handle _any so
            if (pmkb_so == "_any" or pmkb_so == "any"):
                self.add_to_qr(achange_pmkb, gene, qr)
            # if variant description has start and end pos - check matching with the pmkb database
            elif (
                (
                    input_end_pos != ""
                    and pmkb_start_pos == input_start_pos
                    and (
                        pmkb_end_pos == input_end_pos
                        and pmkb_alt_allele == input_alt_allele_seq1
                    )
                )
                or (
                    pmkb_start_pos == input_start_pos
                    and pmkb_alt_allele == input_alt_allele_seq1
                )
                or (pmkb_ref_allele == "_any" and pmkb_alt_allele == "_any")
                or (pmkb_kind == "_exon" and exonno == pmkb_start_pos)
            ):
                self.add_to_qr(achange_pmkb, gene, qr)

    def annotate(self, input_data: dict, secondary_data: Optional[dict] = None):
        assert input_data is not None
        if not self.cursor:
            return None
        # get data for querying the pmkb database from input_data
        gene: str = input_data.get("hugo", "")
        if not gene:
            return None
        self.cursor.execute(self.gene_query, {"gene": gene})
        pmkb_variants = self.cursor.fetchall()
        exonno = input_data.get("exonno", -1)
        if exonno is not None:
            exonno = int(exonno)
        achange: str = input_data.get("achange", "")
        if not achange:
            raise Exception(f"achange is absent in input_data.")
        qr = []
        # exonno = input_data['exonno']
        variant_type = input_data["so"]
        # Handle missense variants
        if variant_type in ["MIS", "missense_variant"]:
            self.annotate_mis(achange, gene, exonno, qr, pmkb_variants)
        # handle frameshift mutations
        elif (
            variant_type == "FSD"
            or variant_type == "FSI"
            or variant_type == "frameshift_insertion"
            or variant_type == "frameshift_deletion"
        ):
            self.annotate_frameshift(achange, gene, exonno, qr, pmkb_variants)
        # handle insertion mutation cases
        elif variant_type in ["INI", "inframe_insertion"]:
            self.annotate_inframe_insertion(achange, gene, exonno, qr, pmkb_variants)
        # Handle deletion category
        elif variant_type in ["IND", "inframe_deletion"]:
            self.annotate_inframe_deletion(achange, gene, exonno, qr, pmkb_variants)
        elif variant_type in ["CSS", "complex_substitution"]:
            self.annotate_css(achange, gene, exonno, qr, pmkb_variants)
            # Handle deletion cases
        if qr and None not in qr:
            return {
                "gene_name": qr[0][0],
                "tumor_type": qr[0][1],
                "tissue_type": qr[0][2],
                "pmkb_url_interpretation": qr[0][3],
                "interpretations": qr[0][4],
                "citations": qr[0][5],
                "achange": qr[0][6],
                "pmkb_url_variants": qr[0][7],
            }
        _ = secondary_data


def codon_trans(seq) -> List[int]:
    if "-" in seq:
        seq = [codon for codon in list(map(int, seq.split("-")))]
        seq = [i for i in range(seq[0], seq[1] + 1)]
        return seq
    if "," in seq:
        seq = [int(codon) for codon in seq.split(",")]
        return seq
    else:
        if seq != "_any":
            return [int(seq)]
        else:
            return [-1]


def seq1(seq, custom_map=None):
    protein_letters_3to1 = {
        "Ala": "A",
        "Cys": "C",
        "Asp": "D",
        "Glu": "E",
        "Phe": "F",
        "Gly": "G",
        "His": "H",
        "Ile": "I",
        "Lys": "K",
        "Leu": "L",
        "Met": "M",
        "Asn": "N",
        "Pro": "P",
        "Gln": "Q",
        "Arg": "R",
        "Ser": "S",
        "Thr": "T",
        "Val": "V",
        "Trp": "W",
        "Tyr": "Y",
    }
    if custom_map is None:
        custom_map = {"Ter": "*"}
    onecode = {k.upper(): v for k, v in protein_letters_3to1.items()}
    onecode.update((k.upper(), v) for k, v in custom_map.items())
    seqlist: List[str] = [seq[3 * i : 3 * (i + 1)] for i in range(len(seq) // 3)]
    return "".join(onecode.get(aa.upper(), "?") for aa in seqlist)
