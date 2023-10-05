from typing import Any
from typing import Optional
from typing import List
from oakvar import BaseConverter
import re
from collections import defaultdict
from io import StringIO
import copy
from pathlib import Path
from math import isnan
from collections import OrderedDict


class Converter(BaseConverter):

    len_NC001807 = 16571
    len_NC012920 = 16569
    hg38_code = "hg38"
    hg19_code = "hg19"
    hg18_code = "hg18"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.format_name = "vcf"
        self._in_header = True
        self._first_variant = True
        self._buffer = StringIO()
        self._reader = None
        self.ex_info_writer = None
        self.curvar = None
        self.csq_fields = None

    def check_format(self, input_path):
        if input_path.endswith(".vcf"):
            return True
        if input_path.endswith(".vcf.gz"):
            return True
        with open(input_path) as f:
            first_line = f.readline()
            if first_line.startswith("##fileformat=VCF"):
                return True
        return False

    def detect_genome_assembly_from_str(self, s):
        if "GRCh37" in s or "hg19" in s:
            return "hg19"
        if "NCBI36" in s or "hg18" in s:
            return "hg18"
        return None

    def detect_genome_assembly(self, reader, input_path: str):
        genome_assembly = self.detect_genome_assembly_from_metadata(reader)
        if genome_assembly:
            return genome_assembly
        return self.detect_genome_assembly_from_contigs(input_path)

    def get_value_from_contig_line(self, value, line):
        if f"{value}=" not in line:
            return None
        return line.split("assembly=")[1].strip().strip(">").split(",")[0].lower()

    def get_assembly_from_contig_line(self, line):
        return self.get_value_from_contig_line("assembly", line)

    def get_length_from_contig_line(self, line):
        return self.get_value_from_contig_line("length", line)

    def detect_genome_assembly_from_contigs(self, input_path: str):
        import gzip

        if input_path.endswith(".vcf"):
            f = open(input_path)
        elif input_path.endswith(".vcf.gz"):
            f = gzip.open(input_path, "rt")
        else:
            return
        f.seek(0)
        for line in f:
            if line.startswith("##contig="):
                assembly = self.get_assembly_from_contig_line(line)
                if assembly == "b37":
                    return self.hg19_code
                elif assembly == "b36":
                    return self.hg18_code
            if line.startswith("#CHROM"):
                break
        return None

    def detect_genome_assembly_from_metadata(self, reader):
        from collections import OrderedDict

        reference = reader.metadata.get("reference")
        if reference:
            assembly = self.detect_genome_assembly_from_str(reference)
            if assembly:
                return assembly
        for k in reader.metadata.keys():
            v = reader.metadata[k]
            if type(v) == list:
                for vv in v:
                    if type(vv) == str:
                        assembly = self.detect_genome_assembly_from_str(vv)
                        if assembly:
                            return assembly
                    elif type(vv) == OrderedDict:
                        for _, vvv in vv.items():
                            assembly = self.detect_genome_assembly_from_str(vvv)
                            if assembly:
                                return assembly
            elif type(v) == str:
                assembly = self.detect_genome_assembly_from_str(v)
                if assembly:
                    return assembly
        return None

    def setup(self, input_path: str, encoding: str="utf-8"):
        import vcf
        import gzip

        _ = encoding
        if hasattr(self, "conf") == False:
            self.conf = {}
        if type(self.conf.get("exclude_info")) == str:
            self.exclude_info = set(self.conf["exclude_info"].split(","))
        else:
            self.exclude_info = set()
        if type(self.conf.get("include_info")) == str:
            self.include_info = set(self.conf["include_info"].split(","))
        else:
            self.include_info = set()
        if input_path.endswith(".vcf"):
            f = open(input_path)
        elif input_path.endswith(".vcf.gz"):
            f = gzip.open(input_path)
        else:
            return
        reader = vcf.Reader(f, compressed=False)
        self.open_extra_info(reader)
        self.input_assembly = self.detect_genome_assembly(reader, input_path)

    def get_do_liftover_chrM(self, genome_assembly, input_path: str, do_liftover):
        return self.chrM_needs_liftover(genome_assembly, input_path, do_liftover)

    def chrM_needs_liftover(self, genome_assembly, input_path: str, do_liftover):
        import vcf
        import gzip

        if genome_assembly != "hg19":
            return do_liftover
        if input_path.endswith(".vcf"):
            f = open(input_path)
        elif input_path.endswith(".vcf.gz"):
            f = gzip.open(input_path)
        else:
            return
        reader = vcf.Reader(f)
        if not reader.contigs:
            return do_liftover
        for contig in reader.contigs.values():
            if contig.id in ["M", "MT", "Mt"]:
                if contig.length == 16571:
                    return True
                elif contig.length == 16569:
                    f.seek(0)
        return do_liftover

    def open_extra_info(self, reader):
        try:
            from oakvar.lib.util.inout import FileWriter  # type: ignore
        except:
            from oakvar.util.inout import FileWriter  # type: ignore
        assert self.input_paths is not None
        if not self.output_dir or not self.run_name:
            return
        writer_path = Path(self.output_dir) / (self.run_name + ".extra_vcf_info.var")
        if self.input_path == self.input_paths[0]:
            self.mode = "w"
        else:
            self.mode = "a"
        self.ex_info_writer = FileWriter(str(writer_path), mode=self.mode)
        info_cols: list[dict[str, Any]] = [
            {"name": "uid", "title": "UID", "type": "int"}
        ]
        info_cols.append(
            {
                "name": "pos",
                "title": "VCF Position",
                "desc": "",
                "type": "int",
                "width": 60,
            }
        )
        info_cols.append(
            {
                "name": "ref",
                "title": "VCF Ref Allele",
                "desc": "",
                "type": "string",
                "width": 60,
            }
        )
        info_cols.append(
            {
                "name": "alt",
                "title": "VCF Alt Allele",
                "desc": "",
                "type": "string",
                "width": 60,
            }
        )
        typemap = {"Integer": "int", "Float": "float"}
        if reader.infos:
            for info in reader.infos.values():
                # Ensure no duplicate column names exist (case-insensitive)
                if info.id.lower() in [x["name"].lower() for x in info_cols]:
                    info_id = info.id + "_"
                else:
                    info_id = info.id
                info_cols.append(
                    {
                        "name": info_id,
                        "title": info_id,
                        "desc": info.desc,
                        "type": typemap.get(info.type, "string"),
                        "hidden": True,
                    }
                )
            if "CSQ" in reader.infos:
                csq_info = reader.infos["CSQ"]
                fields_match = re.search(r"Format: ([^\s]+)", csq_info.desc)
                if fields_match:
                    self.csq_fields = [
                        "CSQ_" + x for x in fields_match.group(1).split("|")
                    ]
                    for cname in self.csq_fields:
                        info_cols.append(
                            {
                                "name": cname,
                                "title": cname.replace("_", " "),
                                "type": "string",
                                "hidden": True,
                            }
                        )
        if self.include_info:
            self.include_info.update([c["name"] for c in info_cols[:3]])
            temp = info_cols
            info_cols = []
            for col in temp:
                if (
                    col["name"] in self.include_info
                    and col["name"] not in self.exclude_info
                ):
                    col["hidden"] = False
                    info_cols.append(col)
            del temp
        else:
            temp = info_cols
            info_cols = []
            for col in temp:
                if col["name"] not in self.exclude_info:
                    col["hidden"] = False
                    info_cols.append(col)
            del temp
        self.info_cols = [v["name"] for v in info_cols]
        self.ex_info_writer.add_columns(info_cols)
        if self.mode != "a":
            self.ex_info_writer.write_definition()
            self.ex_info_writer.write_meta_line("name", "extra_vcf_info")
            self.ex_info_writer.write_meta_line(
                "displayname", "Extra VCF INFO Annotations"
            )

    def convert_line(self, l):
        import vcf

        try:
            from oakvar.lib.exceptions import NoAlternateAllele  # type: ignore
        except:
            from oakvar.exceptions import NoAlternateAllele  # type: ignore

        if l.startswith("#"):
            if self._in_header:
                self._buffer.write(f"{l}\n")
            return self.IGNORE
        if self._first_variant:
            self._first_variant = False
            self._in_header = False
            self._buffer.seek(0)
            self._reader = vcf.Reader(self._buffer)
        if not self._reader:
            return
        self._buffer.seek(0)
        self._buffer.truncate()
        self._buffer.write(f"{l}\n")
        self._buffer.seek(0)
        try:
            variant = next(self._reader)
        except StopIteration:
            return self.IGNORE
        wdict_blanks = {}
        wdicts = []
        for alt_index, alt in enumerate(variant.ALT):
            if alt is None:
                alt_base = variant.REF
            elif alt.type == "NON_REF":
                alt_base = None
            elif alt.type == "*":
                alt_base = "*"
            else:
                alt_base = alt.sequence
            new_pos, new_ref, new_alt = self.trim_variant(
                variant.POS, variant.REF, alt_base
            )
            if variant.FILTER is None:
                filter_val = None
            elif len(variant.FILTER) == 0:
                filter_val = "PASS"
            else:
                filter_val = ";".join(variant.FILTER)
            wdict_blanks[alt_index + 1] = {
                "chrom": variant.CHROM,
                "pos": new_pos,
                "ref_base": new_ref,
                "alt_base": new_alt,
                "tags": variant.ID,
                "phred": variant.QUAL,
                "filter": filter_val,
                "var_no": alt_index,
            }
        self.gt_occur: List[int] = []
        if len(variant.samples) > 0:
            all_gt_zero = True
            for call in variant.samples:
                # Dedup gt but maintain order
                for gt in list(OrderedDict.fromkeys(call.gt_alleles)):
                    if gt in [None, "0", "."]:
                        continue
                    all_gt_zero = False
                    gt = int(gt)
                    wdict = copy.copy(wdict_blanks[gt])
                    if wdict["alt_base"] == "*":
                        continue
                    wdict["sample_id"] = call.sample
                    if call.is_het == True:
                        wdict["zygosity"] = "het"
                    elif call.is_het == False:
                        wdict["zygosity"] = "hom"
                    elif call.is_het is None:
                        wdict["zygosity"] = None
                    else:
                        wdict["zygosity"] = None
                    (
                        wdict["tot_reads"],
                        wdict["alt_reads"],
                        wdict["af"],
                    ) = self.extract_read_info(variant, call, gt)
                    wdict["hap_block"] = None  # TODO
                    wdict["hap_strand"] = None  # TODO
                    wdict["genotype"] = variant.genotype(call.sample).gt_bases
                    wdicts.append(wdict)
                    self.gt_occur.append(gt)
            if all_gt_zero:
                raise NoAlternateAllele()
        else:
            for gt in wdict_blanks:
                wdict = copy.copy(wdict_blanks[gt])
                if wdict["alt_base"] == "*":
                    continue
                wdicts.append(wdict)
                self.gt_occur.append(gt)
        self.curvar = variant
        self.cur_csq = {}
        if self.csq_fields and "CSQ" in variant.INFO:
            csq_entries = defaultdict(list)
            for gt_csq in variant.INFO["CSQ"]:
                l = gt_csq.split("|")
                csq_entries[l[0]].append(l)
            for allele, entries in csq_entries.items():
                transpose = zip(*entries)
                self.cur_csq[allele] = dict(
                    [
                        (cname, self.csq_format(value))
                        for cname, value in zip(self.csq_fields, transpose)
                    ]
                )
        return wdicts

    @staticmethod
    def extract_read_info(call, variant, sample_num: int, gt):
        tot_reads: Optional[int] = None
        alt_reads: Optional[int] = None
        # AD is depth for each allele
        if hasattr(call.data, "AD"):
            # tot_reads
            if hasattr(call.data.AD, "__iter__"):
                tot_reads = sum([0 if x is None else int(x) for x in call.data.AD])
            elif call.data.AD is None:
                tot_reads = 0
            else:
                tot_reads = int(call.data.AD)
            # alt_reads
            if call.data.AD:
                try:
                    alt_reads = int(call.data.AD[gt])
                except IndexError:  # Wrong length
                    alt_reads = None
                except TypeError:  # Not indexable
                    alt_reads = int(call.data.AD)
            else:
                alt_reads = None
        else:
            if hasattr(call.data, "SOMATIC"):
                alt_reads = 0
                for alt in variant.ALT:  # Collect Strelka reads from AU, CU, GU, and TU
                    ualt: str = alt + "U"
                    if hasattr(call.data, ualt):
                        alt_strelka = getattr(call.data, ualt)
                        tier_1 = alt_strelka[sample_num][0]  # Tier 1 reads is the first index
                        alt_reads += tier_1
            else:
                alt_reads = variant.gt_alt_depths[sample_num]
        # DP is total depth
        if hasattr(call.data, "DP"):
            if call.data.DP:
                variant_dp = call.data.DP
                tot_reads = variant_dp[sample_num][0]
            else:
                tot_reads = None
        if (
            tot_reads not in [-1, None]
            and alt_reads not in [-1, None]
            and tot_reads
            and alt_reads is not None
        ):
            try:
                alt_freq = alt_reads / tot_reads
            except ZeroDivisionError:
                alt_freq = None
        else:
            alt_freq = None
        return tot_reads, alt_reads, alt_freq

    @staticmethod
    def csq_format(l):
        # Format a list of CSQ values into it's representation in OC
        # Each value comes from a VEP transcript mapping
        if all([x == "" for x in l]):
            return None
        else:
            return ";".join(l)

    def trim_variant(self, pos, ref, alt):
        if alt is None:
            return pos, ref, alt
        if len(ref) == 1 and len(alt) == 1:
            return pos, ref, alt
        ref = list(ref)
        alt = list(alt)
        adj = 0
        while ref and alt and ref[0] == alt[0]:
            adj += 1
            ref.pop(0)
            alt.pop(0)
        while ref and alt and ref[-1] == alt[-1]:
            ref.pop()
            alt.pop()
        ref = "".join(ref) if ref else "-"
        alt = "".join(alt) if alt else "-"
        return pos + adj, ref, alt

    @staticmethod
    def oc_info_val(info_type, val, force_str=False):
        if val is None or val == ".":
            oc_val = None
        if info_type in ("Integer", "Float"):
            if val is None:
                oc_val = None
            elif isnan(val):
                oc_val = None
            else:
                oc_val = val
        else:
            oc_val = val
        if force_str and oc_val is None:
            return "."
        elif force_str:
            return str(oc_val)
        else:
            return oc_val

    def addl_operation_for_unique_variant(self, wdict, wdict_no):
        if self.ex_info_writer is None:
            return
        gt = self.gt_occur[wdict_no]
        alt_index = gt - 1
        row_data = {"uid": wdict["uid"]}
        if not self.curvar:
            return
        somatic = True if self.curvar.INFO.get("SOMATIC") else False
        for info_name, info_val in self.curvar.INFO.items():
            if info_name not in self.info_cols:
                continue
            if not self._reader or not self._reader.infos:
                continue
            info_desc = self._reader.infos[info_name]
            if info_desc.num == 0:
                oc_val = self.oc_info_val(info_desc.type, info_val)
            elif info_desc.num in [-1, "A"]:  # Number=A
                val = info_val[alt_index] if hasattr(info_val, "__iter__") else info_val  # If not a list
                oc_val = self.oc_info_val(info_desc.type, val)
            elif info_desc.num in [-2, "G"]:  # Number=G
                oc_val = None  # TODO handle Number=G
            elif info_desc.num in [-3, "R"]:  # Number=R
                val = info_val[gt]  # TODO find an example and make sure this is right
                oc_val = self.oc_info_val(info_desc.type, val)
            elif info_desc.num in [None, "."]:  # Number=.
                tmp = lambda val: self.oc_info_val(info_desc.type, val, force_str=True)
                oc_val = ",".join(map(tmp, info_val))
            elif info_desc.num == 1:
                oc_val = self.oc_info_val(info_desc.type, info_val)
            else:  # Number>1
                tmp = lambda val: self.oc_info_val(info_desc.type, val, force_str=True)
                oc_val = ",".join(map(tmp, info_val))
            row_data[info_name] = oc_val
        alt = self.curvar.ALT[gt - 1].sequence  # pyright: ignore
        row_data["pos"] = self.curvar.POS
        row_data["ref"] = self.curvar.REF
        row_data["alt"] = alt
        if "genotype" in wdict:
            row_data["genotype"] = wdict["genotype"]
        if self.cur_csq:
            # TODO csq alts and multiallelic sites
            if len(self.curvar.ALT) == 1:
                if self.cur_csq:
                    row_data.update(next(iter(self.cur_csq.values())))
            else:
                row_data.update(self.cur_csq.get(alt, {}))
        self.ex_info_writer.write_data(row_data)
