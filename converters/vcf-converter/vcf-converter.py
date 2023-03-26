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
from cyvcf2 import VCF


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

    def check_format(self, f):
        if f.name.endswith(".vcf"):
            return True
        if f.name.endswith(".vcf.gz"):
            return True
        first_line = f.readline()
        if first_line.startswith("##fileformat=VCF"):
            return True

    def detect_genome_assembly_from_str(self, s):
        if "GRCh37" in s or "hg19" in s:
            return "hg19"
        if "NCBI36" in s or "hg18" in s:
            return "hg18"
        return None

    def detect_genome_assembly(self, reader, f):
        genome_assembly = self.detect_genome_assembly_from_metadata(reader)
        if genome_assembly:
            return genome_assembly
        return self.detect_genome_assembly_from_contigs(f)

    def get_value_from_contig_line(self, value, line):
        if f"{value}=" not in line:
            return None
        return line.split("assembly=")[1].strip().strip(">").split(",")[0].lower()

    def get_assembly_from_contig_line(self, line):
        return self.get_value_from_contig_line("assembly", line)

    def get_length_from_contig_line(self, line):
        return self.get_value_from_contig_line("length", line)

    def detect_genome_assembly_from_contigs(self, f):
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

    def setup(self, f):
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
        import vcf

        reader = vcf.Reader(f)
        self.open_extra_info(reader)
        self.input_assembly = self.detect_genome_assembly(reader, f)

    def get_do_liftover_chrM(self, genome_assembly, f, do_liftover):
        return self.chrM_needs_liftover(genome_assembly, f, do_liftover)

    def chrM_needs_liftover(self, genome_assembly, f, do_liftover):
        import vcf

        if genome_assembly != "hg19":
            return do_liftover
        fpath = f.name
        reader = vcf.Reader(fpath)
        if not reader.contigs:
            return do_liftover
        for contig in reader.contigs.values():
            if contig.id in ["M", "MT", "Mt"]:
                if contig.length == 16571:
                    f.seek(0)
                    return True
                elif contig.length == 16569:
                    f.seek(0)
                    return False
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

    def convert_file(
        self,
        file,
        *__args__,
        exc_handler=None,
        ignore_sample: bool = False,
        **__kwargs__,
    ):
        """Reads a VCF file as a list of dictionaries using cyvcf2.

        Args:
            file (TextIOWrapper): A VCF file read as a text stream.
        """
        self._reader = VCF(file)
        line_no = (
            self._reader.raw_header.count("\n") - 1
        )  # Count number of lines in the header
        self.cyvcf_samples = self._reader.samples
        self.len_cyvcf_samples = len(self.cyvcf_samples)
        for variant in self._reader:
            line_no += 1
            try:
                yield line_no, self.convert_line(variant, ignore_sample=ignore_sample)
            except Exception as e:
                if exc_handler:
                    exc_handler(line_no, e)
                else:
                    raise e
        return None

    def convert_line(self, variant, ignore_sample: bool = False):
        """Converts a variant from a VCF into a list of dictionaries.

        Args:
            variant (cyvcf2.Variant): A single variant parsed from the VCF file using cyvcf2.
        """
        try:
            from oakvar.lib.exceptions import NoAlternateAllele  # type: ignore
        except:
            from oakvar.exceptions import NoAlternateAllele  # type: ignore

        wdict_blanks = {}
        wdicts = []
        # Empty alternate alleles do not get read by cyvcf2
        var_alt = variant.ALT
        if len(var_alt) == 0:
            var_alt.append(None)
        for alt_index, alt in enumerate(variant.ALT):
            if alt is None:
                alt_base = variant.REF
            else:
                alt_base = alt
            new_pos, new_ref, new_alt = self.trim_variant(
                variant.POS, variant.REF, alt_base
            )
            if variant.FILTER is None:
                filter_val = "PASS"
            else:
                filter_val = variant.FILTER
            wdict_blanks[alt_index + 1] = {
                "chrom": variant.CHROM,
                "pos": new_pos,
                "ref_base": new_ref,
                "alt_base": new_alt,
                "tags": variant.ID,
                "phred": variant.QUAL,
                "filter": filter_val,
            }
        self.gt_occur: List[int] = []
        try:
            # Preload data since cyvcf2 is slow otherwise
            var_genotypes = variant.genotypes
            var_gt_types = variant.gt_types
            var_gt_bases = variant.gt_bases
            no_genotype = False
        except Exception:
            no_genotype = True  # Somatic variants with no genotypes
        no_sample = len(self.cyvcf_samples) == 0
        if no_sample or no_genotype:
            for gt in wdict_blanks:
                wdicts.append(wdict_blanks[gt])
                self.gt_occur.append(gt)
        else:
            all_gt_zero = True
            for i in range(self.len_cyvcf_samples):
                gt_sample = var_genotypes[i]
                var_ploidy = (
                    variant.ploidy
                )  # Length of genotype array is (ploidy + 1), where the last value is phase
                # Dedup gt but maintain order
                for gt in list(OrderedDict.fromkeys(gt_sample[:var_ploidy])):
                    if gt in [None, 0, -1]:
                        continue
                    all_gt_zero = False
                    self.gt_occur.append(gt)
                    if ignore_sample and gt not in self.gt_occur:
                        wdicts.append(wdict_blanks[gt])
                    else:
                        wdict = copy.copy(wdict_blanks[gt])
                        wdict["sample_id"] = self.cyvcf_samples[i]
                        if var_gt_types[i] == 1:
                            wdict["zygosity"] = "het"
                        elif var_gt_types[i] in [0, 3]:
                            wdict["zygosity"] = "hom"
                        else:
                            wdict["zygosity"] = None
                        (
                            wdict["tot_reads"],
                            wdict["alt_reads"],
                            wdict["af"],
                        ) = self.extract_read_info(variant, i, gt)
                        wdict["hap_block"] = None  # TODO
                        wdict["hap_strand"] = None  # TODO
                        wdict["genotype"] = var_gt_bases[i]
                        wdicts.append(wdict)
                if ignore_sample and not all_gt_zero:
                    break
            if all_gt_zero:
                raise NoAlternateAllele()
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
    def extract_read_info(variant, sample_num: int, gt):
        """Extracts read information from a variant for a particular sample.

        Args:
            variant (cyvcf2.Variant): A single variant parsed from the VCF file using cyvcf2.
            sample_num (int): The index of the relevant sample in the list of samples.
            gt (int): The genotype of the sample.
        """
        tot_reads: Optional[int] = None
        alt_reads: Optional[int] = None
        # AD is depth for each allele
        if "AD" in variant.FORMAT:
            var_depths = variant.format("AD")
            # Compute total reads if possible
            if hasattr(var_depths, "__iter__"):
                tot_reads = sum([0 if x is None else x for x in var_depths[sample_num]])
            elif var_depths is None:
                tot_reads = 0
            else:
                tot_reads = var_depths
            # Compute alternate allele reads if possible
            if -1 in var_depths[sample_num]:
                alt_reads = None
            else:
                alt_reads = var_depths[sample_num][gt]
        else:
            alt_reads = variant.gt_alt_depths[sample_num]
        # DP is total depth
        if "DP" in variant.FORMAT:
            variant_dp = variant.format("DP")
            tot_reads = variant_dp[sample_num][0]
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
        gt: int = self.gt_occur[wdict_no]
        alt_index = gt - 1
        row_data = {"uid": wdict["uid"]}
        if not self.curvar:
            return
        somatic = True if self.curvar.INFO.get("SOMATIC") else False
        for info_name, info_val in self.curvar.INFO:
            if info_name not in self.info_cols:
                continue
            if not self._reader:
                continue
            info_desc = self._reader.get_header_type(info_name)
            if info_desc["Number"] == "0":
                oc_val = self.oc_info_val(info_desc["Type"], info_val)
            elif info_desc["Number"] in ["-1", "A"]:  # Number=A
                oc_val = self.oc_info_val(info_desc["Type"], info_val[alt_index])
            elif info_desc["Number"] in ["-2", "G"]:  # Number=G
                oc_val = None  # TODO handle Number=G
            elif info_desc["Number"] in ["-3", "R"]:  # Number=R
                val = info_val[gt] if not somatic else info_val[gt - 1]  # First index is REF if not somatic VCF
                oc_val = self.oc_info_val(info_desc["Type"], val)
            elif info_desc["Number"] == ".":  # Number=.
                tmp = lambda val: self.oc_info_val(info_desc.type, val, force_str=True)
                oc_val = ",".join(map(tmp, info_val))
            elif info_desc["Number"] == "1":
                oc_val = self.oc_info_val(info_desc.type, info_val)
            else:  # Number>1
                tmp = lambda val: self.oc_info_val(info_desc.type, val, force_str=True)
                oc_val = ",".join(map(tmp, info_val))
            row_data[info_name] = oc_val
        alt = self.curvar.ALT[gt - 1]  # pyright: ignore
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
