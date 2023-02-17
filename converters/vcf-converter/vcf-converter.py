from oakvar import BaseConverter
from oakvar import IgnoredVariant
import re
from collections import defaultdict
from oakvar.util.inout import CravatWriter
from io import StringIO
import copy
from pathlib import Path
from math import isnan
from collections import OrderedDict
from typing import Any

class Converter(BaseConverter):

    def __init__(self):
        self.format_name = 'vcf'
        self._in_header = True
        self._first_variant = True
        self._buffer = StringIO()
        self._reader = None
        self.addl_cols = [
            {'name':'phred','title':'Phred','type':'string'},
            {'name':'filter','title':'VCF filter','type':'string'},
            {'name':'zygosity','title':'Zygosity','type':'string'},
            {'name':'alt_reads','title':'Alternate reads','type':'int'},
            {'name':'tot_reads','title':'Total reads','type':'int'},
            {'name':'af','title':'Variant allele frequency','type':'float'},
            {'name':'hap_block','title':'Haplotype block ID','type':'int'},
            {'name':'hap_strand','title':'Haplotype strand ID','type':'int'},
        ]
        self.ex_info_writer = None
        self.curvar = None
        self.csq_fields = None
        self.input_assembly = None

    def check_format(self, f): 
        if f.name.endswith('.vcf'):
            return True
        if f.name.endswith('.vcf.gz'):
            return True
        first_line = f.readline()
        if first_line.startswith('##fileformat=VCF'):
            return True

    def detect_genome_assembly_from_str(self, s):
        if "GRCh37" in s or "hg19" in s:
            return "hg19"
        if "NCBI36" in s or "hg18" in s:
            return "hg18"
        return None

    def detect_genome_assembly(self, reader):
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
                        for kkk, vvv in vv.items():
                            assembly = self.detect_genome_assembly_from_str(vvv)
                            if assembly:
                                return assembly
            elif type(v) == str:
                assembly = self.detect_genome_assembly_from_str(v)
                if assembly:
                    return assembly
        return None

    def setup(self, f):
        if hasattr(self, 'conf') == False:
            self.conf = {}
        if type(self.conf.get('exclude_info')) == str:
            self.exclude_info = set(self.conf['exclude_info'].split(','))
        else:
            self.exclude_info = set()
        if type(self.conf.get('include_info')) == str:
            self.include_info = set(self.conf['include_info'].split(','))
        else:
            self.include_info = set()
        import vcf
        reader = vcf.Reader(f, compressed=False)
        self.fix_formats(reader)
        self.open_extra_info(reader)
        self.input_assembly = self.detect_genome_assembly(reader)

    def fix_formats(self, reader):
        return
        # A user had AD number=1
        if 'AD' in reader.formats:
            reader.formats['AD']._replace(num = -3)

    def open_extra_info(self, reader):
        #if not reader.infos:
        #    return
        if not self.output_dir or not self.run_name:
            return
        writer_path = Path(self.output_dir)/(self.run_name+'.extra_vcf_info.var')
        self.ex_info_writer = CravatWriter(str(writer_path))
        info_cols: list[dict[ str, Any ] ] = [{'name':'uid','title':'UID','type':'int'}]
        info_cols.append({
            'name': 'pos',
            'title': 'VCF Position',
            'desc': '',
            'type': 'int',
            'width': 60,
        })
        info_cols.append({
            'name': 'ref',
            'title': 'VCF Ref Allele',
            'desc': '',
            'type': 'string',
            'width': 60,
        })
        info_cols.append({
            'name': 'alt',
            'title': 'VCF Alt Allele',
            'desc': '',
            'type': 'string',
            'width': 60,
        })
        typemap = {'Integer':'int','Float':'float'}
        if reader.infos:
            for info in reader.infos.values():
                info_cols.append({
                    'name': info.id,
                    'title': info.id,
                    'desc': info.desc,
                    'type': typemap.get(info.type,'string'),
                    'hidden': True,
                })
            if 'CSQ' in reader.infos:
                csq_info = reader.infos['CSQ']
                fields_match = re.search(r'Format: ([^\s]+)', csq_info.desc)
                if fields_match:
                    self.csq_fields = ['CSQ_'+x for x in fields_match.group(1).split('|')]
                    for cname in self.csq_fields:
                        info_cols.append({
                            'name': cname,
                            'title': cname.replace('_',' '),
                            'type': 'string',
                            'hidden': True,
                        })
        if self.include_info:
            self.include_info.update([c['name'] for c in info_cols[:3]])
            temp = info_cols
            info_cols = []
            for col in temp:
                if col['name'] in self.include_info and col['name'] not in self.exclude_info:
                    col['hidden'] = False
                    info_cols.append(col)
            del temp
        else:
            temp = info_cols
            info_cols = []
            for col in temp:
                if col['name'] not in self.exclude_info:
                    col['hidden'] = False
                    info_cols.append(col)
            del temp
        self.info_cols = [v['name'] for v in info_cols]
        self.ex_info_writer.add_columns(info_cols)
        self.ex_info_writer.write_definition()
        self.ex_info_writer.write_meta_line('name', 'extra_vcf_info')
        self.ex_info_writer.write_meta_line('displayname', 'Extra VCF INFO Annotations')

    def convert_line(self, l):
        import vcf
        if l.startswith('#'):
            if self._in_header:
                self._buffer.write(l)
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
        self._buffer.write(l)
        self._buffer.seek(0)
        try:
            variant = next(self._reader)
        except StopIteration:
            return self.IGNORE
        wdict_blanks = {}
        for alt_index, alt in enumerate(variant.ALT):
            if alt is None:
                alt_base = variant.REF
            elif alt.type == 'NON_REF':
                alt_base = None
            else:
                alt_base = alt.sequence
            new_pos, new_ref, new_alt = self.trim_variant(variant.POS, variant.REF, alt_base)
            #new_pos, new_ref, new_alt = variant.POS, variant.REF, alt_base
            if variant.FILTER is None:
                filter_val = None
            elif len(variant.FILTER) == 0:
                filter_val = 'PASS'
            else:
                filter_val = ';'.join(variant.FILTER)
            wdict_blanks[alt_index+1] = {
                'chrom': variant.CHROM,
                'pos': new_pos,
                'ref_base': new_ref,
                'alt_base': new_alt,
                'tags': variant.ID,
                'phred': variant.QUAL,
                'filter': filter_val,
            }
        wdicts = []
        self.gt_occur = []
        if len(variant.samples) > 0:
            all_gt_zero = True
            for call in variant.samples:
                # Dedup gt but maintain order
                for gt in list(OrderedDict.fromkeys(call.gt_alleles)):
                    if gt == '0' or gt is None:
                        continue
                    all_gt_zero = False
                    gt = int(gt)
                    wdict = copy.copy(wdict_blanks[gt])
                    if wdict['alt_base'] == '*':
                        continue
                    wdict['sample_id'] = call.sample
                    if call.is_het == True:
                        wdict['zygosity'] = 'het'
                    elif call.is_het == False:
                        wdict['zygosity'] = 'hom'
                    elif call.is_het is None:
                        wdict['zygosity'] = None
                    else:
                        wdict['zygosity'] = None
                    wdict['tot_reads'], wdict['alt_reads'], wdict['af'] = self.extract_read_info(call, gt)
                    wdict['hap_block'] = None #TODO
                    wdict['hap_strand'] = None #TODO
                    wdicts.append(wdict)
                    self.gt_occur.append(gt)
            if all_gt_zero:
                e = IgnoredVariant('All samples have the reference genotype.')
                e.traceback = False
                raise e
        else:
            for gt in wdict_blanks:
                wdict = copy.copy(wdict_blanks[gt])
                if wdict['alt_base'] == '*':
                        continue
                wdicts.append(wdict)
                self.gt_occur.append(gt)
        self.curvar = variant
        self.cur_csq = {}
        if self.csq_fields and 'CSQ' in variant.INFO:
            csq_entries = defaultdict(list)
            for gt_csq in variant.INFO['CSQ']:
                l = gt_csq.split('|')
                csq_entries[l[0]].append(l)
            for allele, entries in csq_entries.items():
                transpose = zip(*entries)
                self.cur_csq[allele] = dict([(cname, self.csq_format(value)) for cname, value in zip(self.csq_fields, transpose)])
        return wdicts

    @staticmethod
    def extract_read_info(call, gt):
        tot_reads = None
        alt_reads = None
        # AD is depth for each allele
        if hasattr(call.data,'AD'):
            # tot_reads
            if hasattr(call.data.AD, '__iter__'):
                tot_reads = sum([0 if x is None else int(x) for x in call.data.AD])
            elif call.data.AD is None:
                tot_reads = 0
            else:
                tot_reads = int(call.data.AD)
            # alt_reads
            try:
                alt_reads = int(call.data.AD[gt])
            except IndexError: # Wrong length
                alt_reads = None
            except TypeError: # Not indexable
                alt_reads = int(call.data.AD)
        # DP is total depth
        if hasattr(call.data,'DP'):
            tot_reads = int(call.data.DP)
        if tot_reads is not None and alt_reads is not None:
            try:
                alt_freq = alt_reads/tot_reads
            except ZeroDivisionError:
                alt_freq = None
        else:
            alt_freq = None
        return tot_reads, alt_reads, alt_freq

    
    @staticmethod
    def csq_format(l):
        # Format a list of CSQ values into it's representation in OC
        # Each value comes from a VEP transcript mapping
        if all([x=='' for x in l]):
            return None
        else:
            return ';'.join(l)

    def trim_variant(self, pos, ref, alt):
        if alt is None:
            return pos, ref, alt
        if len(ref) == 1 and len(alt) == 1:
            return pos, ref, alt
        ref = list(ref)
        alt = list(alt)
        adj = 0
        while ref and alt and ref[0]==alt[0]:
            adj += 1
            ref.pop(0)
            alt.pop(0)
        while ref and alt and ref[-1]==alt[-1]:
            ref.pop()
            alt.pop()
        ref = ''.join(ref) if ref else '-'
        alt = ''.join(alt) if alt else '-'
        return pos+adj, ref, alt

    @staticmethod
    def oc_info_val(info_type, val, force_str=False):
        if val is None or val=='.':
            oc_val = None
        if info_type in ('Integer','Float'):
            if isnan(val):
                oc_val = None
            else:
                oc_val = val
        else:
            oc_val = val
        if force_str and oc_val is None:
            return '.'
        elif force_str:
            return str(oc_val)
        else:
            return oc_val

    def addl_operation_for_unique_variant (self, wdict, wdict_no):
        if self.ex_info_writer is None:
            return
        gt = self.gt_occur[wdict_no]
        alt_index = gt-1
        row_data = {'uid':wdict['uid']}
        if not self.curvar:
            return
        for info_name, info_val in self.curvar.INFO.items():
            if info_name not in self.info_cols:
                continue
            if not self._reader or not self._reader.infos:
                continue
            info_desc = self._reader.infos[info_name]
            if info_desc.num == 0:
                oc_val = self.oc_info_val(info_desc.type, info_val)
            elif info_desc.num == -1: # Number=A
                oc_val = self.oc_info_val(info_desc.type, info_val[alt_index])
            elif info_desc.num == -2: # Number=G
                oc_val = None #TODO handle Number=G
            elif info_desc.num == -3: # Number=R
                val = info_val[gt] #TODO find an example and make sure this is right
                oc_val = self.oc_info_val(info_desc.type, val)
            elif info_desc.num is None: # Number=.
                tmp = lambda val: self.oc_info_val(info_desc.type, val, force_str=True)
                oc_val = ','.join(map(tmp, info_val))
            elif info_desc.num == 1:
                oc_val = self.oc_info_val(info_desc.type, info_val)
            else: # Number>1
                tmp = lambda val: self.oc_info_val(info_desc.type, val, force_str=True)
                oc_val = ','.join(map(tmp, info_val))
            row_data[info_name] = oc_val
        alt = self.curvar.ALT[gt-1].sequence # pyright: ignore
        row_data['pos'] = self.curvar.POS
        row_data['ref'] = self.curvar.REF
        row_data['alt'] = alt
        if self.cur_csq:
            #TODO csq alts and multiallelic sites
            if len(self.curvar.ALT) == 1:
                if self.cur_csq:
                    row_data.update(next(iter(self.cur_csq.values())))
            else:
                row_data.update(self.cur_csq.get(alt,{}))
        self.ex_info_writer.write_data(row_data)
