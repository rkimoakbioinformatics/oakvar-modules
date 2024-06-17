from oakvar import BaseAnnotator

class Annotator(BaseAnnotator):

    def annotate(self, input_data):
        if self.cursor is None:
            return None
        af_col_names = ['af','af_afr','af_amr','af_asj','af_eas',
                        'af_fin','af_nfe','af_oth', 'af_sas']
        out = {x:'' for x in af_col_names}

        chrom = input_data['chrom']
        pos = input_data['pos']
        ref = input_data['ref_base']
        alt = input_data['alt_base']
        
        q = 'select {} from {} where pos=? and ref=? and alt=?'.format(
            ', '.join(af_col_names),
            chrom
        )
        self.cursor.execute(q, (pos, ref, alt))
        qr = self.cursor.fetchone()
        if qr:
            for i, k in enumerate(af_col_names):
                out[k] = qr[i]
            # Patch in latino for american
            out['af_lat'] = out['af_amr']
            del out['af_amr']
            return out
