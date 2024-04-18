from oakvar import BaseAnnotator


class Annotator(BaseAnnotator):
    def annotate(self, input_data, secondary_data=None):
        _ = secondary_data
        if not self.cursor:
            return None
        q = 'select samples from cancer where chrom = "{chrom}" and pos = {pos} and ref = "{ref}" and alt = "{alt}"'.format(
            chrom=input_data["chrom"],
            pos=int(input_data["pos"]),
            ref=input_data["ref_base"],
            alt=input_data["alt_base"],
        )
        self.cursor.execute(q)
        row = self.cursor.fetchone()
        if row:
            samples_tmp = str(row[0]).split("|")
            samples = []
            for row in samples_tmp:
                toks = row.split(":")
                samples.append([toks[0], int(toks[1])])
            out = {"samples": samples}
        else:
            out = None
        return out
