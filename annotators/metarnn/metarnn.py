from oakvar import BaseAnnotator


class Annotator(BaseAnnotator):
    valid_chroms = [
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chr20",
        "chr21",
        "chr22",
        "chrX",
        "chrY",
        "chrM",
    ]

    def annotate(self, input_data):
        if not self.cursor:
            return None
        if input_data["chrom"] not in self.valid_chroms:
            return None
        chrom = input_data["chrom"]
        self.cursor.execute(
            f"select score, rankscore, pred from {chrom} where pos=? and ref=? and alt=?;",
            (input_data["pos"], input_data["ref_base"], input_data["alt_base"]),
        )
        qr = self.cursor.fetchone()
        if qr is not None:
            return {
                "score": qr[0],
                "rankscore": qr[1],
                "pred": qr[2]
            }
