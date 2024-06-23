from oakvar import BaseAnnotator

class Annotator(BaseAnnotator):
    def setup(self):
        pass
    
    def split_pval_bits(self, pval_bits):
        pval_coeff = (pval_bits >> 8) & 0xFF
        pval_exp = -(pval_bits & 0xFF)
        return pval_coeff, pval_exp
        
    def format_p_value(self, p_value):
        return format(p_value, '.9g')

    def annotate(self, input_data, secondary_data):
        rsid = None
        
        if 'ID' in input_data and input_data['ID'].startswith('rs'):
            rsid = int(input_data['ID'].strip('rs'))
        elif 'dbsnp' in secondary_data and secondary_data['dbsnp'] is not None and len(secondary_data['dbsnp']) > 0:
            rsid = int(secondary_data['dbsnp'][0]['rsid'].lstrip('rs'))
        
        if rsid:
            all_results = []
            bin_index = rsid % 1000
            table_name = 'bin_{}'.format(bin_index)
            q = f"SELECT data FROM {table_name} WHERE rsid = ?"
            self.cursor.execute(q, [rsid])
            rows = self.cursor.fetchall()
            for row in rows:
                data = row[0]
                if data:
                    for annotation in data:
                        allele = annotation['allele']
                        trait = annotation['trait']
                        pval_bits = annotation['pval_bits']
                        pval_info = annotation['pval_info']
                        or_beta = annotation['or_beta']
                        pubmed_id = annotation['pubmed_id']
                        pval_coeff, pval_exp = self.split_pval_bits(pval_bits)
                        p_value = pval_coeff * (10 ** pval_exp)
                        if or_beta is not None:
                            or_beta = round(or_beta, 6)
                        all_results.append({
                            'allele': allele,
                            'trait': trait,
                            'p_value': self.format_p_value(p_value),
                            'or_beta': round(or_beta, 9) if or_beta is not None else None,
                            'pubmed_id': pubmed_id
                        })
            result = {'all': all_results} if all_results else None
            return result
        
        return None
    
    def cleanup(self):
        pass


if __name__ == '__main__':
    annotator = Annotator()
    annotator.run()
