# PMKB: Open source archive of clinical interpretations of cancer variants in a structured way

PMKB(Precision Medicine Knowledge Base) is an open source database for interpretation of cancer variants. PMKB will take a name of the variant as an input and match it to all variants available and return an interpretation for that specific variant


## Data Granularity and Variant Description

- Variant categories include SNVs, Indels, CNV, and Gene fusions.
- Description of small, localized mutations such as SNVs, and indels can be broken down into to groups:
    1. Specific Mutation described using HGVS protein change and DNA change notation
    2. Gene-region-based description:
        - When a variant is entered; PMKB retreives specific gene region information from ENSEMBL based on its canonical transcript for a gene and its GRCh37-based API
        - PMKB is set up to take variant's HGVS protein notation as input and match that variant (eg. KRAS p.G12A) against multiple level of variant descriptions (eg. KRAS mutations), and return all relevant interpretations
- PMKB is under the lincense of Creative Common 4.0 Internation

## You can always refer to PMKB published paper:
    Linda Huang and others, The cancer precision medicine knowledge base for structured clinical-grade mutations and interpretations, Journal of the American Medical Informatics Association, Volume 24, Issue 3, May 2017, Pages 513–519, https://doi.org/10.1093/jamia/ocw148


