# MAF Reporter

Creates a MAF report by mapping values contained in OakVar database from multiple previously-run OakVar modules.

The MAF Reporter is created according to the 
[Mutation Annotation Format (MAF)](https://docs.gdc.cancer.gov/Encyclopedia/pages/Mutation_Annotation_Format_TCGAv2/) 
documentation and details from [GDC MAF Format](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/#protected-maf-file-structure) specifications.

The MAF Reporter works by mapping values from previously run OakVar modules, such as clinvar, dbSNP, ensembl etc.
If these modules were not previously run some columns might result in being empty.