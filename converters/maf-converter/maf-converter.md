# maf-converter
Converts [MAF](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/#:~:text=Mutation%20Annotation%20Format%20(MAF)%20is,through%20the%20Somatic%20Aggregation%20Workflow.) format files



A standard MAF file consists of 126 columns which contain tab-delimited values with aggregated mutation information from VCF Files.
A 10 column example of a MAF file
```
Hugo_Symbol	Entrez_Gene_Id	Center	NCBI_Build	Chromosome	Start_Position	End_Position	Strand	Variant_Classification	Variant_Type
FGFR2			GRCh38	chr10	121593817	121593817	+	Start_Lost	SNP
```

For a file to be regarded as a MAF file it must pass 12 steps of verification.
