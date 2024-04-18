# gvf-converter
Converts [GVF](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gvf.md) and [GFF/GFF3](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) format files


A standard GVF/GFF/GFF3 file consists of 9 columns which contain tab-delimited or space-delimited values. With Refrence_seq and Variant_seq in attributes. Below is an example of the format:

A 9 column example of a GVF/GFF/GFF3 file
```
Chromosome	source	type	start	end	score	strand	phase	attributes
chr16	samtools	SNV	49291141	49291141	.	+	.	ID=ID_1;Variant_seq=A,G;Reference_seq=G;
```