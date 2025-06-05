# Variant calling at the single-cell level

1. Umi extraction 
2. Alignment 
3. Variant calling by bc (single-cell level)
4. Integrating into Seurat object


## Umi extraction
*Run time:* 3-12 hours to run

It takes many hours to run 
It takes as input a file called fastq.job that contains fastq paths for R1. 
Keep in mind, if using a different protocol, change --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN C = barcode ; N = UMI

*input: fastq file*

*output: fastq file*

## Alignment 
*Run time:* 2-3 hours to run

You can use the same command as bulk alignment. 

*input: fastq file*

*output: bam file*


## Variant calling by bc 
*Run time:* a few minutes

The script calculates the allelic counts for each variant in the provided BAM files, at the single-cell level. It can handle both short read (sr) and long read (lr) sequencing data.

*input: a lis of bam files and a variants file*

*output: a mutation/variant table*

```
python run_pysam_varcall.umi_bc.copy.py -m example_metadata.tsv -v variants.tsv -o variants.table.out
```


Metadata file should contain the following columns:
```
BAM                 SM  DATA_TYPE
/path/to/file1.bam  sm1 sr
/path/to/file2.bam  sm2 sr
/path/to/file3.bam  sm3 lr
```
where DATA_TYPE can be either `sr` (short read) or `lr` (long read) sequencing data.

NOTE: The bam files must contain the UMI and BC in the read names (using umi_tools extract) or as tags `UB` and `CB` respectively. They should also be sorted by coordinate and indexed.

The variants file should contain the following columns:
```
gene	transcript	    mut         strand	REF	ALT	chrom	pos
SRSF2   NM_003016.4     c.284C>G        -       G       C       chr17   76736877
TET2	NM_001127208.2	c.3909C>G	+       C       G       chr4    105259724
TET2	NM_001127208.2	c.2599T>C	+       T       C       chr4    105236541
TET2	NM_001127208.2	c.5167C>T	+       C       T       chr4    105275677
CBL     NM_005188.3     c.1112A>C       +       A       C       chr11   119278182
```
You can use the website: https://genebe.net/tools/hgvs to convert from HGVS to genomic coordinates.
e.g. NM_003016.4:c.284C>G converts to chr17:76736877 G>C

If missing, you can fill in `transcript` and `gene` with placeholder values.

*Output*
The output will be a file with the following columns:
``` 
        var1_sr_cov var1_sr_ratio   var2_sr_cov var2_sr_ratio   ... varN_sr_cov varN_sr_ratio
sm1 BC1 0.0         NA              1.0         1.0             ... 2.0         0.5   
sm2 BC2 0.0         NA              1.0         1.0             ... 2.0         0.5 
sm3 BC3 0.0         NA              1.0         1.0             ... 2.0         0.5
```

Where `var1`, `var2`, etc. are the variants from the variants file, and `BC1`, `BC2`, etc. are the barcodes from the metadata file.
The entries are the allelic ratio or coverage for each variant

## Integrating into Seurat object
*Run time:*  a few minutes

refer to the jupyter notebook merge_mecom_samples.ipynb 

*input: a mutation/variant table, seurat object*

*output: a seurat object with mutations integrated*

## Denovo mutation calling
*Run time:*  a few minutes

*input: a bam file sorted and indexed*

*output: a vcf file*

script: denovo_variantcall_bcftools.sh