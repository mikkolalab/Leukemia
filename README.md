# Requirements

- Python 3.10.0
- STAR aligner v2.7.10a
- R 4.2.2
- Rstudio 2024.04.1+748

## Variant calling

This script runs the `run_pysam_varcall` to call variants from a list of bam files.
It uses the pysam library to call variants.
Usage:
```
python run_pysam_varcall.py \
       --bam_fofn bam.fofn.txt \
       --output_file_name output.tsv \
```
Where the fofn file is a tab-separated file with the bam files you wish to use for the variant calling:

```    
bam   sample   cluster
sample1.cls1.bam   sample1       cluster1
sample1.cls2.bam   sample1       cluster2
```

If the bam files are not spliut by cluster you can put NA in the cluster column.


## 2-pass STAR alignment 

This script runs the STAR aligner in two passes to improve the aligmnet of split reads. 
In the first pass, it skips novel junctions longer than 100 kb. It also applies more stringent 
parameters for the alignment such as the overhang length and the number of mismatches.

Usage:
```
qsub -t 1-N remap.two_pass_start.sh
```

This script is designed to run in the Hoffman2 cluster (or other SGE clusters). 
It reads the input from a file called fastq.job:

```
sample1.R1.fastq.gz sample1.R2.fastq.gz sample1
sample2.R1.fastq.gz sample2.R2.fastq.gz sample2
```

because there are two lines, the script will run two jobs (N=2).
The STAR index, GTF file and output directory are hardcoded in the script. Please change to run
new samples.


## Subset bam file by cell clusters 

This script splits the reads from the bulk bam file by cell clusters. 
It uses the Seurat object to link each barcode to its cluster. 
It uses the Fastq file to link the barcode to read id.
Finally, for each read, it looks up the cluster assigment and writes a new bam file for that cluster.


Usage:
```
python subset_metacell_mecom_sampels.py \
	-s sample1 \
	-b sample1.bam \
	-f sample1.R1.fastq.gz \
	-m seurat.metadata.csv 
```

You can convert the seurat object into a csv file in R
```
df <- readRDS('seurat_obj.rds')
write.csv(df@meta.data, 'seurat.metadata.csv')
```

Seurat csv:
```
"","orig.ident","nCount_RNA","nFeature_RNA","percent.mito","RNA_snn_res.0.6","Cell_type"
"AAACCTGAGAATCTCC-1","BM-1",10189,3388,0.0104033761900088,"6","6"
"AAACCTGAGAGCCCAA-1","BM-1",10200,3079,0.0282352941176471,"3","3"
"AAACCTGAGATGCCAG-1","BM-1",6996,2389,0.0118639222412807,"0","0"
"AAACCTGAGCACACAG-1","BM-1",4071,1456,0.0203881110292311,"0","0"
"AAACCTGAGCGATTCT-1","BM-1",4230,1561,0.0191489361702128,"0","0"
```

The first column should contain the barcodes (it can include dash).

The column `orig.ident` contains all the samples in the object. Select one of these values to input in the 
python script (-s option).  

The script will split the bam file based on the `Cell_type` column. 

You may need to edit the Seurat obhect to match this format.


## Make density plots

This script makes a density plot of the reads in a gene of interest.

First, you will need to generate a sashimi file using the `ggsashimi.py` script. 
Then, you can use the `plot_density.R` script to generate the plot.

Usage:
```
python ggsashimi.py \
       -b FOFN.txt \
       -o sashimi.gene1 \
       --gtf gencode.v44.comprehensive.annotation.gtf.gz \
       --gene-gtf gencode.v44.basic.annotation.gene_only.gtf \
       --gene gene1 
```
`--gene` is the gene of interest.
`--gtf` and `--gene-gtf` are the gtf files. You can use the same file for both options, or subset a GTF 
to only include records for the gene of interest for the --gene-gtf option.

The FOFN file (`--b` option) is a tab-separated file with the bam files you wish to use for the density plot:

```
id1    sample1.bam   sample1
```
where the columns are id, path to bam file and sample name. The id is only used to tell samples apart. 
The sample name is used to label the samples in the plot.

Second, the `plot_density.R` script will generate the plot. It takes as input the sashimi file and the `--gene-gtf` file.
See the .Rmd file for instructions on how to run the script.


