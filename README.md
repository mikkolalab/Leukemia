## 2-pass STAR alignment 

remap.two_pass_start.sh


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


python ggsashimi.py \
       -b data/plots/densities.${dataset}_samples.fofn.txt \
       -o data/plots/densities.${dataset}_${GENE} \
       --gtf gencode.v44.comprehensive.annotation.gtf.gz \
       --gene-gtf gencode.v44.basic.annotation.gene_only.gtf \
       --gene ${GENE}
