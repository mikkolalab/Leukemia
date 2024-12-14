import Bio
import Bio.SeqIO as SeqIO
import Bio.SeqRecord
import Bio.Seq
import sys 
import pandas as pd
import argparse
import gzip
import time 
import pysam
import pickle
import os
from collections import defaultdict
from glob import glob
import time



def print_time_stamp(message):
    sys.stderr.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - {message}\n")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Subset metacell samples')
    parser.add_argument('-s', '--sm', type=str, help='TargetSample')
    parser.add_argument('-f', '--fastq', type=str, help='Fastq file R1')
    parser.add_argument('-b', '--bam', type=str, help='Bam file')
    parser.add_argument('-m', '--metadata', type=str, help='Seurat metadata file')
    args = parser.parse_args()


    print_time_stamp(f"TargetSample: {args.sm}\n")

    ## read Seurat metadata file

    barcodes_split = defaultdict(dict)

    df = pd.read_csv(args.metadata, sep = ',', quotechar = '"')
    df.columns = ['barcodes'] + list(df.columns[1:])

    for row in df.itertuples():
        bc = row.barcodes.split('-')[0] 
        sample = row._2 #orig.ident
        cluster = row.Cell_type # split by cell type 
        barcodes_split[sample][bc] = cluster

    ## get read ids from each barcode

    print_time_stamp(f"barcodes split done.\n")
    
    for seurat_sm, BC_to_Cls in barcodes_split.items():
        
        if args.sm == seurat_sm:
                       
            ReadId_to_BC = defaultdict()

            with gzip.open(args.fastq, "rt") as handle:
                for j, record in enumerate(SeqIO.parse(handle, "fastq"), 1):
      
                    if j % int(1e6) == 0:
                        print_time_stamp(f"\t{j:,} records processed")   

                    bc = str(record.seq)[:16]
                    ReadId_to_BC[str(record.id)] = bc

            print_time_stamp(f"ReadId to BC done. {seurat_sm} {Fastq_File}\n")

            ## Create bam handles for each cluster

            main_bam_handle = pysam.AlignmentFile(args.bam, 'rb')

            Cls_to_Bam = {} 
            for cluster in set(df['CCell_type'].tolist()):
                cls_name = str(cluster).replace(' ', '_').replace('/', '-')
                bam_file_name = Bam_File.replace('.bam', f'_cls_{cls_name}.bam')
                Cls_to_Bam[cluster] = pysam.AlignmentFile(bam_file_name, "wb", template = main_bam_handle)

            print_time_stamp(f"Bam handles created. {seurat_sm} {Bam_File}\n")


            ## Iterate over the main bam file and write to the cluster bam files

            for j, record in enumerate(main_bam_handle.fetch(), 1):
                if j % int(1e6) == 0:
                    print_time_stamp(f"\t{j:,} records processed") 

                BC = ReadId_to_BC.get(record.query_name, None)
                cluster = BC_to_Cls.get(BC, None)

                if cluster is not None:
                    new_record = record
                    new_record.set_tag('CB', BC)
                    new_record.set_tag('CL', cluster)
                    Cls_to_Bam[cluster].write(new_record)

            print_time_stamp(f"Sample done. {seurat_sm}\n")

