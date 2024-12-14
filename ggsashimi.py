#!/usr/bin/env python
import subprocess as sp
import sys, re, copy, os, codecs, gzip
from argparse import ArgumentParser, Action as ArgParseAction
from collections import OrderedDict, Counter, defaultdict
import pysam
import tabix
import numpy as np
import pandas as pd
import HTSeq
import tabix


def read_is_bad(read):
        if read.is_unmapped:
                return True
        if read.is_duplicate:
                return True
        if read.is_secondary:
                return True
        if read.is_supplementary:
                return True
        if read.mapping_quality < 20:
                return True


def find_high_cov_regions(bam_file, c, max_cov_cutoff = 100000):
        cov =  Counter()
        r_count = 0

        if not os.path.exists(bam_file + ".bai"):
                sys.stderr.write(f"Indexing {bam_file}\n")
                bam_handle = pysam.AlignmentFile(bam_file)
                pysam.index(bam_file)
                bam_handle.close()
        
        bam_handle = pysam.AlignmentFile(bam_file)
        
        for read in bam_handle.fetch(c.chrom, c.start, c.end):
                r_count += 1 
                if read_is_bad(read):
                        continue
                for bs, be in read.get_blocks():
                        cov[bs] += 1
                        cov[be] -= 1
        bam_handle.close()

        cum_dp = 0 
        cov_cum = defaultdict(int)

        for pos, dp in sorted(cov.items()):
                cum_dp += dp
                cov_cum[pos] = cum_dp 

        return cov_cum


def read_bam_input(f):
        ANNOT = {}
        with codecs.open(f, encoding='utf-8') as openf:
                for line in openf:
                        line_sp = line.strip().split("\t")
                        bam_index = line_sp[0]
                        bam_file  = line_sp[1]
                        bam_label = line_sp[2] 
                        ANNOT[bam_index] = [bam_file, bam_label]
        return ANNOT


def get_gene_coordinates(gtf_file, gene_name):
        iv = None
        f = HTSeq.GFF_Reader(gtf_file)
        for feature in f:
                if feature.type == "gene" and feature.attr["gene_name"] == gene_name:
                        return feature.iv
        return iv 

def write_subregion_gtf(gtf_file, iv, new_gtf_file):
        gtf = tabix.open(gtf_file)
        for record in gtf.query(iv.chrom, iv.start, iv.end):
                record_line = "\t".join(record)
                print(record_line, file = new_gtf_file)


def main():
        print(' '.join(sys.argv))
        parser = ArgumentParser(description='Create sashimi plot for a given genomic region')
        
        parser.add_argument("-b", "--bam", 
                type = str, 
                required = True,
                help = """
                Individual bam file or file with a list of bam files.
                In the case of a list of files the format is tsv:
                1col: IDX for bam file,
                2col: path of bam file,
                3+col: additional columns
                """)
        
        parser.add_argument("-o", "--out-prefix", 
                type = str, 
                dest = "out_prefix", 
                required = True,
                help = "Prefix for plot file name [default = %(default)s]")

        parser.add_argument("--gene", 
                type = str,
                required = True,
                help = "Gene Name")

        parser.add_argument("--gtf",
                type = str,
                required = True,
                help = "GTF file")

        parser.add_argument("--gene-gtf",
                type = str,
                dest = "gene_gtf",
                required = True,
                help = "GTF file")
       
        args = parser.parse_args()

        METADATA = read_bam_input(args.bam)

        

        coordinates = get_gene_coordinates(args.gene_gtf, args.gene)

        target_gene_gtf = open(f"{args.out_prefix}.{args.gene}.gtf", "w")
        write_subregion_gtf(args.gtf, coordinates, target_gene_gtf)

        outf = open(args.out_prefix + ".sashimi", "w")

        bam_dict = defaultdict(dict)
        for IDX, (bam_file, label_text) in METADATA.items():    
                sys.stderr.write(f"Processing {IDX} {label_text}\n")            
                cum_cov = find_high_cov_regions(bam_file, coordinates, max_cov_cutoff = 100000)
                for pos, dp in sorted(cum_cov.items()):
                        outf.write(f"{IDX}\t{label_text}\t{pos - 0.5}\t0\n")
                        outf.write(f"{IDX}\t{label_text}\t{pos}\t{dp}\n")
                        outf.write(f"{IDX}\t{label_text}\t{pos + 0.5}\t0\n")

        outf.close()
        target_gene_gtf.close()


if __name__ == '__main__':
        main()
