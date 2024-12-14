import vcf
import pandas as pd
import numpy as np
import pysam
from collections import defaultdict
from glob import glob
import re
import sys
import os
import multiprocessing as mp
import vcf 
import logging



if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)

    VARS = [("ABL1",  "NM_005157", "+", 763,  "G", "A", "chr9", 130862976),
            ("ABL1",  "NM_005157", "+", 944,  "C", "T", "chr9", 130872896),
            ("ABL1",  "NM_005157", "+", 1075, "T", "G", "chr9", 130873027),

            ('ASXL1', "NM_015338", "+", 2564, "TGATT", "T", "chr20", 32435275),
            ('ASXL1', "NM_015338", "+", 1762, "C", "T", "chr20", 32434474),
            ("CUX1",  "NM_181552", "+", 1,    "A", "G", "chr7", 101817640),
            ("CALR",  "NM_004343", "+", 1154, "A", "AG", "chr19", 12943813),
            ("CBL",   "NM_005188", "+", 1211, "G", "A", "chr11", 119278281),
            
            ("ETV6",  "NM_001987", "+", 613,  "C", "CC", "chr12", 11869573),
            
            ("GATA2", "NM_032638", "-", 58,   "G", "A", "chr3", 128486974),
            
            ("KRAS",  "NM_004985", "-", 38,   "C", "T", "chr12", 25245347),
            ("KRAS",  "NM_004985", "-", 108,  "T", "C", "chr12", 25245277), 
            
            ("NRAS",  "NM_002524", "-", 38,   "C", "T", "chr1", 114716123),
            ("NRAS",  "NM_002524", "-", 35,   "C", "T", "chr1", 114716126),
            
            ("PTPN11", "NM_002834", "+", 226, "G", "A", "chr12", 112450406),
            ("PTPN11", "NM_002834", "+", 1508, "G", "C", "chr12", 112489084),
            ("PTPN11", "NM_002834", "+", 174, "C", "A", "chr12", 112450354),
                
            ("RUNX1", "NM_001754", "-", 973,  "G", "A", "chr21", 34792605),
            ("RUNX1", "NM_001754", "-", 496,  "G", "GG", "chr21", 34880569),
            ("RUNX1", "NM_001754", "-", 592,  "C", "T", "chr21", 34859495),

            ("SF3B1", "NM_012433", "-", 2111, "A", "T", "chr2", 197402097), 
            ("SF3B1", "NM_012433", "-", 2098, "T", "C", "chr2", 197402110),  
            ("SF3B1", "NM_012433", "-", 2219, "C", "T", "chr2", 197401989),
            
            ("SRSF2", "NM_003016", "-", 284,  "G", "C", "chr17", 76736877),
            
            ("TET2",  "NM_001127208", "+", 5167, "C", "T", "chr4", 105275677),
            ("TET2",  "NM_001127208", "+", 2599, "T", "C", "chr4", 105236541),
            ("TET2",  "NM_001127208", "+", 3909, "C", "G", "chr4", 105259724),
            
            ("CBL",   "NM_005188", "+", 1112,  "A", "C", "chr11", 119278182),
            ("DNMT3A",  "NM_022552", "-", 2645, "C", "T", "chr2", 25234373),
            ("FLT3",  "NM_004119", "-", 1780, "N", "NN", "chr13", 28034139),
            ("FLT3",  "NM_004119", "-", 1747, "N", "NN", "chr13", 28034172),
            ("NPM1",  "NM_002520", "+", 860, "N", "NN", "chr5", 171410540 ),
            ("SF3B1", "NM_012433", "-", 1998, "C", "G", "chr2", 197402635), 
            ("PTPN11", "NM_002834", "+", 1504, "T", "C", "chr12", 112489080), 
            ("NRAS",  "NM_002524", "-", 182,   "T", "C", "chr1", 114713908), 
            ("GATA2", "NM_032638", "-", 953,   "G", "C", "chr3", 128483924), 
            ("JAK2",  "NM_004972", "+", 1849,  "G", "T", "chr9", 5073770), ]




    # make a vcf file with these variants
    vcf_file = "data/shortread_10X_5UTR_singlecell_mecom/variants_of_interest.vcf"
    with open(vcf_file, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
        f.write("##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele Count\">\n")
        f.write("##INFO=<ID=AB,Number=1,Type=Float,Description=\"Allele Balance\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for var in VARS:
            gene, transcript, strand, pos, REF, ALT, chrom, POS = var
            f.write(f"{chrom}\t{POS}\t.\t{REF}\t{ALT}\t.\tPASS\tDP={0};AC={0};AB={0}\n")



    rev = {"A" : "T", "T" : "A", "C" : "G", "G" : "C"}

    bam_files = glob('data/shortread_10X_5UTR_singlecell_mecom*/bam/*2passAligned.out.sorted_cls_*.bam')
    bam_files = list(bam_files)

    sys.stderr.write(f"Found {len(bam_files)} \n")
    MAT = defaultdict(dict)

    for bi, bam_file in enumerate(bam_files):
        try:
            pattern = r'data/shortread_10X_5UTR_singlecell_mecom/bam/(.+).star_remap.2passAligned.out.sorted_cls_(.+).bam'
            sample  = re.search(pattern, bam_file).group(1)
            cluster = re.search(pattern, bam_file).group(2)
        except:
            pattern = r'data/shortread_10X_5UTR_singlecell_mecom2/bam/(.+).2passAligned.out.sorted_cls_(.+).bam'
            sample  = re.search(pattern, bam_file).group(1)
            cluster = re.search(pattern, bam_file).group(2)

        logging.info(f"Processing {bi}, {sample} {cluster}")
        
        bh = pysam.AlignmentFile(bam_file, "rb")

        if not bh.has_index():
            pysam.index(bam_file)
            bh.close()
            bh = pysam.AlignmentFile(bam_file, "rb")

        for var in VARS:
            gene, transcript, strand, pos, REF, ALT, chrom, gPOS = var
            
            try: annot = var_annot[(chrom, gPOS)][transcript]
            except: annot = "NA"

            pileup = bh.pileup(chrom, gPOS-2, gPOS+1, truncate=True, stepper="nofilter")
            base_counts = defaultdict(int)

            for pileupcolumn in pileup:
                # for snp
                if len(REF) == 1 and len(ALT) == 1:
                    if pileupcolumn.pos == gPOS-1:
                        for pileupread in pileupcolumn.pileups:
                            if pileupread.query_position is None:
                                continue
                            if pileupread.is_del:
                                continue
                            base = pileupread.alignment.query_sequence[pileupread.query_position]
                            base_counts[base] += 1
                
                elif len(REF) > len(ALT):
                    # for del 
                    if pileupcolumn.pos == gPOS-1:
                        for pileupread in pileupcolumn.pileups:
                            if "D" in pileupread.alignment.cigarstring:
                                base_counts[ALT] += 1
                            else:
                                base_counts[REF] += 1

                elif len(REF) < len(ALT):
                    # for ins 
                    if pileupcolumn.pos == gPOS-1:
                        for pileupread in pileupcolumn.pileups:
                            if "I" in pileupread.alignment.cigarstring:
                                base_counts[ALT] += 1
                            else:
                                base_counts[REF] += 1
            
            DP = sum(base_counts.values())
            AC = base_counts[ALT]

            if DP >= 1:
                AB = AC/DP
            else:
                AB = np.nan
            
            var_str = f'{gene}:{gPOS}:{REF}>{ALT},{transcript}:{pos},{annot}'
            MAT[var_str][f'{sample}-CLS{cluster}'] = f'{DP}|{AB:.2f}'
        
        bh.close()

    df = pd.DataFrame(MAT).T
    output_file_name = "data/shortread_10X_5UTR_singlecell_mecom/variants_of_interest.tsv"
    df.to_csv(output_file_name, sep = "\t")

    print(output_file_name)
    print(df)

