import pandas as pd
import pysam
from collections import defaultdict, Counter
import os
import logging
import pyfaidx
import argparse as parser


rev = {"A" : "T", "T" : "A", "C" : "G", "G" : "C"}


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)

    parser = parser.ArgumentParser(description="Run variant calling on BAM files using pysam")
    parser.add_argument('-v', '--variants', type=str, required=True, help='Path to the variants file (TSV format)')
    parser.add_argument('-o', '--output',   type=str, required=True, help='Output file name for the variants summary')
    parser.add_argument('-m', '--metadata', type=str, required=True, help='Path to the metadata file (TSV format)')
    args = parser.parse_args()


    VARS = pd.read_csv(args.variants, sep="\t")
    metadata = pd.read_csv(args.metadata , sep = "\t")

    logging.info(f"Starting variant calling: {VARS.shape[0]} variants")
    
    var_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    for i, row in metadata.iterrows():
    
        bam_file, sample, data_type = row['BAM'], row['SM'], row['DATA_TYPE']

        if os.path.exists(bam_file + ".bai") == False:
            logging.error(f"Index file not found for {bam_file}")
            continue

        try:
            bh = pysam.AlignmentFile(bam_file)
        except:
            logging.error(f"Error opening {bam_file}")
            continue


        for var in VARS.itertuples():

            var_annot = f'{var.chrom}_{var.pos}_{var.REF}_{var.ALT}_{var.mut}'
            
            pileup = bh.pileup(var.chrom, var.pos - 2, var.pos + 1, truncate=True, stepper="nofilter")

            for pileupcolumn in pileup:
                # for snp
                if pileupcolumn.pos == var.pos - 1:
                    for pileupread in pileupcolumn.pileups:

                        read_id = pileupread.alignment.query_name

                        try:
                            umi = pileupread.alignment.get_tag('UB')
                            bc  = pileupread.alignment.get_tag('CB')
                        except:
                            read_name, bc, umi = read_id.split("_")

                        if len(var.REF) == 1 and len(var.ALT) == 1:
                            if pileupread.query_position is None:
                                continue
                            if pileupread.is_del:
                                continue
                            base = pileupread.alignment.query_sequence[pileupread.query_position]

                        elif len(var.REF) > len(var.ALT):
                            # for del 
                            if pileupread.indel == len(var.ALT) - len(var.REF):
                                base = var.ALT
                            else:
                                base = var.REF

                        elif len(var.REF) < len(var.ALT):
                            # for ins 
                            if pileupread.indel == len(var.ALT) - len(var.REF):
                                base = var.ALT
                            else:
                                base = var.REF
                
                        var_dict[(sample, bc, var_annot)][data_type][umi].append(base)
            
        bh.close()
        logging.info(f"\tDone processing {sample} {i}")


    outline = defaultdict(dict)                

    for (sm, bc, var), val in var_dict.items():

        summary = {}

        CHROM, POS, REF, ALT = var.split("_")[:4]

        for data_type, base_lst in val.items():

            base_counter = Counter()
            for umi, base_lst in base_lst.items():
                umi_base, count = Counter(base_lst).most_common(1)[0]
                base_counter[umi_base] += 1

            DP = sum(base_counter.values())
            AB = base_counter[ALT] / DP

            print(sm, bc, var, Counter(base_lst), DP, AB)

            outline[sm, bc][f"{var}_{data_type}_cov"]   = DP #(ab_lr + ab_sr)
            outline[sm, bc][f"{var}_{data_type}_ratio"] = AB

    table = pd.DataFrame(outline).T
    table = table.fillna('NA')
    table.to_csv(args.output , sep = "\t")

    print(table)
