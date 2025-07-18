import pandas as pd
import pysam
from collections import defaultdict, Counter
import os
import logging
import pyfaidx
import argparse as parser
import HTSeq


rev = {"A" : "T", "T" : "A", "C" : "G", "G" : "C"}


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)

    parser = parser.ArgumentParser(description="Run psi calling on BAM files using pysam")
    parser.add_argument('-b', '--bedfile', type=str, required=True, help='Path to the bed file with exons of interest (BED format)')
    parser.add_argument('-o', '--output',   type=str, required=True, help='Output file name for the variants summary')
    parser.add_argument('-ov', '--ov', type=int, default=5, help='Overhang region size around the exon (default: 5)')
    parser.add_argument('-m', '--metadata', type=str, required=True, help='Path to the metadata file (TSV format)')
    args = parser.parse_args()


    regions  = pd.read_csv(args.bedfile, sep="\t")
    metadata = pd.read_csv(args.metadata , sep = "\t")


    logging.info(f"Starting psi calculation: {regions.shape[0]} regionss")
    
    region_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    for i, row in metadata.iterrows():
    
        bam_file, sample, data_type = row['BAM'], row['SM'], row['DATA_TYPE']

        if data_type == "lr":
            continue

        if os.path.exists(bam_file + ".bai") == False:
            logging.error(f"Index file not found for {bam_file}")
            continue

        try:
            bh = pysam.AlignmentFile(bam_file)
        except:
            logging.error(f"Error opening {bam_file}")
            continue


        for region in regions.itertuples():

            region_annot = f'{region.name}'
            region_iv = HTSeq.GenomicInterval(region.chrom, region.start, region.end)

            for read in bh.fetch(region.chrom, region.start - args.ov, region.end + args.ov):

                read_id = read.query_name
                RL = len(read.query_sequence)

                try:
                    umi = reads.get_tag('UB')
                    bc  = reads.get_tag('CB')
                except:
                    read_name, bc, umi = read_id.split("_")

                blocks = list(read.get_blocks())

                exons_ov, introns_ov = [], []

                for j in range(len(blocks)):

                    exon = HTSeq.GenomicInterval(region.chrom, blocks[j][0], blocks[j][1])

                    try: 
                        intron = HTSeq.GenomicInterval(region.chrom, blocks[j][1], blocks[j + 1][0])
                    except: 
                        intron = None
                        
                    exon_ov, intron_ov = 0, 0

                    # check which one overlaps more
                    if exon.overlaps(region_iv):
                        exon_ov = min(exon.end - region_iv.start, region_iv.end - exon.start)/region_iv.length
                    elif intron is not None and intron.overlaps(region_iv):
                        intron_ov = min(intron.end - region_iv.start, region_iv.end - intron.start)/region_iv.length

                    exons_ov.append(exon_ov)
                    introns_ov.append(intron_ov)

                if sum(exons_ov) > 0 or sum(introns_ov) > 0:
                    print(sample, read_id, exons_ov, introns_ov)

                if max(exons_ov) == 0 and max(introns_ov) == 0:
                    continue
                elif max(exons_ov) > max(introns_ov):
                    inc_or_exc = "inc"
                elif max(introns_ov) > max(exons_ov):
                    inc_or_exc = "exc"
                else:
                    print("Equal overlap", sample, read_id, exons_ov, introns_ov)
                
                region_dict[(sample, bc, region_annot)][data_type][umi].append(inc_or_exc)
            
        bh.close()
        logging.info(f"\tDone processing {sample} {i}")


    outline = defaultdict(dict)                

    for (sm, bc, region), val in region_dict.items():

        summary = {}

        for data_type, base_lst in val.items():

            base_counter = Counter()
            for umi, base_lst in base_lst.items():
                umi_base, count = Counter(base_lst).most_common(1)[0]
                base_counter[umi_base] += 1

            # I think because it's UMI, it's ok to not normalize by exon or read length
            # as it's typically done in short read data.
           
            DP = sum(base_counter.values())
            AC = base_counter['inc'] 
            PSI = base_counter['inc'] / DP if DP > 0 else np.nan

            outline[sm, bc][f"{region}_{data_type}_cov"] = AC 
            outline[sm, bc][f"{region}_{data_type}_psi"] = PSI

    table = pd.DataFrame(outline).T
    table = table.fillna('NA')
    table.to_csv(args.output , sep = "\t")

    print(table)
