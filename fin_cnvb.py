import pandas as pd
from intervaltree import IntervalTree

def load_cnv_data(filepath):
    df = pd.read_excel(filepath)
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    return df

def extract_cnvr_regions(df):
    cnvr_list = []
    region_intervals = {}
    cnvr_id = 1

    for chrom in df['chr'].unique():
        chrom_df = df[df['chr'] == chrom].sort_values(by='start')
        tree = IntervalTree()

        for idx, row in chrom_df.iterrows():
            tree[row['start']:row['end']] = {
                'sample': row['Sample'],
                'cn': row['cn'],
                'start': row['start'],
                'end': row['end']
            }

        tree.merge_overlaps(strict=False)
        merged = sorted(tree)

        region_intervals[chrom] = IntervalTree()

        for region in merged:
            overlapping = chrom_df[
                (chrom_df['start'] <= region.end) & (chrom_df['end'] >= region.begin)
            ]

            if overlapping.empty:
                continue

            unique_cns = sorted(overlapping['cn'].unique())
            if all(cn < 2 for cn in unique_cns):
                region_type = "DEL CNV Region"
            elif all(cn > 2 for cn in unique_cns):
                region_type = "DUP CNV Region"
            else:
                region_type = "Mixed CNV Region"

            cnvr_id_str = f'CNVR_{cnvr_id}'

            cnvr_list.append({
                'CNVR_ID': cnvr_id_str,
                'chr': chrom,
                'region_start': region.begin,
                'region_end': region.end,
                'num_cnvs': len(overlapping),
                'region_type': region_type,
                'unique_cn_values': ','.join(map(str, unique_cns)),
                'length': region.end - region.begin
            })

            region_intervals[chrom].addi(region.begin, region.end, cnvr_id_str)
            cnvr_id += 1

    cnvr_df = pd.DataFrame(cnvr_list)
    return cnvr_df, region_intervals

def extract_cnv_blocks(df, region_intervals):
    all_blocks = []

    for chrom in df['chr'].unique():
        chrom_df = df[df['chr'] == chrom]

        breakpoints = set(chrom_df['start']).union(set(chrom_df['end']))
        sorted_breaks = sorted(breakpoints)

        for i in range(len(sorted_breaks) - 1):
            start = sorted_breaks[i]
            end = sorted_breaks[i + 1]

            overlapping = chrom_df[
                (chrom_df['start'] <= start) & (chrom_df['end'] >= end)
            ]

            if overlapping.empty:
                continue

            samples = overlapping['Sample'].unique().tolist()
            cnvr_matches = region_intervals.get(chrom, IntervalTree()).overlap(start, end)
            cnvr_id = next(iter(cnvr_matches)).data if cnvr_matches else 'NA'

            all_blocks.append({
                'chr': chrom,
                'CNVR_ID': cnvr_id,
                'START': start,
                'END': end,
                'NUM_SAMPLES': len(samples),
                'SAMPLES': ','.join(map(str, samples)),
                'LENGTH': end - start
            })

    return pd.DataFrame(all_blocks)

def process_cnv_file(input_path):
    df = load_cnv_data(input_path)
    cnvr_df, region_intervals = extract_cnvr_regions(df)
    blocks_df = extract_cnv_blocks(df, region_intervals)

    blocks_df.to_excel('CNV_Blocks_Output.xlsx', index=False)
    cnvr_df.to_excel('CNVR_Summary_Output.xlsx', index=False)

    print("✅ Output files generated:")
    print("  → CNV_Blocks_Output.xlsx")
    print("  → CNVR_Summary_Output.xlsx")

# 🔁 Example call
# Replace this path with your input .xlsx file path:
process_cnv_file("/gpfs/data/user/yuvraj/cnv/sanscog_allbatch.goodcnv.merge.regions_filtered_segdup_filtered_conf10.xlsx")

