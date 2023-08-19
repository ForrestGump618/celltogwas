import pandas as pd
import numpy as np
import gzip
import re
import time
import bisect

def extract_gene_info(gtf_file):
    time0 = time.time()
    data = [["GENE","CHR","START","END"]]
    with gzip.open(gtf_file, 'rt') as file:
        i = 0
        for line in file:
            line = line.strip()
            if not line.startswith('#'):
                columns = line.split('\t')
                if len(columns) >= 9:
                    gene_info = columns[8]
                    match = re.search(r'gene_name "(.*?)"; transcript_type', gene_info)
                    if match:
                        i += 1
                        if i%500000 == 0:
                            print(f"{i} items completed!")
                        gene_name = match.group(1)
                        data.append([gene_name,columns[0], columns[3], columns[4]])
    print(f"runing time: {time.time()-time0}s")
    return pd.DataFrame(data[1:], columns=data[0])
                        
def merge_intervals(intervals):
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    merged = []
    for interval in sorted_intervals:
        if not merged or interval[0] > merged[-1][1]:
            merged.append(interval)
        else:
            merged[-1] = [merged[-1][0], max(merged[-1][1], interval[1])]
    merged = sorted(merged, key=lambda x: x[0])
    return merged

def celltype_genome_range(cell_gene_dict, gene_info, windowsize):
    time0 = time.time()
    celltype_range = {}
    
    for celltype, genes in cell_gene_dict.items():
        gene_set = pd.DataFrame(genes, columns=["GENE"])
        df = pd.merge(gene_set, gene_info, on='GENE', how='inner')
        df['START'] = np.maximum(1, df['START'].astype(int) - windowsize)
        df['END'] = df['END'].astype(int) + windowsize
        
        iter_df_pre = [['chr' + str(x1).lstrip('chr'), x2 - 1, x3] for x1, x2, x3 in np.array(df[['CHR', 'START', 'END']])]
        iter_df_pre_dic = {}
        
        for sublist in iter_df_pre:
            key = sublist[0]
            value = sublist[1:]
            iter_df_pre_dic.setdefault(key, []).append(value)
        
        celltype_range_sub = {}
        
        for chrom, intervals in iter_df_pre_dic.items():
            range_list = merge_intervals(intervals)
            element = [[range_list[i][0], range_list[i][1]] for i in range(len(range_list))]
            celltype_range_sub[chrom] = element
        
        print(f"{celltype} completed！")
        celltype_range[celltype] = celltype_range_sub
    
    print(f"running：{time.time() - time0}s")
    return celltype_range


def check_value_in_intervals(intervals, value):
    left_boundaries = [interval[0] for interval in intervals]
    index = bisect.bisect_left(left_boundaries, value)
    
    if index < len(intervals):
        interval = intervals[index]
        if interval[0] <= value <= interval[1]:
            return True
    
    return False

def gwas_annotation(gwas, celltype_range):
    time0 = time.time()
    chr_list = set(gwas["CHR"].tolist())
    all_annotation = {}
    for celltype, annotation_ref in celltype_range.items():
        sub_gwas_annotation = []
        for chrom in chr_list:
            sub_gwas = gwas[gwas["CHR"] == chrom].values.tolist()
            sub_annotation_ref = annotation_ref[chrom]
            for item in sub_gwas:
                if check_value_in_intervals(sub_annotation_ref, item[1]):
                    item.append("1")
                else:
                    item.append("0")
                sub_gwas_annotation.append(item)
            print(f"{celltype} {chrom}: completed!",sep = " ")

        all_annotation[celltype] = sub_gwas_annotation
        print(f"{celltype} completed!")
        print(f"running time: {time.time() - time0}s")
        time0 = time.time()

    return all_annotation