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

def celltype_genome_range(cell_gene_dict,gene_info,windowsize):
    time0 = time.time()
    # input file: a dictionary: 
                    # key(str): name of celltype;
                    # value(list):list of genes' names.
    celltype_range = {}
    for celltype in list(cell_gene_dict.keys()):
        gene_set = pd.DataFrame(cell_gene_dict[celltype],columns=["GENE"])
        df = pd.merge(gene_set, gene_info, on = 'GENE', how = 'inner')
        df['START'] = np.maximum(1, int(df['START']) - windowsize)
        df['END'] = df['END'] + windowsize
        iter_df_pre = [['chr'+(str(x1).lstrip('chr')), x2 - 1, x3] for (x1,x2,x3) in np.array(df[['CHR', 'START', 'END']])]
        iter_df_pre_dic = {}
        for sublist in iter_df_pre:
            key = sublist[0]
            value = sublist[1:]
            if key in iter_df_pre_dic:
                iter_df_pre_dic[key].append(value)
            else:
                iter_df_pre_dic[key] = [value]
        celltype_range_sub = {}
        for chrom in list(iter_df_pre_dic.keys()):
            range_list = merge_intervals(iter_df_pre_dic[chrom])
            element = [[range_list[i][0],range_list[i][1]] for i in range(len(range_list))]
            celltype_range_sub[chrom] = element
        print(f"{celltype} completed!")
        celltype_range[celltype] = element
    print(f"runing time: {time.time()-time0}s")
    return celltype_range


def check_value_in_intervals(intervals, value):
    left_boundaries = [interval[0] for interval in intervals]
    index = bisect.bisect_left(left_boundaries, value)
    
    if index < len(intervals):
        interval = intervals[index]
        if interval[0] <= value <= interval[1]:
            return True
    
    return False

def gwas_annotaion(gwas,celltype_range):
    time0 = time.time()
    # gwas file must contain:chr("chr*"), pos(int), and it should be the first two cols.
    # option: rsid and so on
    chr_list = list(set(gwas["chr"].tolist()))
    all_annotation = {}
    for celltype in list(celltype_range.keys()):
        annotation_ref = celltype_range[celltype]
        sub_gwas_annotation = []
        for chrom in chr_list:
            sub_gwas = gwas[gwas["chr"==chrom]].values().tolist()
            sub_annotation_ref = annotation_ref[chrom]
            for item in sub_gwas:
                if check_value_in_intervals(sub_annotation_ref, item[1]):
                    item.append("1")
                else:
                    item.append("0")
                sub_gwas_annotation.append([chrom]+item)
        all_annotation[celltype] = sub_gwas_annotation
    print(f"runing time: {time.time()-time0}s")
    return all_annotation