import pandas as pd
import numpy as np
import csv

def load_fa(path):
    """a function to read fasta file from the path and store in a dict"""
    genes_seq = {}  #将序列存入字典
    with open(path,"r") as sequences:  #以读取方式打开文件
        lines = sequences.readlines()


    for line in lines:
        if line.startswith(">"):
            genename = line.split()[0]  #这个地方需要灵活调整
            genes_seq[genename] = ''  #序列为字符串
        else:
            genes_seq[genename] += line.strip()
    return genes_seq

def build_cs(seq):
    """a function to calculate cs from seq"""
    cs = []  # cs存储在列表中
    n_cs = len(seq)

    for i in range(n_cs):
        c = seq[i]
        cs.append(c)

    return cs


from collections import Counter
def summary_ACGT(cs):
    """a function to summarize the kmers"""
    ACGT_stat = dict(Counter(cs))
    return ACGT_stat

genes_seq = load_fa(path="../data/DNA_f/phage.fasta")

genes_cs = {}
for gene in genes_seq.keys():
  genes_cs[gene] = summary_ACGT(build_cs(seq=genes_seq[gene]))

cs_stat_frame = pd.DataFrame(genes_cs)
cs_freq = lambda x : x/np.sum(x)
cs_freq_frame =cs_stat_frame.apply(cs_freq, axis=0)


cs_freq_frame = pd.DataFrame(cs_freq_frame.values.T, index = cs_freq_frame.columns, columns = cs_freq_frame.index)

print(cs_freq_frame.shape)
#print(cs_freq_frame.apply(lambda x:x.sum(),axis =1))
#print(cs_freq_frame)

cs_freq_frame.to_csv("../data/DNA_f/onlycs/phage_cs.csv",index = True, sep=',',header = True,na_rep='0')





