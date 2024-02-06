#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   plugins.py
@Time    :   2023/12/14 14:02:56
@Author  :   biolxy
@Version :   1.0
@Contact :   biolxy@aliyun.com
@Desc    :   None
'''

# here put the import lib
import sys
import os
import pandas as pd
from typing import List
# from icecream import ic
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# 生成 野生型 fasta   WILDTYPE


def ConstructWildtypeFa(infa, outfa):
    # /home/report/lixy/neo/tc1/fastaFiles/test1.tmp.8.fasta
    """
    fasta format:
        :>line7;NM_001099852;c.G247A;p.D83N;protein-altering;;(position;83;changed;from;D;to;N)
        :ETFQAVLNGLDALLT
    Args:
        infa (_type_): _description_
        outfa (_type_): _description_
    Add:
        vcf_manipulate.py line 191
    """
    records = []
    for seq_record in SeqIO.parse(infa, "fasta"):
        seqid = str(seq_record.id)
        tag = seqid.split(";")
        aa_wildType = tag[-3]
        aa_mutant = tag[-1]
        half_len = int(len(seq_record.seq) / 2)
        # ic(seq_record.id)
        # ic(seq_record.seq)
        # ic(aa_wildType, aa_mutant, seq_record.seq[half_len])
        # 赋值
        tag[4] = "wildType"
        nid = ";".join(tag)
        seq = str(seq_record.seq)
        seq = seq[:half_len] + aa_wildType + seq[half_len+1:]
        tmp_rec = SeqRecord(Seq(seq), id=nid, description="")
        records.append(tmp_rec)
        # ic(seq_record.seq)
        # ic(tmp_rec.seq)
    SeqIO.write(records, outfa, "fasta")


def getdf(infile):
    """_summary_

    Args:
        infile (_type_): 读取class1 的 df 表格 TestRun.MT.neoantigens.txt

    Returns:
        _type_: _description_
    """
    columns = ["Sample", "R1", "R2", "Line", "chr", "allelepos", "ref", "alt", "GeneName:RefSeqID", "pos", "hla", "peptide", "Core", "Of",
               "Gp", "Gl", "Ip", "Il", "Icore", "Identity", "Score_EL", "%Rank_EL", "Score_BA", "%Rank_BA", "Aff(nM)", "Candidate", "BindLevel"]
    df = pd.read_csv(infile, sep="\t", names=columns)
    df["peplen"] = len(df["peptide"].astype(str))
    df = df.drop(["R1", "R2", "Core", "Of", "Gp",
                 "Gl", "Ip", "Il", "Icore"], axis=1)
    return df


def getdf_class2(infile):
    """_summary_

    Args:
        infile (_type_): 读取class2 的 df 表格 TestRun.MT.neoantigens.txt

    ../netMHCIIpan -f example.pep -inptype 1 -a DRB1_0101 -BA

    INFO: Running Epitope Predictions for test2 on epitopes of length 13

    Returns:
        _type_: _description_
    """
    columns = ["Sample", "R1", "R2", "Line", "chr", "allelepos", "ref", "alt", "GeneName:RefSeqID", "pos", "hla", "peptide", "Of",
               "Core", "Core_Rel", "Identity", "Score_EL", "%Rank_EL", "Exp_Bind", "Score_BA", "Aff(nM)", "%Rank_BA", "BindLevel"]
    df = pd.read_csv(infile, sep="\t", names=columns)
    df["peplen"] = len(df["peptide"].astype(str))
    df = df.drop(["R1", "R2", "Core", "Of", "Core_Rel",
                 "Exp_Bind"], axis=1)
    return df


def mergeMTWT(infile1, infile2, outfile):
    """
    infile1: MT
    infile2: WT
    """
    df1 = getdf(infile1)
    df2 = getdf(infile2)
    on_columns = ["Sample", "Line", "chr", "allelepos",
                  "ref", "alt", "GeneName:RefSeqID", "pos", "hla", "Identity"]  # , "peplen"
    df = pd.merge(df1, df2, on=on_columns,
                  how='inner', suffixes=('_MT', '_WT'))
    df["Aff(nM)_FC"] = df["Aff(nM)_WT"] / df["Aff(nM)_MT"]
    # df.to_excel(outfile, index = True)
    df.to_csv(outfile, encoding='utf-8', sep='\t', index=False)
    return df


def mergeMTWT_class2(infile1, infile2, outfile):
    """
    infile1: MT
    infile2: WT
    """
    df1 = getdf_class2(infile1)
    df2 = getdf_class2(infile2)
    on_columns = ["Sample", "Line", "chr", "allelepos",
                  "ref", "alt", "GeneName:RefSeqID", "pos", "hla", "Identity"]  # , "peplen"
    df = pd.merge(df1, df2, on=on_columns,
                  how='inner', suffixes=('_MT', '_WT'))
    df["Aff(nM)_FC"] = df["Aff(nM)_WT"] / df["Aff(nM)_MT"]
    # df.to_excel(outfile, index = True)
    df.to_csv(outfile, encoding='utf-8', sep='\t', index=False)
    return df

def main():
    # infile = "/home/report/lixy/neo/tc1/fastaFiles/test1.tmp.9.fasta"
    # outfile = "xx.fasta"
    # ConstructWildtypeFa(infile, outfile)

    # infile1 = "/home/report/lixy/neo/tc1/TestRun.MT.neoantigens.txt"
    # infile2 = "/home/report/lixy/neo/tc1/TestRun.WT.neoantigens.txt"
    # outfile = 'aa.txt'
    # df = mergeMTWT(infile1, infile2, outfile)
    # ic(df)


    infile1 = "/home/report/lixy/neo/tc2/TestRun.MT.neoantigens.txt"
    infile2 = "/home/report/lixy/neo/tc2/TestRun.WT.neoantigens.txt"
    outfile = 'aa2.txt'

    df = mergeMTWT_class2(infile1, infile2, outfile)
    # ic(df)

    pass


if __name__ == "__main__":
    main()
