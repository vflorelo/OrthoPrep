#!/usr/bin/python3
import sys
import pandas as pd
prep_dir   = sys.argv[1]
query      = sys.argv[2]
subject    = sys.argv[3]
len_file   = sys.argv[4]
len_file   = prep_dir + "/" + len_file
len_df     = pd.read_csv(len_file,header=None,sep="\t",names=["species","seqid","len"])
len_types  = {"species":str,
              "seqid"  :str,
              "len"    :int}
len_df     = len_df.astype(len_types)
short_frac = float(sys.argv[5])
long_frac  = float(sys.argv[6])
blast_file = prep_dir+"/Blast"+query+"_"+subject+".txt.gz"
blast_df   = pd.read_csv(blast_file,
                         header=None,
                         sep="\t",
                         names = ["qseqid","sseqid",  "pident",
                                  "length","mismatch","gapopen",
                                  "qstart","qend",    "sstart",
                                  "send",  "evalue",  "bitscore"])
def filter_by_len(qlen,slen,sf,lf):
    if(qlen <= slen):
        low_bound = int(float(qlen*sf))
        up_bound  = int(float(slen*lf))
    else:
        low_bound = int(float(slen*sf))
        up_bound  = int(float(qlen*lf))
    low_frac = abs(qlen - slen) / low_bound
    up_frac  = abs(qlen - slen) / up_bound
    if ((low_frac   <= 1 ) and (up_frac <= 1 )):
        pair_cat = "hq"
    else:
        pair_cat = "lq"
    return pair_cat
q_len_df = len_df.copy()
q_len_df = q_len_df[q_len_df["species"]==query]
q_len_df = q_len_df.rename(columns={"seqid":"qseqid","len":"qlen"})
q_len_df = q_len_df.drop(columns=['species'],axis=1)
s_len_df = len_df.copy()
s_len_df = s_len_df[s_len_df["species"]==subject]
s_len_df = s_len_df.rename(columns={"seqid":"sseqid","len":"slen"})
s_len_df = s_len_df.drop(columns=['species'],axis=1)
blast_df = pd.merge(blast_df,q_len_df,on="qseqid")
blast_df = pd.merge(blast_df,s_len_df,on="sseqid")
blast_df["pair_cat"] = blast_df.apply(lambda x : filter_by_len(x["qlen"],x["slen"],short_frac,long_frac), axis=1)
op_blast_df = blast_df.copy()
op_blast_df = op_blast_df[op_blast_df["pair_cat"]=="hq"]
blast_columns = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"]
op_blast_df = op_blast_df[blast_columns]
op_blast_file = prep_dir+"/tmp/Blast"+query+"_"+subject+".txt.gz"
op_blast_df.to_csv(op_blast_file,sep="\t",header=None,index=False,compression='gzip')