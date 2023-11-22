#!/usr/bin/python3
import sys
import pandas as pd
prep_dir   = sys.argv[1]
query      = sys.argv[2]
subject    = sys.argv[3]
len_file   = sys.argv[4]
short_frac = float(sys.argv[5])
long_frac  = float(sys.argv[6])
blast_file = prep_dir+"/Blast"+query+"_"+subject+".txt.gz"
len_file   = prep_dir + "/" + len_file
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
    elif ((low_frac >= 1 ) and (up_frac <= 1 )):
        pair_cat = "mq"
    else:
        pair_cat = "lq"
    return pair_cat
sequence_len_df = pd.read_csv(len_file,header=None,sep="\t")
sequence_len_df.columns = ["species","seqid","len"]
blast_columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
q_len_df = sequence_len_df.copy()
q_len_df = q_len_df[q_len_df["species"]==query]
q_len_df.rename(columns={"seqid":"qseqid"})
q_len_df.drop(columns=['species'])
s_len_df = sequence_len_df.copy()
s_len_df = s_len_df[s_len_df["species"]==subject]
s_len_df.rename(columns={"seqid":"qseqid"})
s_len_df.drop(columns=['species'])
op_blast_file = prep_dir+"/WorkingDirectory/Blast"+query+"_"+subject+".txt.gz"
blast_df   = pd.read_csv(blast_file,header=None,sep="\t")
blast_df.columns = blast_columns
blast_df = pd.merge(blast_df,q_len_df,on="qseqid")
blast_df = pd.merge(blast_df,s_len_df,on="sseqid")
blast_df["pair_cat"] = blast_df.apply(lambda x : filter_by_len(x["qlen"],x["slen"],short_frac,long_frac), axis=1)
op_blast_df = blast_df.copy()
op_blast_df = op_blast_df[op_blast_df["pair_cat"]=="hq"]
op_blast_df = op_blast_df[blast_columns]
op_blast_df.to_csv(op_blast_file,sep="\t",header=None,index=False,compression='gzip')