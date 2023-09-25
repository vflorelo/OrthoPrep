#!/usr/bin/python3
import sys
import pandas as pd
work_dir      = sys.argv[1]
blast_file    = sys.argv[2]
hq_blast_file = work_dir+"/hq/"+blast_file
mq_blast_file = work_dir+"/mq/"+blast_file
q_size_file   = sys.argv[3]
s_size_file   = sys.argv[4]
blast_file    = work_dir+"/"+blast_file
q_size_file   = work_dir+"/"+q_size_file
s_size_file   = work_dir+"/"+s_size_file
short_frac    = float(sys.argv[5])
long_frac     = float(sys.argv[6])
def filter_by_size(qlen,slen,sf,lf):
    if(qlen <= slen):
        low_bound = int(float(qlen*sf))
        up_bound  = int(float(slen*lf))
    else:
        low_bound = int(float(slen*sf))
        up_bound  = int(float(qlen*lf))
    low_frac = abs(qlen - slen) / low_bound
    up_frac  = abs(qlen - slen) / up_bound
    if ((low_frac   <= 1 ) and (up_frac <= 1)):
        pair_cat = "hq"
    elif ((low_frac >= 1 ) and (up_frac <= 1)):
        pair_cat = "mq"
    else:
        pair_cat = "lq"
    return pair_cat
blast_columns        = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
q_size_df            = pd.read_csv(q_size_file,header=None,sep="\t")
q_size_df.columns    = ["qseqid","qlen"]
s_size_df            = pd.read_csv(s_size_file,header=None,sep="\t")
s_size_df.columns    = ["sseqid","slen"]
blast_df             = pd.read_csv(blast_file,header=None,sep="\t")
blast_df.columns     = blast_columns
blast_df             = pd.merge(blast_df,q_size_df,on="qseqid")
blast_df             = pd.merge(blast_df,s_size_df,on="sseqid")
blast_df["pair_cat"] = blast_df.apply(lambda x : filter_by_size(x["qlen"],x["slen"],short_frac,long_frac), axis=1)
hq_blast_df          = blast_df.copy()
hq_blast_df          = hq_blast_df[hq_blast_df["pair_cat"]=="hq"]
hq_blast_df          = hq_blast_df[blast_columns]
hq_blast_df.to_csv(hq_blast_file,sep="\t",header=None,index=False,compression='gzip')
mq_blast_df          = blast_df.copy()
mq_blast_df          = mq_blast_df[mq_blast_df["pair_cat"]!="lq"]
mq_blast_df          = mq_blast_df[blast_columns]
mq_blast_df.to_csv(mq_blast_file,sep="\t",header=None,index=False,compression='gzip')