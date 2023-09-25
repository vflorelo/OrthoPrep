#!/usr/bin/python3
import sys
import pandas as pd
work_dir   = sys.argv[1]
short_frac = float(sys.argv[2])
long_frac  = float(sys.argv[3])
def filter_by_size(qlen,slen,sf,lf):
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
species_list_file = work_dir + "/SpeciesIDs.txt"
with open(species_list_file, "r") as file:
    species_list = file.read().split('\n')
species_num = [species.split(":")[0] for species in species_list]
species_num = list(filter(None, species_num))
blast_columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
for query in species_num:
    q_size_file = work_dir+"/Species"+query+".sizes.tsv"
    q_size_df = pd.read_csv(q_size_file,header=None,sep="\t")
    q_size_df.columns = ["qseqid","qlen"]
    for subject in species_num:
        s_size_file = work_dir+"/Species"+subject+".sizes.tsv"
        s_size_df = pd.read_csv(s_size_file,header=None,sep="\t")
        s_size_df.columns = ["sseqid","slen"]
        blast_file = work_dir+"/Blast"+query+"_"+subject+".txt.gz"
        hq_blast_file = work_dir+"/hq/Blast"+query+"_"+subject+".txt.gz"
        mq_blast_file = work_dir+"/mq/Blast"+query+"_"+subject+".txt.gz"
        blast_df   = pd.read_csv(blast_file,header=None,sep="\t")
        blast_df.columns = blast_columns
        blast_df = pd.merge(blast_df,q_size_df,on="qseqid")
        blast_df = pd.merge(blast_df,s_size_df,on="sseqid")
        blast_df["pair_cat"] = blast_df.apply(lambda x : filter_by_size(x["qlen"],x["slen"],short_frac,long_frac), axis=1)
        hq_blast_df = blast_df.copy()
        hq_blast_df = hq_blast_df[hq_blast_df["pair_cat"]=="hq"]
        hq_blast_df = hq_blast_df[blast_columns]
        hq_blast_df.to_csv(hq_blast_file,sep="\t",header=None,index=False,compression='gzip')
        mq_blast_df = blast_df.copy()
        mq_blast_df = mq_blast_df[mq_blast_df["pair_cat"]!="lq"]
        mq_blast_df = mq_blast_df[blast_columns]
        mq_blast_df.to_csv(mq_blast_file,sep="\t",header=None,index=False,compression='gzip')