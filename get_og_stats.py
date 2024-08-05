#!/usr/bin/python3
import pandas as pd
import statistics
import sys
from   scipy import stats
in_tsv_file  = sys.argv[1]
out_tsv_file = sys.argv[2]
df = pd.read_csv(in_tsv_file,sep="\t")
df["len_list"] = df["lengths"].apply(lambda  x: list(map(int,x.split(","))))
df["median"]   = df["len_list"].apply(lambda x: statistics.median(x))
df["mad"]      = df["len_list"].apply(lambda x: stats.median_abs_deviation(x)   )
df["size"]     = df["len_list"].apply(lambda x: len(x))
col_list       = ["og_id","median","mad","size"]
df[col_list].to_csv(out_tsv_file,sep="\t",index=False,header=False)