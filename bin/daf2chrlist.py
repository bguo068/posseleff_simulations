#! /usr/bin/env python3
import argparse
import pandas as pd


p = argparse.ArgumentParser()
p.add_argument("--daf", type=str, required=True)
p.add_argument("--min_afreq", type=float, required=True)
args = p.parse_args()
daf = pd.read_csv(args.daf, sep="\t")
# print(daf.head())

cols = []
if daf.shape[0] == 0:
    for colname in daf.columns:
        if colname == "GEN":
            continue
        else:
            cols.append(colname)
else:
    for colname in daf.columns:
        if colname == "GEN":
            continue
        if daf[colname].iat[-1] > args.min_afreq:
            cols.append(colname)

cols = [n.replace("DAF_CHR", "") for n in cols]
print(" ".join(cols))
