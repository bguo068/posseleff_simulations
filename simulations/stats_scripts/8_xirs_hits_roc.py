#! /usr/bin/env python3
import pandas as pd

df = pd.read_parquet("./num_xirs_hits.pq")
df_orig = df.copy()

df["Model"] = df.Label.str.replace("_.*$", "", regex=True)
df["Sel"] = "neu"
df.loc[df.Label.str.contains("03"), "Sel"] = "s03"


lst = []
for nhits in [0, 1, 2, 3, 5, 10]:
    label = "Unfilt" if nhits == 0 else f"H{nhits}"

    s = (
        df[df.NumXirsHits >= nhits]
        .groupby(["Model", "Sel", "PeakFindMeth", "XirsCorrMeth"])["End"]
        .count()
        .rename(label)
    )

    lst.append(s)

df = pd.concat(lst, axis=1).fillna(0)
pd.options.display.float_format = "{:.1f}".format
print(df.round(1).style.format(precision=1).to_latex())

# Index(['Chromosome', 'Start', 'End', 'Median', 'Thres', 'GwChromStart',
#        'GwStart', 'GwEnd', 'NumXirsHits', 'GsId', 'Label', 'Rep',
#        'PeakFindMeth', 'XirsCorrMeth'],
#       dtype='object')


# single pop FN
_ = (420 - df.transpose()["sp_s03"]) / 420
print(_)

# single pop FP
_ = df.transpose()["sp_neu"] / 420
print(_)

# multi-pop FN
_ = (420 - df.transpose()["mp_s03"]) / 420
print(_)
