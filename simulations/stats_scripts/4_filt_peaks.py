#! /usr/bin/env python3

from ibdutils.utils.ibdutils import IBD
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybedtools as pb


params = [
    dict(gs_id=10003, label="sp_s03", rep=0),
    dict(gs_id=10000, label="sp_neu", rep=0),
    dict(gs_id=20003, label="mp_s03", rep=0),
    dict(gs_id=20000, label="mp_s00", rep=0),
]


def proc_data(params_id: int):
    gs_id = params[params_id]["gs_id"]
    label = params[params_id]["label"]
    rep = params[params_id]["rep"]

    if label.startswith("sp"):
        fn = f"resdir/{gs_id}_{label}_rep{rep}/ibdne_ibd/{gs_id}_orig.ibdne.ibdobj.gz"
    else:
        fn = f"resdir/{gs_id}_{label}_rep{rep}/ifm_input/{gs_id}_orig.ifm.ibdobj.gz"

    ibd = IBD.pickle_load(fn)

    return ibd


res = []

for i in range(0, len(params)):
    res_i = proc_data(i)
    print(i)
    res.append(res_i)

fig, axes = plt.subplots(nrows=6, constrained_layout=True)
ibd = res[0]
xirs = ibd._xirs_df.copy()
ibd.plot_coverage(axes[0], label="s=0.3, unfilt", which="unfilt", plot_proportions=True)
ibd.plot_coverage(
    axes[1], label="s=0.3, xirsfilt", which="xirsfilt", plot_proportions=True
)
ibd.plot_xirs(axes[2], label="s=0.3")
ibd = res[1]
xirs = ibd._xirs_df.copy()
ibd.plot_coverage(axes[3], label="s=0.0, unfilt", which="unfilt", plot_proportions=True)
ibd.plot_coverage(
    axes[4], label="s=0.0, xirsfilt", which="xirsfilt", plot_proportions=True
)
ibd.plot_xirs(axes[5])
# axes[5].sharey(axes[2])
fig.set_size_inches(11, 7)
fig.savefig("sim_xirs_single_pop_all.png", dpi=600)


fig, axes = plt.subplots(nrows=6, constrained_layout=True)
ibd = res[2]
xirs = ibd._xirs_df.copy()
ibd.plot_coverage(axes[0], label="s=0.3, unfilt", which="unfilt", plot_proportions=True)
ibd.plot_coverage(
    axes[1], label="s=0.3, xirsfilt", which="xirsfilt", plot_proportions=True
)
ibd.plot_xirs(axes[2], label="s=0.3")
ibd = res[3]
xirs = ibd._xirs_df.copy()
ibd.plot_coverage(axes[3], label="s=0.0, unfilt", which="unfilt", plot_proportions=True)
ibd.plot_coverage(
    axes[4], label="s=0.0, xirsfilt", which="xirsfilt", plot_proportions=True
)
ibd.plot_xirs(axes[5])
# axes[5].sharey(axes[2])
fig.set_size_inches(11, 7)
fig.savefig("sim_xirs_multi_pop_all.png", dpi=600)
