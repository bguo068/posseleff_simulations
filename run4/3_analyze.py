#! /usr/bin/env python3

import igraph
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from scipy import stats

ms = pd.read_parquet("member.pq")
ne = pd.read_parquet("ne.pq")


def run_ne_stats(grp):
    sel = (ne.Gen <= 100) & (ne.Grp == grp) & (ne.OrigRmMode == "orig")
    x = ne[sel].pivot(columns="Rep", index=["Gen"], values="Ne").sort_index()
    x = np.log10(x)
    sel = (ne.Gen <= 100) & (ne.Grp == grp) & (ne.OrigRmMode == "rmpeaks")
    y = ne[sel].pivot(columns="Rep", index=["Gen"], values="Ne").sort_index()
    y = np.log10(y)
    res1 = stats.wilcoxon(x, y, axis=1)
    res2 = stats.ttest_rel(x, y, axis=1)
    x_mean = x.mean(axis=1)
    y_mean = y.mean(axis=1)
    x_l95, x_u95 = stats.t.interval(
        0.95, x.shape[1] - 1, loc=x_mean, scale=stats.sem(x, axis=1)
    )
    y_l95, y_u95 = stats.t.interval(
        0.95, y.shape[1] - 1, loc=y_mean, scale=stats.sem(y, axis=1)
    )
    non_overlapp = (x_u95 < y_l95) | (x_l95 > y_u95)
    M = np.vstack(
        [
            np.power(10, x_mean),
            np.power(10, x_l95),
            np.power(10, x_u95),
            np.power(10, y_mean),
            np.power(10, y_l95),
            np.power(10, y_u95),
            res1.pvalue,
            res2.pvalue,
        ]
    )

    out = pd.DataFrame(M.transpose())
    out.columns = [
        "Orig",
        "OrigL95",
        "OrigU95",
        "Rmpeaks",
        "RmpeaksL95",
        "RmpeaksU95",
        "PvalueSignedRank",
        "PvaluePairedtest",
    ]
    out["Ci95ONotOverlap"] = non_overlapp
    return out


def plot_ne(ne, ne_stats_dict, grp, ax1, ax2, ax3):
    sel = (ne.Gen <= 100) & (ne.Grp == grp) & (ne.OrigRmMode == "orig")
    x = ne[sel].pivot(columns="Rep", index=["Gen"], values="Ne")
    sel = (ne.Gen <= 100) & (ne.Grp == grp) & (ne.OrigRmMode == "rmpeaks")
    y = ne[sel].pivot(columns="Rep", index=["Gen"], values="Ne")

    ncols = x.shape[1]
    for col in range(ncols):
        xx = x.iloc[:, col]
        yy = y.iloc[:, col]
        gen = xx.index
        label = "orig" if col == 0 else None
        ax1.plot(gen, xx, alpha=0.05, color="red", label=label)
        label = "rmpeaks" if col == 0 else ""
        ax1.plot(gen, yy, alpha=0.05, color="blue", label=label)
        ax1.set_xlim(0, 100)
        ax1.set_ylim(3e2, 3e4)
        ax1.set_yscale("log")
    ax1.plot(gen, x.median(axis=1), color="red", lw=2, linestyle="--")
    ax1.plot(gen, y.median(axis=1), color="blue", lw=2, linestyle="--")
    ax1.set_title(grp)
    ax1.set_ylabel("Ne")
    ax1.set_xlabel("generations ago")
    ax1.legend()

    stats = ne_stats_dict[grp]
    ax2.plot(
        stats.index, -np.log10(stats.PvaluePairedtest), "k-", label="paired-t test"
    )
    ax2.plot(
        stats.index, -np.log10(stats.PvalueSignedRank), "k+", label="signed rank test"
    )
    ax2.axhline(y=-np.log10(0.05 / x.shape[0]), color="k", linestyle="--")
    ax2.set_ylabel("-log10(p)")
    ax2.set_xlim(0, 100)
    ax2.legend()

    ax3.plot(stats.index, np.median(x - y, axis=1), "k*")
    ax3.axhline(y=0, color="k", linestyle="--")
    ax3.set_ylabel("median of difference")
    ax3.set_xlim(0, 100)


def calc_ms_adj_rand():
    grp_lst = []
    rep_lst = []
    mode_lst = []
    adjrand_lst = []
    for grp in ms.Grp.unique():
        for rep in ms.Rep.unique():
            for mode in ms.OrigRmMode.unique():
                sel = (ms.Grp == grp) & (ms.Rep == rep) & (ms.OrigRmMode == mode)
                df = ms[sel]
                true = df.Population.to_list()
                infer = df.Rank.to_list()
                adjrand = igraph.compare_communities(
                    true, infer, method="adjusted_rand"
                )
                grp_lst.append(grp)
                rep_lst.append(rep)
                mode_lst.append(mode)
                adjrand_lst.append(adjrand)
    df = pd.DataFrame(
        {"Grp": grp_lst, "Rep": rep_lst, "OrigRmMode": mode_lst, "AdjRand": adjrand_lst}
    )
    return df


def run_ms_stats(ms_adjrand):
    sel = ms_adjrand.OrigRmMode == "orig"
    x = ms_adjrand[sel].pivot(columns="Rep", values="AdjRand", index="Grp").sort_index()
    sel = ms_adjrand.OrigRmMode == "rmpeaks"
    y = ms_adjrand[sel].pivot(columns="Rep", values="AdjRand", index="Grp").sort_index()

    return pd.DataFrame(
        {
            "OrigMean": x.mean(axis=1),
            "OrigStd": x.std(axis=1),
            "RmpeaskMean": y.mean(axis=1),
            "RmpeaskStd": y.std(axis=1),
            "PvaluePairedtest": stats.ttest_rel(x, y, axis=1).pvalue,
            "PvalueTtest": stats.ttest_ind(x, y, axis=1, equal_var=False).pvalue,
        }
    )


sp_grps = ["sp_neu", "sp_s01", "sp_s02", "sp_s03"]
mp_grps = ["mp_s00", "mp_s01", "mp_s02", "mp_s03"]

ne_stats_dict = {}
with pd.ExcelWriter("ne_summary.xlsx") as writer:
    for grp in sp_grps:
        df = run_ne_stats(grp)
        ne_stats_dict[grp] = df
        df.to_excel(writer, sheet_name=grp)

fig, axes = plt.subplots(
    ncols=len(sp_grps), nrows=3, constrained_layout=True, sharey="row", figsize=(11, 6)
)

for i, grp in enumerate(sp_grps):
    ax1 = axes[0, i]
    ax2 = axes[1, i]
    ax3 = axes[2, i]
    plot_ne(ne, ne_stats_dict, grp, ax1, ax2, ax3)

fig.savefig("sim_ne.png", dpi=600)

ms_adjrand = calc_ms_adj_rand()
ms_stats = run_ms_stats(ms_adjrand)
ms_stats.to_excel("membership_summary.xlsx")

with open("membership_summary.tex", "w") as f:
    s = ms_stats.style.to_latex()
    f.write(s)
