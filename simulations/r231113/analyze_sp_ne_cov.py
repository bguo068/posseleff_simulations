import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ibdutils.utils.ibdutils import IBD

RESDIR: str = "./resdir_sp"
SIMULATIONS = json.loads(Path("./genome_sets_sp.json").read_text())
S_LIST = [0.0, 0.3]
SELG_LIST = [50, 80]
DELTA_LIST = [0.5, 0.3, 0.1, 0.01, 0.001]
RELG_LIST = [100, 50, 25]
POWER_LIST = [8, 16]
N_DELTA = len(DELTA_LIST)
N_RELG = len(RELG_LIST)
N_POWER = len(POWER_LIST)
N_S = len(S_LIST)
N_SELG = len(SELG_LIST)


def get_ne_path(resdir, label, gsid, rmmode):
    p = f"{resdir}/{gsid}_{label}/ne_output/{gsid}_2.0_10.0_none_{rmmode}.ne"
    assert Path(p).exists(), f"{p} does not exist"
    return p


def get_ne_path_rep(resdir, label, gsid, repno, rmmode):
    gsid = gsid + 10 * repno
    p = f"{resdir}/{gsid}_{label}_rep{repno}/ne_output/{gsid}_2.0_10.0_none_{rmmode}.ne"
    assert Path(p).exists(), f"{p} does not exist"
    return p


def get_true_ne_path(resdir, label, gsid):
    p = f"{resdir}/{gsid}_{label}/true_ne/{gsid}_1.true_ne"
    assert Path(p).exists(), f"{p} does not exist"
    return p


def get_true_ne_path_rep(resdir, label, gsid, repno):
    gsid = gsid + 10 * repno
    p = f"{resdir}/{gsid}_{label}_rep{repno}/true_ne/{gsid}_1.true_ne"
    assert Path(p).exists(), f"{p} does not exist"
    return p


def read_ne(ne_file):
    df = pd.read_csv(ne_file, sep="\t")
    df.columns = ["GEN", "NE", "L95", "U95"]
    return df


def read_true_ne(ne_file):
    df = pd.read_csv(ne_file, sep="\t")
    return df


def get_ibd_obj_dist_path(resdir, label, gsid):
    # "resdir/10000_sp_neu_rep0/ibddist_ibd/10000_2.0_10.0_none.ibddist.ibdobj.gz"
    p = f"{resdir}/{gsid}_{label}/ibddist_ibd/{gsid}_2.0_10.0_none.ibddist.ibdobj.gz"
    assert Path(p).exists(), f"{p} does not exist"
    return p


def get_ibd_obj_dist_path_repo(resdir, label, gsid, repno):
    # "resdir/10000_sp_neu_rep0/ibddist_ibd/10000_2.0_10.0_none.ibddist.ibdobj.gz"
    gsid = gsid + 10 * repno
    p = f"{resdir}/{gsid}_{label}_rep{repno}/ibddist_ibd/{gsid}_2.0_10.0_none.ibddist.ibdobj.gz"
    assert Path(p).exists(), f"{p} does not exist"
    return p


def get_daf_path(resdir, label, gsid):
    # resdir/10000_sp_neu_rep0/daf/10000.daf
    p = f"{resdir}/{gsid}_{label}/daf/{gsid}.daf"
    assert Path(p).exists(), f"{p} does not exist"
    return p


def get_daf_path_rep(resdir, label, gsid, repno):
    # resdir/10000_sp_neu_rep0/daf/10000.daf
    gsid = gsid + 10 * repno
    p = f"{resdir}/{gsid}_{label}_rep{repno}/daf/{gsid}.daf"
    assert Path(p).exists(), f"{p} does not exist"
    return p


def read_daf(daf_file):
    df = pd.read_csv(daf_file, sep="\t")
    return df


def get_restart_count(resdir, label, gsid, repno):
    # resdir/10000_sp_neu_rep0/restart_count/10000.restart_count
    p = f"{resdir}/{gsid}_{label}/restart_count/{gsid}.restart_count"
    assert Path(p).exists(), f"{p} does not exist"
    return int(Path(p).read_text())


def get_restart_count_rep(resdir, label, gsid, repno):
    # resdir/10000_sp_neu_rep0/restart_count/10000.restart_count
    gsid = gsid + 10 * repno
    p = f"{resdir}/{gsid}_{label}_rep{repno}/restart_count/{gsid}.restart_count"
    assert Path(p).exists(), f"{p} does not exist"
    return int(Path(p).read_text())


def get_axes_position(label):
    sim = SIMULATIONS[label]
    if sim["sim_relatedness"] == 0:
        row = 0
    else:
        i = DELTA_LIST.index(sim["sim_relatedness_delta"])
        ii = RELG_LIST.index(sim["sim_relatedness_g"])
        iii = POWER_LIST.index(sim["sim_relatedness_power"])
        row = i * N_RELG * N_POWER + ii * N_POWER + iii + 1

    j = SELG_LIST.index(sim["g_sel_start"])
    jj = S_LIST.index(sim["s"])
    col = j * N_S + jj
    return (row, col)


def get_axes_position2(label, rmmode):
    sim = SIMULATIONS[label]
    if sim["sim_relatedness"] == 0:
        row = 0
    else:
        i = DELTA_LIST.index(sim["sim_relatedness_delta"])
        ii = RELG_LIST.index(sim["sim_relatedness_g"])
        iii = POWER_LIST.index(sim["sim_relatedness_power"])
        row = i * N_RELG * N_POWER + ii * N_POWER + iii + 1

    j = SELG_LIST.index(sim["g_sel_start"])
    jj = S_LIST.index(sim["s"])
    jjj = ["orig", "rmpeaks"].index(rmmode)
    col = j * N_S * 2 + jj * 2 + jjj
    return (row, col)


def get_title(icol):
    i_g = icol // N_S
    i_s = icol % N_S
    return f"g_sel_start: {SELG_LIST[i_g]}\ns={S_LIST[i_s]}"


def get_title2(icol):
    x = icol
    i = x // (N_S * 2)
    x -= i * (N_S * 2)
    ii = x // 2
    x -= ii * 2
    iii = x
    rmmode = ["orig", "rmpeaks"][iii]
    return f"g_sel_start: {SELG_LIST[i]}\ns={S_LIST[ii]}\n{rmmode}"


def get_ylabel(irow):
    if irow == 0:
        return "Noninbreeding"
    else:
        x = irow - 1
        j = x // (N_RELG * N_POWER)
        x -= j * (N_RELG * N_POWER)
        jj = x // N_POWER
        x -= jj * N_POWER
        jjj = x
        lines = [
            f"Delta={DELTA_LIST[j]}",
            f"RelG={RELG_LIST[jj]}",
            f"Power={POWER_LIST[jjj]}",
        ]
        return "\n".join(lines)


def plot_ne(
    ax,
    true_ne_path,
    neutral_ne_path,
    orig_ne_path,
    rmpeaks_ne_path,
):
    df_true = read_true_ne(true_ne_path)
    df_neutral = read_ne(neutral_ne_path)
    df_orig = read_ne(orig_ne_path)
    df_rmpeaks = read_ne(rmpeaks_ne_path)
    ax.plot(df_neutral.GEN, df_neutral.NE / 4, "b", label="neutral")
    ax.plot(df_orig.GEN, df_orig.NE / 4, "r", label="Orig")
    ax.plot(df_rmpeaks.GEN, df_rmpeaks.NE / 4, "r--", label="Rmpeaks")
    ax.plot(df_true.GEN, df_true.NE, "k--", label="True")
    ax.set_xlim(0, 100)
    ax.set_ylim(100, 1e6)
    ax.set_yscale("log")
    ax.legend()


def plot_daf(ax, daf):
    nchrom = daf.shape[1] - 1
    for i in range(nchrom):
        ichr = i + 1
        colname = f"DAF_CHR{ichr}"
        ax.plot(daf.GEN, daf[colname])


ncols = N_S * N_SELG
nrows = 1 + N_DELTA * N_RELG * N_POWER
fig, axes = plt.subplots(
    nrows=nrows,
    ncols=ncols,
    figsize=(ncols * 7, nrows * 2),
    sharex=True,
    sharey=True,
    constrained_layout=True,
)
for label in SIMULATIONS.keys():
    (i, j) = get_axes_position(label)
    ax = axes[i, j]
    gsid = SIMULATIONS[label]["genome_set_id"]
    p = get_ibd_obj_dist_path(RESDIR, label, gsid)
    ibd = IBD.pickle_load(p)
    print("cov_df", ibd._cov_df.head(1))
    x = ibd._cov_df.GwStart
    y = ibd._cov_df.Coverage
    ax.plot(x, y)
    # ax.set_ylim(-100, 30000)
    if j == 0:
        ax.set_ylabel(get_ylabel(i))
    if i == 0:
        ax.set_title(get_title(j))
fig.savefig("coverage_sp.png")


ncols = N_S * N_SELG
nrows = 1 + N_DELTA * N_RELG * N_POWER
fig, axes = plt.subplots(
    nrows=nrows,
    ncols=ncols,
    figsize=(ncols * 2, nrows * 2),
    sharex=True,
    sharey=True,
    constrained_layout=True,
)
for label in SIMULATIONS.keys():
    (i, j) = get_axes_position(label)
    ax = axes[i, j]
    gsid = SIMULATIONS[label]["genome_set_id"]
    if j == 0:
        ax.set_ylabel(get_ylabel(i))
    if i == 0:
        ax.set_title(get_title(j))
    df_true = read_true_ne(get_true_ne_path(RESDIR, label, gsid))
    df_orig = read_ne(get_ne_path(RESDIR, label, gsid, "orig"))
    df_rmpeaks = read_ne(get_ne_path(RESDIR, label, gsid, "rmpeaks"))
    ax.plot(df_orig.GEN, df_orig.NE / 4, "r", label="Orig")
    ax.plot(df_rmpeaks.GEN, df_rmpeaks.NE / 4, "r--", label="Rmpeaks")
    ax.plot(df_true.GEN, df_true.NE, "k--", label="True")
    ax.set_xlim(0, 100)
    ax.set_ylim(100, 1e6)
    ax.set_yscale("log")
    ax.legend()
fig.savefig("ne_sp.png")
