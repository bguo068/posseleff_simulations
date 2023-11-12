import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from ibdutils.utils.ibdutils import IBD
import json

# RESDIR = "/local/scratch/bing/posseleff_simulations/simulations/r231111/resdir"
RESDIR = "./resdir"
SIMULATIONS = json.loads(Path("./genome_sets.json").read_text())
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


def read_ibdobj_ifm(resdir, label, gsid, rmmode):
    # eg. resdir/70000_mp_rels03a/ifm_input/70000_orig.ifm.ibdobj.gz
    p = f"{resdir}/{gsid}_{label}/ifm_input/{gsid}_{rmmode}.ifm.ibdobj.gz"
    assert Path(p).exists()
    return IBD.pickle_load(p)


def get_cluster_ibdmat_and_avg_totalibd(ibd: IBD):
    mat = ibd.make_ibd_matrix()
    n = mat.shape[0]
    sample_reordered, _lnk_matrix = ibd.clust_ibd_matrix(mat)
    reorder_mat = mat.loc[sample_reordered, sample_reordered]
    avg_totalibd = mat.sum().sum() / (n * n - n)
    return reorder_mat, avg_totalibd


def get_ifm_member(resdir, label, gsid, rmmode):
    # "resdir/70001_mp_rels00a/ifm_output/"
    p = f"{resdir}/{gsid}_{label}/ifm_output/{gsid}_{rmmode}_member.pq"
    assert Path(p).exists()
    return pd.read_parquet(p)


def conv_ifm_table_to_matrix(df, keep_inferred_comms=5):
    nrows = keep_inferred_comms
    df = (
        df.groupby(["Population", "Rank"])["Sample"]
        .count()
        .unstack(0)
        .fillna(0)
        .astype(int)
    )
    # make sure the matrix has enough rows
    if df.shape[0] >= nrows:
        df = df.iloc[:nrows, :]
    else:
        nrows_add = nrows - df.shape[0]
        values = np.zeros((nrows_add, df.shape[1]), dtype=int)
        index = list(range(df.shape[0], df.shape[0] + nrows_add))
        df2 = pd.DataFrame(values, index=index, columns=df.columns)
        df = pd.concat([df, df2], axis=0).rename_axis("Rank", axis=0)
    df = df.sort_values(list(df.columns), ascending=False).transpose()
    return df


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
    ibd = read_ibdobj_ifm(RESDIR, label, gsid, "orig")
    x = ibd._cov_df.GwStart
    y = ibd._cov_df.Coverage
    ax.plot(x, y)
    ax.set_ylim(-100, 3000)
    if j == 0:
        ax.set_ylabel(get_ylabel(i))
    if i == 0:
        ax.set_title(get_title(j))
fig.savefig("coverage.png", dpi=600)

ncols = N_S * N_SELG * 2  # times 2 because we need to plot before/after remove peaks
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
    for rmmode in ["orig", "rmpeaks"]:
        (i, j) = get_axes_position2(label, rmmode)
        ax = axes[i, j]
        gsid = SIMULATIONS[label]["genome_set_id"]
        try:
            df = get_ifm_member(RESDIR, label, gsid, rmmode)
            df2 = conv_ifm_table_to_matrix(df)
            ax.imshow(df2, aspect="auto", cmap="Blues", vmin=0, vmax=100)
        except Exception as e:
            print(e)
        if j == 0:
            ax.set_ylabel(get_ylabel(i))
        if i == 0:
            ax.set_title(get_title2(j))
fig.savefig("membership.png", dpi=600)
