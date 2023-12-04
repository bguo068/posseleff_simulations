import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from ibdutils.utils.ibdutils import IBD
import json

# RESDIR = "/local/scratch/bing/posseleff_simulations/simulations/r231113/resdir"
RESDIR = "/local/scratch/bing/posseleff_simulations/simulations/r231125/resdir_mp"
SIMULATIONS_ALL = json.loads(Path("./mp_genome_sets.json").read_text())


SELMIG_LIST = [0.01, 0.015, 0.03]

S_LIST = [0.0, 0.3]
SELG_LIST = [80]
SFR_LIST = [0.0, 0.1, 0.2, 0.5]
SFG_LIST = [10, 25, 50, 100, 200]
N_SFR = len(SFR_LIST)
N_SFG = len(SFG_LIST)
N_S = len(S_LIST)
N_SELG = len(SELG_LIST)
# REPNO = [0, 1, 2][0]
REPNO = None


def read_ibdobj_ifm(resdir, label, gsid, rmmode, repno=None):
    # eg. resdir/70000_mp_rels03a/ifm_input/70000_orig.ifm.ibdobj.gz
    if repno is None:
        p = f"{resdir}/{gsid}_{label}/ifm_input/{gsid}_{rmmode}.ifm.ibdobj.gz"
    else:
        gsid = gsid + repno * 10
        p = f"{resdir}/{gsid}_{label}_rep{repno}/ifm_input/{gsid}_{rmmode}.ifm.ibdobj.gz"
    assert Path(p).exists(), f"{p} does not exist"
    return IBD.pickle_load(p)


def get_cluster_ibdmat_and_avg_totalibd(ibd: IBD):
    mat = ibd.make_ibd_matrix()
    n = mat.shape[0]
    sample_reordered, _lnk_matrix = ibd.clust_ibd_matrix(mat)
    reorder_mat = mat.loc[sample_reordered, sample_reordered]
    avg_totalibd = mat.sum().sum() / (n * n - n)
    return reorder_mat, avg_totalibd


def get_ifm_member(resdir, label, gsid, rmmode, repno):
    # "resdir/70001_mp_rels00a/ifm_output/"
    if repno is None:
        p = f"{resdir}/{gsid}_{label}/ifm_output/{gsid}_{rmmode}_member.pq"
    else:
        gsid = gsid + 10 * repno
        p = f"{resdir}/{gsid}_{label}_rep{repno}/ifm_output/{gsid}_{rmmode}_member.pq"

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
    j_s = S_LIST.index(sim["s"])
    j_sfg = SFG_LIST.index(sim["sim_selfing_g"]) if "sim_selfing_g" in sim else 0

    i_ssg = SELG_LIST.index(sim["g_sel_start"])
    i_sfr = (
        SFR_LIST.index(sim["sim_selfing_rate"])
        if "sim_selfing_rate" in sim.keys()
        else 0
    )

    row = 0 if i_sfr == 0 else i_ssg * (N_SFR) + i_sfr
    col = j_sfg * N_S + j_s

    if row == 0:
        # expand current plot to multiple axes
        return [(row, icol) for icol in range(N_S * N_SFG) if icol % (N_S) == j_s]
    else:
        return [(row, col)]


def get_axes_position2(label, rmmode):
    sim = SIMULATIONS[label]
    j_rm = ["orig", "rmpeaks"].index(rmmode)
    j_s = S_LIST.index(sim["s"])
    j_sfg = SFG_LIST.index(sim["sim_selfing_g"]) if "sim_selfing_g" in sim else 0

    i_ssg = SELG_LIST.index(sim["g_sel_start"])
    i_sfr = (
        SFR_LIST.index(sim["sim_selfing_rate"])
        if "sim_selfing_rate" in sim.keys()
        else 0
    )

    row = 0 if i_sfr == 0 else i_ssg * (N_SFR) + i_sfr
    col = j_sfg * (N_S * 2) + j_s * 2 + j_rm

    if row == 0:
        # expand current plot to multiple axes
        return [
            (row, icol)
            for icol in range(N_S * N_SFG * 2)
            if icol % (N_S * 2) == j_rm + j_s * 2
        ]
    else:
        return [(row, col)]


def get_title(icol):
    i_g = icol // N_S
    i_s = icol % N_S
    return f"selfing g: {SFG_LIST[i_g]}\ns={S_LIST[i_s]}"


def get_title2(icol):
    x = icol
    i = x // (N_S * 2)
    x -= i * (N_S * 2)
    ii = x // 2
    x -= ii * 2
    iii = x
    rmmode = ["orig", "rmpeaks"][iii]
    return f"selfing g: {SFG_LIST[i]}\ns={S_LIST[ii]}\n{rmmode}"


def get_ylabel(irow):
    if irow == 0:
        return "NonSelfing"
    else:
        i_gselg = irow // N_SFR
        i_sfr = irow % N_SFG

        lines = [
            f"Sel Gen={ SELG_LIST[i_gselg] }",
            f"Selfing Rate={SFR_LIST[i_sfr]}",
        ]
        return "\n".join(lines)


def get_avg_ibd_prop(ibd: IBD):
    cov_df = ibd._cov_df.copy()
    cov_df["InPeak"] = False

    if (ibd._peaks_df is not None) and (ibd._peaks_df.shape[0] > 0):
        for gwstart, gwend in ibd._peaks_df[["GwStart", "GwEnd"]].itertuples(
            index=None
        ):
            cov_df["InPeak"] |= (cov_df.GwStart < gwend) & (cov_df.GwStart >= gwstart)

    y = ibd._cov_df.Coverage
    n = ibd.get_samples_shared_ibd().shape[0]

    m = y[~cov_df.InPeak].mean()
    avg_prop_wo = m / (n * (n - 1) / 2)
    m = y.mean()
    avg_prop_wi = m / (n * (n - 1) / 2)
    return avg_prop_wi, avg_prop_wo


for selmig in SELMIG_LIST:
    SIMULATIONS = {k: v for k, v in SIMULATIONS_ALL.items() if v["sel_mig"] == selmig}
    group_label = f"selmig{selmig}"

    for i_selmig, selmig in enumerate(SELG_LIST):
        ncols = N_S * N_SFG
        nrows = N_SFR * N_SELG
        fig, axes = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            figsize=(ncols * 7, nrows * 2),
            sharex=True,
            sharey=True,
            constrained_layout=True,
        )
        for label in SIMULATIONS.keys():
            gsid = SIMULATIONS[label]["genome_set_id"]
            for i, j in get_axes_position(label):
                ax = axes[i, j]
                try:
                    ibd = read_ibdobj_ifm(RESDIR, label, gsid, "orig", repno=REPNO)
                    avg_prop_wi, avg_prop_wo = get_avg_ibd_prop(ibd)
                    x = ibd._cov_df.GwStart
                    y = ibd._cov_df.Coverage
                    ax.plot(x, y)
                    ax.text(
                        0.99,
                        0.99,
                        f"ex peaks:{avg_prop_wo:.4f}\nwith peaks: {avg_prop_wi:.4f}",
                        transform=ax.transAxes,
                        va="top",
                        ha="right",
                    )
                except Exception as e:
                    print(e)
                ax.set_ylim(-100, 10000)
                if j == 0:
                    ax.set_ylabel(get_ylabel(i))
                if i == 0:
                    ax.set_title(get_title(j))
        fig.savefig(f"tmp_{group_label}_coverage_mp.png")

        ncols = N_S * N_SFG * 2
        nrows = N_SFR * N_SELG
        fig, axes = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            figsize=(ncols * 2, nrows * 2),
            sharex=True,
            sharey=True,
            constrained_layout=True,
        )
        for label in SIMULATIONS.keys():
            gsid = SIMULATIONS[label]["genome_set_id"]
            for rmmode in ["orig", "rmpeaks"]:
                for i, j in get_axes_position2(label, rmmode):
                    ax = axes[i, j]
                    try:
                        df = get_ifm_member(RESDIR, label, gsid, rmmode, repno=REPNO)
                        df2 = conv_ifm_table_to_matrix(df)
                        ax.imshow(df2, aspect="auto", cmap="Blues", vmin=0, vmax=100)
                    except Exception as e:
                        print(e)
                    if j == 0:
                        ax.set_ylabel(get_ylabel(i))
                    if i == 0:
                        ax.set_title(get_title2(j))
        fig.savefig(f"tmp_{group_label}_membership_mp.png")
