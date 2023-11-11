import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

RESDIR: str = "/local/scratch/bing/posseleff_simulations/simulations/r231109/resdir"
SIMULATIONS = {
    "sp_neu": 10000,
    "sp_s03": 10003,
    "sp_rels00_24c": 240005,
    "sp_rels03_24c": 240004,
    "sp_rels00_23c": 230005,
    "sp_rels03_23c": 230004,
    "sp_s03_g80": 10004,
    "sp_s03_g120": 10005,
}


def get_ne_path(resdir, label, gsid, repno, rmmode):
    gsid = gsid + 10 * repno
    p = f"{resdir}/{gsid}_{label}_rep{repno}/ne_output/{gsid}_2.0_10.0_none_{rmmode}.ne"
    assert Path(p).exists(), f"{p} does not exist"
    return p


def get_true_ne_path(resdir, label, gsid, repno):
    gsid = gsid + 10 * repno
    p = f"{resdir}/{gsid}_{label}_rep{repno}/true_ne/{gsid}.true_ne"
    assert Path(p).exists(), f"{p} does not exist"
    return p


def read_ne(ne_file):
    df = pd.read_csv(ne_file, sep="\t")
    df.columns = ["GEN", "NE", "L95", "U95"]
    return df


def read_true_ne(ne_file):
    df = pd.read_csv(ne_file, sep="\t")
    return df


def get_ibd_obj_dist_path(resdir, label, gsid, repno):
    # "resdir/10000_sp_neu_rep0/ibddist_ibd/10000_2.0_10.0_none.ibddist.ibdobj.gz"
    gsid = gsid + 10 * repno
    p = f"{resdir}/{gsid}_{label}_rep{repno}/ibddist_ibd/{gsid}_2.0_10.0_none.ibddist.ibdobj.gz"
    assert Path(p).exists(), f"{p} does not exist"
    return p


def get_daf_path(resdir, label, gsid, repno):
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
    gsid = gsid + 10 * repno
    p = f"{resdir}/{gsid}_{label}_rep{repno}/restart_count/{gsid}.restart_count"
    assert Path(p).exists(), f"{p} does not exist"
    return int(Path(p).read_text())


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
    # ax.plot(df_true.GEN, df_true.NE, "k--", label="True")
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


# def test1():
#     print(get_ne_path(RESDIR, "sp_neu", 10000, 0, "orig"))
#     print(get_true_ne_path(RESDIR, "sp_neu", 10000, 0))


def plot_ne_group(neu_idx, sel_idx, out_prefix, comments):
    """plot a 5 x 4 block of ne plots, using specified set of neutral and
    selection simulations"""
    fig, axes = plt.subplots(
        nrows=4, ncols=5, sharex=True, sharey=True, figsize=(20, 16)
    )
    for i in range(20):
        # 240034_sp_rels03_24c_rep3/
        row = i // 5
        col = i % 5
        repno = i
        neu_label = list(SIMULATIONS.keys())[neu_idx]
        neu_gsid = SIMULATIONS[neu_label]
        sel_label = list(SIMULATIONS.keys())[sel_idx]
        sel_gsid = SIMULATIONS[sel_label]
        try:
            neutral_ne_path = get_ne_path(RESDIR, neu_label, neu_gsid, repno, "orig")
            true_ne_path = get_true_ne_path(RESDIR, sel_label, sel_gsid, repno)
            orig_ne_path = get_ne_path(RESDIR, sel_label, sel_gsid, repno, "orig")
            rmpeaks_ne_path = get_ne_path(RESDIR, sel_label, sel_gsid, repno, "rmpeaks")
            plot_ne(
                axes[row, col],
                true_ne_path,
                neutral_ne_path,
                orig_ne_path,
                rmpeaks_ne_path,
            )
        except Exception as e:
            print(e)
            raise e
    Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{out_prefix}_neplot.png")
    Path(f"{out_prefix}_comments.txt").write_text(comments)


def plot_daf_group(sel_idx, out_prefix, comments):
    fig, axes = plt.subplots(
        nrows=4, ncols=5, sharex=True, sharey=True, figsize=(20, 16)
    )
    label = list(SIMULATIONS.keys())[sel_idx]
    gsid0 = SIMULATIONS[label]
    for i in range(20):
        row = i // 5
        col = i % 5
        repno = i
        gsid = gsid0
        try:
            daf_path = get_daf_path(RESDIR, label, gsid, repno)
            daf = read_daf(daf_path)
            plot_daf(axes[row, col], daf)
        except Exception as e:
            print(e)
    Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{out_prefix}_daf.png")
    Path(f"{out_prefix}_comments.txt").write_text(comments)


plot_ne_group(
    0,
    1,
    "pngs/noninbreeding/selg50",
    "selection sweeps for 50 generation might be too short to see selection effect on Ne",
)

plot_ne_group(
    0,
    6,
    "pngs/noninbreeding/selg80",
    "selection sweeps for 80 generation is best to see selection effect on Ne",
)
plot_ne_group(
    0,
    7,
    "pngs/noninbreeding/selg120",
    "selection sweeps for 120 generation is too long to see selection effect on Ne",
)

plot_ne_group(
    2,
    3,
    "pngs/inbreeding/g25power1delta0.001",
    "0.001 causes very large change in Ne pattern",
)

plot_ne_group(
    4,
    5,
    "pngs/inbreeding/g25power1delta0.01",
    "delta 0.01 causes moderate change in Ne pattern",
)

plot_daf_group(
    1, "pngs/noninbreeding/selg50", "afreq is still increasing; not close to plateau"
)

plot_daf_group(6, "pngs/noninbreeding/selg80", "close to plateau")
plot_daf_group(7, "pngs/noninbreeding/selg120", "over-plateaued")

plot_daf_group(
    3,
    "pngs/inbreeding/g25power1delta0.001",
    "delta 0.001, there seems to be a drop when inbreeding starts",
)
plot_daf_group(5, "pngs/inbreeding/g24power1delta0.01", "delta0.01, more smooth")
