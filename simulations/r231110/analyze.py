import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

RESDIR: str = "/local/scratch/bing/posseleff_simulations/simulations/r231110/resdir"
SIMULATIONS = dict(
    mp_s00=dict(s=0.0, genome_set_id=20000),
    mp_s03=dict(s=0.3, genome_set_id=20003),
    mp_rels03a=dict(
        sim_relatedness=1,
        s=0.3,
        sim_relatedness_power=1,
        sim_relatedness_delta=0.01,
        sim_relatedness_g=100,
        genome_set_id=70000,
    ),
    mp_rels00a=dict(
        sim_relatedness=1,
        s=0.0,
        sim_relatedness_power=1,
        sim_relatedness_delta=0.01,
        sim_relatedness_g=100,
        genome_set_id=70001,
    ),
    mp_rels03b=dict(
        sim_relatedness=1,
        s=0.3,
        sim_relatedness_power=1,
        sim_relatedness_delta=0.01,
        sim_relatedness_g=50,
        genome_set_id=70002,
    ),
    mp_rels00b=dict(
        sim_relatedness=1,
        s=0.0,
        sim_relatedness_power=1,
        sim_relatedness_delta=0.01,
        sim_relatedness_g=50,
        genome_set_id=70003,
    ),
    mp_rels03c=dict(
        sim_relatedness=1,
        s=0.3,
        sim_relatedness_power=1,
        sim_relatedness_delta=0.01,
        sim_relatedness_g=25,
        genome_set_id=70004,
    ),
    mp_rels00c=dict(
        sim_relatedness=1,
        s=0.0,
        sim_relatedness_power=1,
        sim_relatedness_delta=0.01,
        sim_relatedness_g=25,
        genome_set_id=70005,
    ),
    mp_rels03d=dict(
        sim_relatedness=1,
        s=0.3,
        sim_relatedness_power=1,
        sim_relatedness_delta=0.001,
        sim_relatedness_g=100,
        genome_set_id=80000,
    ),
    mp_rels00d=dict(
        sim_relatedness=1,
        s=0.0,
        sim_relatedness_power=1,
        sim_relatedness_delta=0.001,
        sim_relatedness_g=100,
        genome_set_id=80001,
    ),
    mp_rels03e=dict(
        sim_relatedness=1,
        s=0.3,
        sim_relatedness_power=1,
        sim_relatedness_delta=0.001,
        sim_relatedness_g=50,
        genome_set_id=80002,
    ),
    mp_rels00e=dict(
        sim_relatedness=1,
        s=0.0,
        sim_relatedness_power=1,
        sim_relatedness_delta=0.001,
        sim_relatedness_g=50,
        genome_set_id=80003,
    ),
    mp_rels03f=dict(
        sim_relatedness=1,
        s=0.3,
        sim_relatedness_power=1,
        sim_relatedness_delta=0.001,
        sim_relatedness_g=25,
        genome_set_id=80004,
    ),
    mp_rels00f=dict(
        sim_relatedness=1,
        s=0.0,
        sim_relatedness_power=1,
        sim_relatedness_delta=0.001,
        sim_relatedness_g=25,
        genome_set_id=80005,
    ),
    mp_rels03g=dict(
        sim_relatedness=1,
        s=0.3,
        sim_relatedness_power=1,
        sim_relatedness_delta=1e-4,
        sim_relatedness_g=100,
        genome_set_id=90000,
    ),
    mp_rels00g=dict(
        sim_relatedness=1,
        s=0.0,
        sim_relatedness_power=1,
        sim_relatedness_delta=1e-4,
        sim_relatedness_g=100,
        genome_set_id=90001,
    ),
    mp_rels03h=dict(
        sim_relatedness=1,
        s=0.3,
        sim_relatedness_power=1,
        sim_relatedness_delta=1e-4,
        sim_relatedness_g=50,
        genome_set_id=90002,
    ),
    mp_rels00h=dict(
        sim_relatedness=1,
        s=0.0,
        sim_relatedness_power=1,
        sim_relatedness_delta=1e-4,
        sim_relatedness_g=50,
        genome_set_id=90003,
    ),
    mp_rels03i=dict(
        sim_relatedness=1,
        s=0.3,
        sim_relatedness_power=1,
        sim_relatedness_delta=1e-4,
        sim_relatedness_g=25,
        genome_set_id=90004,
    ),
    mp_rels00i=dict(
        sim_relatedness=1,
        s=0.0,
        sim_relatedness_power=1,
        sim_relatedness_delta=1e-4,
        sim_relatedness_g=25,
        genome_set_id=90005,
    ),
)

DELTA_LIST = [0.01, 1e-3, 1e-4]
RELG_LIST = [100, 50, 25]
S_LIST = [0.0, 0.3]


def get_ifm_member(resdir, label, gsid, rmmode):
    # "resdir/70001_mp_rels00a/ifm_output/"
    p = f"{resdir}/{gsid}_{label}/ifm_output/{gsid}_{rmmode}_member.pq"
    # print(p)
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


labels_list = [
    "mp_s00",
    "mp_s03",
    "mp_rels00a",
    "mp_rels03a",
    "mp_rels00b",
    "mp_rels03b",
    "mp_rels00c",
    "mp_rels03c",
    "mp_rels00d",
    "mp_rels03d",
    "mp_rels00e",
    "mp_rels03e",
    "mp_rels00f",
    "mp_rels03f",
    "mp_rels00g",
    "mp_rels03g",
    "mp_rels00h",
    "mp_rels03h",
    "mp_rels00i",
    "mp_rels03i",
]
nlabels = len(labels_list)

fig, axes = plt.subplots(
    figsize=(20, 8), ncols=nlabels // 2, nrows=4, sharex=True, sharey=True
)
for i, label in enumerate(labels_list):
    col = i // 2
    for irm, rmmode in enumerate(["orig", "rmpeaks"]):
        ylabel = ["NEUTRAL" if i % 2 == 0 else "SELECTION"]
        ylabel = ylabel + [rmmode] + ["True label"]
        ylabel = "\n".join(ylabel)
        row = (i % 2) * 2 + irm
        sim = SIMULATIONS[label]
        title = [
            f"{k.replace('sim_relatedness_', '')}={v}"
            for k, v in sim.items()
            if k not in ["s", "sim_relatedness", "genome_set_id"]
        ]
        title = "\n".join(title)
        df = get_ifm_member(RESDIR, label, sim["genome_set_id"], rmmode)
        df = conv_ifm_table_to_matrix(df)
        ax = axes[row, col]
        ax.imshow(df, aspect="auto", cmap="Blues", vmin=0, vmax=100)
        if col == 0:
            ax.set_ylabel(ylabel)
        if row == 3:
            ax.set_xlabel("Inferred label")
        if row == 0:
            ax.set_title(title)
fig.savefig("ifm_table.png", dpi=600)
