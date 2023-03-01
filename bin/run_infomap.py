#! /usr/bin/env python3

import pandas as pd
import numpy as np
from ibdutils.utils.ibdutils import IBD, Genome
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser


def parse_args():
    p = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument("--ibd_pq", type=str, required=True)
    p.add_argument("--npop", type=int, required=True)
    p.add_argument("--nsam", type=int, required=True)
    p.add_argument("--genome_set_id", type=int, required=True)
    p.add_argument("--cut_mode", type=str, required=True)
    p.add_argument("--ntrails", type=int, default=1000)
    p.add_argument(
        "--transform", type=str, choices=["square", "cube", "none"], default="square"
    )
    args = p.parse_args()
    if args.transform == "none":
        args.transform = None

    return p.parse_args()


def run(args) -> pd.DataFrame:

    g = Genome.get_genome("simu_14chr_100cm")
    ibd = IBD(genome=g, label=f"gsid_{args.genome_set_id}")
    ibd.read_ibd([args.ibd_pq], format="parquet")

    # make meta data
    meta = pd.DataFrame(
        {
            "Sample": np.arange(args.nsam * args.npop),  # use haploid here
            "Population": np.repeat(np.arange(args.npop), args.nsam),
        }
    )

    mat = ibd.make_ibd_matrix()
    member_df = ibd.call_infomap_get_member_df(
        mat, meta, trials=args.ntrails, transform=args.transform
    )

    return member_df


def main():
    args = parse_args()
    member_df = run(args)

    ofs = f"{args.genome_set_id}_{args.cut_mode}_member.pq"
    member_df.to_parquet(ofs)
    print(member_df)


main()
