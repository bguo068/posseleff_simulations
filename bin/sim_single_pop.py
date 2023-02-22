#! /usr/bin/env python3

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import subprocess
import tskit
import msprime
import pyslim
import io
import warnings
import argparse
from pathlib import Path
import sys
import gzip
import json

slim_script_dir = Path(__file__).parents[1] / "slim"


class SlimMsprimeSimulatorForSinglePop:
    def __init__(
        self,
        chrno=1,
        seqlen=100 * 15_000,
        selpos=None,
        num_origins=1,
        N=10000,
        h=0.5,
        s=0.2,
        g_sel_start=80,
        r=0.01 / 15_000,
        sim_relatedness=0,
        #
        ne_change_g=200,
        N0=3000,
        u=1e-8,
        nsam=100,
        slim_script=str(slim_script_dir / "single_pop.slim"),
    ):
        """
        init simulation parameters

        # pyslim doc: https://tskit.dev/pyslim/docs/stable/introduction.html

        """
        local_dict = {k: v for k, v in locals().items() if k != "self"}
        with open(f"config_{chrno}.json", "w") as fp:
            json.dump(local_dict, fp)
        self.chrno = chrno
        self.seqlen = int(seqlen)
        self.selpos = int(seqlen // 2 if selpos is None else selpos)
        self.num_origins = num_origins
        self.N = N
        self.h = h
        self.s = s
        self.g_sel_start = g_sel_start
        self.r = r
        self.sim_relatedness = sim_relatedness

        self.ne_change_g = ne_change_g
        self.N0 = N0
        self.u = u
        self.nsam = nsam
        self.slim_script = slim_script
        assert self.nsam // 2 <= self.N0

    def _run_slim(self, idx, slim_seed) -> str:
        outid = idx

        # run slim
        slim_params = {
            "L": self.seqlen,
            "selpos": self.selpos,
            "num_origins": self.num_origins,
            "N": self.N,
            "h": self.h,
            "s": self.s,
            "g_sel_start": self.g_sel_start,
            "r": self.r,
            "outid": outid,
            "max_restart": 100,
            "sim_relatedness": self.sim_relatedness,
            "N0": self.N0,
            "g_ne_change_start": self.ne_change_g,
        }
        slim_params_str = " ".join([f"-d {k}={v}" for k, v in slim_params.items()])
        seed_str = f"-seed {slim_seed}" if slim_seed is not None else ""
        cmd = f"slim {slim_params_str} {seed_str} {self.slim_script}"
        print(f"simulate chrom with id {outid}")
        print(cmd)
        res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if res.returncode != 0:
            print(cmd)
            sys.stderr.write(res.stderr)
            sys.exit(1)

        self.slim_stdout = res.stdout
        self.slim_tree_fn = f"tmp_slim_out_single_pop_{outid}.trees"

    def _parse_slim_stdout(self):
        """parse stdout to get true ne and restart_count"""
        ne_lines = ["GEN\tNE"]
        daf_lines = ["GEN\tDAF"]
        restart_count = 0
        for line in self.slim_stdout:
            if line.startswith("restart_count"):
                restart_count = int(line.split("\t")[1])
            elif line.startswith("True_Ne"):
                ne_lines.append(line.replace("True_Ne\t", ""))
            elif line.startswith("DAF"):
                daf_lines.append(line.replace("DAF\t", ""))

        self.true_ne_df = pd.read_csv(io.StringIO("\n".join(ne_lines)), sep="\t")
        self.daf_df = pd.read_csv(io.StringIO("\n".join(daf_lines)), sep="\t")

        # Given slim could have restarted the selection process, dedup is needed
        self.daf_df.drop_duplicates("generation", keep="last", inplace=True)
        self.true_ne_df.drop_duplicates("generation", keep="last", inplace=True)

        self.restart_count = restart_count

    def simulate_a_chromosome(
        self,
        idx,
        rm_slim_tree_files=False,
        slim_seed=None,
        recapitate_seed=None,
    ):

        self._run_slim(idx, slim_seed=slim_seed)
        self._parse_slim_stdout()

        # load slim trees
        ts = tskit.load(self.slim_tree_fn)
        if rm_slim_tree_files:
            Path(self.slim_tree_fn).unlink()

        # suppress warning
        warnings.simplefilter("ignore", msprime.TimeUnitsMismatchWarning)

        # simplify slim trees but keep roots
        alive_inds = pyslim.individuals_alive_at(ts, 0)
        keep_inds = np.random.choice(alive_inds, self.nsam // 2, replace=False)
        keep_nodes = []
        for i in keep_inds:
            keep_nodes.extend(ts.individual(i).nodes)

        sts = ts.simplify(keep_nodes, keep_input_roots=True)

        # recpaitate (finish coalescent)
        rts: tskit.TreeSequence = pyslim.recapitate(
            sts,
            ancestral_Ne=self.N,
            recombination_rate=self.r,
            random_seed=recapitate_seed,
        )

        # simplify for the second time
        cur_sample_nodes = np.nonzero(
            (rts.tables.nodes.flags == tskit.NODE_IS_SAMPLE)
            & (rts.tables.nodes.time == 0)
        )[0]
        sts2 = rts.simplify(cur_sample_nodes)

        # remove exising slim mutation (it has different ref/alt coding)
        site_list = list(range(sts2.num_sites))
        ts = sts2.delete_sites(site_list)

        # mutate
        mts = msprime.sim_mutations(
            ts,
            rate=self.u,
            model=msprime.SLiMMutationModel(type=0),
            keep=True,
        )

        self.mutated_trees = mts
        self.simplified_trees = sts2


def write_peudo_homozygous_vcf(ts_mutated, chrno, out_vcf):
    # extract information from ts
    gt_list = []
    pos_list = []
    ref_list = []
    alt_list = []
    for v in ts_mutated.variants():
        if len(v.alleles) != 2:
            continue
        gt_list.append(v.genotypes)
        ref_list.append(v.alleles[0])
        alt_list.append(v.alleles[1])
        pos_list.append(int(v.position))

    # prep header
    header = f"""##fileformat=VCFv4.2
##source=tskit 0.4.0
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID={chrno},length={int(ts_mutated.sequence_length)}>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
"""
    # colnames = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]

    # variant info
    df1 = pd.DataFrame(
        {
            "#CHROM": chrno,
            "POS": pos_list,
            "ID": ".",
            "REF": ref_list,
            "ALT": alt_list,
            "QUAL": ".",
            "FILTER": "PASS",
            "INFO": ".",
            "FORMAT": "GT",
        }
    )

    # gt info
    df2 = pd.DataFrame(gt_list)
    df2.columns = [f"tsk_{n}" for n in df2.columns]
    df2 = df2.astype(str)
    df2 = df2 + "|" + df2

    # combine variant info with gt infor
    df = pd.concat([df1.reset_index(drop=True), df2.reset_index(drop=True)], axis=1)

    # write
    with gzip.open(out_vcf, "wt") as f:
        f.write(header)
        df.to_csv(f, sep="\t", header=True, index=False)


def prepare_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--chrno", type=int, default=1, help="chromosome number")
    parser.add_argument("--s", type=float, default=0.3, help="selection coefficient")
    parser.add_argument("--h", type=float, default=0.5, help="dominance coefficient")
    parser.add_argument("--selpos_0_1", type=float, default=0.33, help="selpos: 0~1")
    parser.add_argument("--num_origins", type=int, default=1, help="dac under sel")
    parser.add_argument("--g_sel_start", type=int, default=80, help="selStartG")
    parser.add_argument("--u", type=float, default=1e-8, help="only sim anc, u ignored")
    parser.add_argument(
        "--bp_per_cm", type=int, default=15000, help="r = 0.01/bp_per_cm"
    )
    parser.add_argument("--seqlen_in_cm", type=int, default=100, help="chrlen in cM")
    parser.add_argument(
        "--ne_change_start_g", type=int, default=200, help="Ne change time in gen ago)"
    )
    parser.add_argument("--N", type=int, default=10000, help="ancient Ne")
    parser.add_argument("--N0", type=int, default=10000, help="Ne at present")
    parser.add_argument("--nsam", type=int, default=1000, help="no. hap sampled")
    parser.add_argument(
        "--test",
        action="store_true",
        dest="test",
        help="when this is set, nsam values is replaced with 50 and seqlen with 100cM for quick testing",
    )
    parser.add_argument(
        "--sim_relatedness",
        type=int,
        choices=[0, 1],
        default=0,
        help="0: normal simulation; 1: simulate highly related",
    )
    parser.add_argument("--genome_set_id", type=int, required=True)

    args = parser.parse_args()
    print(args)

    # parameters -- simulation
    args.r = 0.01 / args.bp_per_cm
    args.seqlen = args.seqlen_in_cm * args.bp_per_cm

    # parameters -- find true ibd
    # args.min_cm = 2.0
    # args.sample_window = int(0.01 * args.bp_per_cm)
    # args.min_bp = args.min_cm * args.bp_per_cm
    # args.remove_hbd = True
    # args.min_tmrca = 1.5
    # args.out_ibd = f"{args.chrno}.ibd"

    # if test using a small number for nsam and seqlen
    if args.test:
        args.nsam = 50
        args.seqlen = 20 * args.bp_per_cm

    args.selpos_bp = args.seqlen * args.selpos_0_1

    return args


if __name__ == "__main__":
    args = prepare_args()

    # setting parameters for simulator wrapper
    simulator = SlimMsprimeSimulatorForSinglePop(
        chrno=args.chrno,
        seqlen=args.seqlen,
        selpos=args.selpos_bp,
        num_origins=args.num_origins,
        N=args.N,
        h=args.h,
        s=args.s,
        g_sel_start=args.g_sel_start,
        r=args.r,
        sim_relatedness=args.sim_relatedness,
        #
        ne_change_g=args.ne_change_start_g,
        N0=args.N0,
        u=args.u,
        nsam=args.nsam,
    )

    # first do slim simulation and tree sequence recording, then using
    # msprime/pyslim to recapitate (finish coalescence)
    simulator.simulate_a_chromosome(
        idx=args.chrno,
        slim_seed=args.chrno + args.genome_set_id * 14,
        recapitate_seed=args.chrno + args.genome_set_id * 14,
    )

    # output files
    prefix = f"{args.genome_set_id}_{args.chrno}"
    ofn_slim_restart_count = f"{prefix}.restart_count"
    ofn_true_ne = f"{prefix}.true_ne"
    ofn_daf = f"{prefix}.daf"
    ofn_tree = f"{prefix}.trees"
    ofn_vcf = f"{prefix}.vcf.gz"

    # write files
    Path(ofn_slim_restart_count).write_text(f"{simulator.restart_count}")
    simulator.true_ne_df.to_csv(ofn_true_ne, sep="\t", index=None)
    simulator.daf_df.to_csv(ofn_daf, index=None, sep="\t")
    simulator.simplified_trees.dump(ofn_tree)  # tree before neutral mutations are added
    write_peudo_homozygous_vcf(simulator.mutated_trees, args.chrno, ofn_vcf)

    print(
        f"""
    output files:
        {ofn_slim_restart_count}
        {ofn_true_ne}
        {ofn_daf}
        {ofn_tree}
        {ofn_vcf}
    """
    )
