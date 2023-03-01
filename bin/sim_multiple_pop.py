#! /usr/bin/env python3

import msprime
import argparse
import demesdraw
import matplotlib.pyplot as plt
from pathlib import Path
import subprocess
import pyslim
import numpy as np
import tskit
import pandas as pd
import io
import sys
import logging
from pathlib import Path
import gzip
from ibdutils.utils.others import CheckableParams

slim_script_dir = Path(__file__).parents[1] / "slim"

# NOTE:
# simulation:
# Two reasons for not using msprime geneicsweep model for simulating selection
# 1. There is no way to specifiy the starting time for selection;
# 2. geneic sweeps model only allows for single population and constant size


class SimParams(CheckableParams):
    def __init__(self) -> None:
        self.chrno = 1
        self.seqlen = 100 * 15000
        self.selpos = int(100 * 15000 * 0.33)
        self.num_origins = 1
        self.N = 10000
        self.h = 0.5
        self.s = 0.3
        self.g_sel_start = 80
        self.r = 0.01 / 15_000
        self.sim_relatedness = 0
        #
        self.mig = 1e-5
        self.sel_mig = 0.01
        self.npop = 5
        self.nsam = 200
        self.Tsplit = 500
        self.u = 1e-8
        self.slim_script = str(slim_script_dir / "multiple_pop.slim")

        # save values to  __defaults__ and clear above attributes
        super().__init__()

    def prepare_args(self):
        d = self.__defaults__
        p = argparse.ArgumentParser(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

        p.add_argument("--chrno", type=int, default=d["chrno"])
        p.add_argument("--seqlen", type=int, default=d["seqlen"])
        p.add_argument("--selpos", type=int, default=d["selpos"])
        p.add_argument("--num_origins", type=int, default=d["num_origins"])
        p.add_argument("--h", type=float, default=d["h"])
        p.add_argument("--s", type=float, default=d["s"])
        p.add_argument("--N", type=int, default=d["N"])
        p.add_argument("--g_sel_start", type=int, default=d["g_sel_start"])
        p.add_argument(
            "--sim_relatedness", choices=[0, 1], type=int, default=d["sim_relatedness"]
        )
        p.add_argument("--r", type=float, default=d["r"])

        p.add_argument("--mig", type=float, default=d["mig"])
        p.add_argument("--sel_mig", type=float, default=d["sel_mig"])
        p.add_argument("--npop", type=int, default=d["npop"])
        p.add_argument("--nsam", type=int, default=d["nsam"])
        p.add_argument("--Tsplit", type=int, default=d["Tsplit"])
        p.add_argument("--u", type=float, default=d["u"])

        p.add_argument("--test", action="store_true", default=False)
        p.add_argument("--genome_set_id", type=int, required=True)

        args = p.parse_args()

        # if test using a small number for nsam and seqlen
        if args.test:
            args.nsam = 50
            args.seqlen = 20 * int(0.01 / args.r)
            args.selpos = args.seqlen // 3

        # save to params
        for k, v in vars(args).items():
            if k not in ["test", "genome_set_id"]:
                self.__setattr__(k, v)
        self.fill_defaults()

        return args


class SlimMsprimeSimulatorForMultiplePop:
    def __init__(self, params: SimParams):
        self.params = params

    def _run_slim(self, idx, slim_seed):
        outid = idx

        # params dict
        slim_params = dict(
            L=self.params.seqlen,
            selpos=self.params.selpos,
            sim_relatedness=self.params.sim_relatedness,
            num_origins=self.params.num_origins,
            N=self.params.N,
            h=self.params.h,
            s=self.params.s,
            g_sel_start=self.params.g_sel_start,
            r=self.params.r,
            outid=outid,
            max_restart=100,
            npop=self.params.npop,
            sel_mig=self.params.sel_mig,
        )

        slim_params_str = " ".join([f"-d {k}={v}" for k, v in slim_params.items()])
        seed_str = f"-seed {slim_seed}" if slim_seed is not None else ""
        cmd = f"slim {slim_params_str} {seed_str} {self.params.slim_script}"
        print(f"simulate chrom with id {outid}")
        print(cmd)
        res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if res.returncode != 0:
            print(cmd)
            sys.stderr.write(res.stderr)
            sys.exit(1)

        self.slim_stdout = res.stdout
        self.slim_tree_fn = f"tmp_slim_out_multiple_pop_{outid}.trees"

    def _parse_slim_stdout(self):
        restart_count = 1

        daf_header_str = "\t".join(
            ["Type", "RestartCount", "GEN"]
            + [f"DAF{i}" for i in range(1, self.params.npop + 1)]
        )
        daf_list = [daf_header_str]

        for line in self.slim_stdout.split("\n"):
            if line.startswith("restart count"):
                tokens = line.split("\t")
                restart_count = int(tokens[1])
            elif line.startswith("DAF\t"):
                daf_list.append(line)

        # read the selected lines into pd dataframe and filter the line for the last
        # restart
        daf_df = pd.read_csv(io.StringIO("\n".join(daf_list)), sep="\t")
        print(daf_df)
        daf_df = daf_df[daf_df.RestartCount == restart_count]

        self.daf_df = daf_df
        self.restart_count = restart_count

    def plot_demography(self, ofn_demog):
        g = self.demog.to_demes()
        w = demesdraw.utils.size_max(g) * 2
        positions = dict(Ancestral=2 * w, p1=0, p2=w, p3=2 * w, p4=3 * w, p5=4 * w)
        fig, ax = plt.subplots()

        demesdraw.tubes(
            g, log_time=True, labels="xticks-legend", positions=positions, ax=ax
        )
        fig.savefig(ofn_demog, dpi=600)

    def _run_msprime(self, recapitate_seed=None, rm_slim_tree_files=False):
        # load the slim trees file
        orig_ts = tskit.load(self.slim_tree_fn)
        if rm_slim_tree_files:
            Path(self.slim_tree_fn).unlink()

        # sample nodes from each population at t = 0
        sample_nodes = []
        for pop in range(1, self.params.npop + 1):
            ind_alive = pyslim.individuals_alive_at(orig_ts, 0, population=pop)
            sampled_ind = np.random.choice(
                ind_alive, size=self.params.nsam // 2, replace=False
            )
            for i in sampled_ind:
                sample_nodes.extend(orig_ts.individual(i).nodes)

        # simplify tree with keep_input_roots on for recapitation
        sts = orig_ts.simplify(sample_nodes, keep_input_roots=True)

        # prepare demography for recapitation with migration
        demog = msprime.Demography.from_tree_sequence(sts)
        pop_names = []

        #       set initial size
        for pop in demog.populations:
            if pop.name == "pop_0":
                continue
            pop.initial_size = self.params.N
            pop_names.append(pop.name)

        #       add ancestral population.
        demog.add_population(
            name="Ancestral",
            initial_size=self.params.N,
            # the slim_id is required meta for slim generated tree
            extra_metadata={"slim_id": self.params.npop + 1},
        )

        #       and migration events
        slim_time = orig_ts.metadata["SLiM"]["cycle"]
        for i in range(1, self.params.npop):
            demog.add_symmetric_migration_rate_change(
                time=slim_time, populations=[f"p{i+1}", f"p{i}"], rate=self.params.mig
            )

        #       set population join event
        demog.add_population_split(
            self.params.Tsplit, derived=pop_names, ancestral="Ancestral"
        )
        self.demog = demog

        # recapitulate
        rts = pyslim.recapitate(
            sts,
            demography=demog,  # defined above for migration in the coalescent part
            recombination_rate=self.params.r,
            random_seed=recapitate_seed,
        )

        # second time simplification
        cur_sample_nodes = np.nonzero(
            (rts.tables.nodes.flags == tskit.NODE_IS_SAMPLE)
            & (rts.tables.nodes.time == 0)
        )[0]
        sts2 = rts.simplify(cur_sample_nodes)

        # mutate
        mts = msprime.sim_mutations(
            sts2,
            rate=self.params.u,
            model=msprime.SLiMMutationModel(type=0),
            keep=True,
        )

        self.mutated_trees = mts
        self.simplified_trees = sts2

    def simulate_a_chromosome(
        self,
        idx,
        rm_slim_tree_files=False,
        slim_seed=None,
        recapitate_seed=None,
    ):
        logging.info("running slim")
        self._run_slim(idx, slim_seed=slim_seed)
        self._parse_slim_stdout()
        logging.info("running msprime")
        self._run_msprime(recapitate_seed, rm_slim_tree_files=rm_slim_tree_files)


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


if __name__ == "__main__":
    params = SimParams()
    args = params.prepare_args()
    print(args)

    logging.basicConfig(level=logging.INFO)

    simulator = SlimMsprimeSimulatorForMultiplePop(params)

    simulator.simulate_a_chromosome(
        idx=args.chrno,
        slim_seed=args.chrno + args.genome_set_id * 14,
        recapitate_seed=args.chrno + args.genome_set_id * 14,
    )

    # output files
    prefix = f"{args.genome_set_id}_{args.chrno}"
    ofn_slim_restart_count = f"{prefix}.restart_count"
    ofn_daf = f"{prefix}.daf"
    ofn_tree = f"{prefix}.trees"
    ofn_vcf = f"{prefix}.vcf.gz"
    ofn_demog = f"{prefix}_demog.png"

    # write files
    Path(ofn_slim_restart_count).write_text(f"{simulator.restart_count}")
    simulator.daf_df.to_csv(ofn_daf, index=None, sep="\t")
    simulator.simplified_trees.dump(ofn_tree)  # tree before neutral mutations are added
    write_peudo_homozygous_vcf(simulator.mutated_trees, args.chrno, ofn_vcf)
    simulator.plot_demography(ofn_demog)

    print(
        f"""
    output files:
        {ofn_slim_restart_count}
        {ofn_daf}
        {ofn_tree}
        {ofn_vcf}
        {ofn_demog}
    """
    )
