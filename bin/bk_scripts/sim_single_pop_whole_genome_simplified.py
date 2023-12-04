#! /usr/bin/env python3

# NOTE: The difference of this script to the `./sim_single_pop_whole_genome.py`
# is that after the forward simulation of multi-chromosome genomes, the genomes
# are split back to individual chromosomes and are recapitated (run coalescent
# simulation) independently to accelerate Initial characterization shows that
# this script does NOT generate expected Ne for neutral, non-assortative mating
# simultaions while the `./sim_single_pop_whole_genome.py` does.  For now, I
# will stick to this script `./sim_single_pop_whole_genome.py`  instead of
# `./sim_single_pop_whole_genome_simplified.py`
#
# NOTE: It would be interesting to explore why
# `./sim_single_pop_whole_genome_simplified.py` does not generate expected Ne
# pattern: (1) coding bug; (2) or this strategy does not work at all

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import subprocess
import tskit
import msprime
import pyslim
import math
import io
import warnings
import argparse
from pathlib import Path
import sys
import gzip
import json
import logging
from ibdutils.utils.others import CheckableParams

slim_script_dir = Path(__file__).parents[1] / "slim"
logging.basicConfig(level=logging.INFO)


class SimParams(CheckableParams):
    def __init__(self) -> None:
        self.nchrs = 14
        self.seqlen = 100 * 15_000
        self.selpos = int(0.33 * 100 * 15_000)
        self.num_origins = 1
        self.N = 10000
        self.h = 0.5
        self.s = 0.3
        self.g_sel_start = 80
        self.r = 0.01 / 15_000
        self.sim_relatedness = 0
        self.sim_relatedness_power = 1.0
        self.sim_relatedness_delta = 0.01
        self.sim_relatedness_g = 40
        #
        self.g_ne_change_start = 200
        self.N0 = 1000
        self.u = 1e-8
        self.nsam = 1000
        self.slim_script = str(slim_script_dir / "single_pop_whole_genome.slim")

        # save values to  __defaults__ and clear above attributes
        super().__init__()

    def prepare_args(self):
        d = self.__defaults__
        p = argparse.ArgumentParser(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        p.add_argument("--nchrs", type=int, default=d["nchrs"])
        p.add_argument("--seqlen", type=int, default=d["seqlen"])
        p.add_argument("--selpos", type=int, default=d["selpos"])
        p.add_argument("--num_origins", type=int, default=d["num_origins"])
        p.add_argument("--N", type=int, default=d["N"])
        p.add_argument("--s", type=float, default=d["s"])
        p.add_argument("--h", type=float, default=d["h"])
        p.add_argument("--g_sel_start", type=int, default=d["g_sel_start"])
        p.add_argument("--r", type=float, default=d["r"])
        p.add_argument(
            "--sim_relatedness", type=int, choices=[0, 1], default=d["sim_relatedness"]
        )
        p.add_argument(
            "--sim_relatedness_power", type=float, default=d["sim_relatedness_power"]
        )
        p.add_argument(
            "--sim_relatedness_delta", type=float, default=d["sim_relatedness_delta"]
        )
        p.add_argument("--sim_relatedness_g", type=int, default=d["sim_relatedness_g"])
        p.add_argument("--g_ne_change_start", type=int, default=d["g_ne_change_start"])
        p.add_argument("--N0", type=int, default=d["N0"])
        p.add_argument("--u", type=float, default=d["u"])
        p.add_argument("--nsam", type=int, default=d["nsam"])

        # extra arguments
        p.add_argument(
            "--test",
            action="store_true",
            dest="test",
            help="when this is set, nsam values is replaced with 50 and seqlen with 100cM for quick testing",
        )
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


class SlimMsprimeSimulatorForSinglePop:
    def __init__(
        self,
        params: SimParams,
    ):
        """
        init simulation parameters

        # pyslim doc: https://tskit.dev/pyslim/docs/stable/introduction.html

        """
        self.params = params
        assert self.params.nsam // 2 <= self.params.N0

    def _run_slim(self, idx, slim_seed) -> str:
        outid = idx

        # run slim
        slim_params = {
            "nchrs": self.params.nchrs,
            "L": self.params.seqlen,
            "selpos": self.params.selpos,
            "num_origins": self.params.num_origins,
            "N": self.params.N,
            "h": self.params.h,
            "s": self.params.s,
            "g_sel_start": self.params.g_sel_start,
            "r": self.params.r,
            "outid": outid,
            "max_restart": 100,
            "sim_relatedness": self.params.sim_relatedness,
            "sim_relatedness_power": self.params.sim_relatedness_power,
            "sim_relatedness_delta": self.params.sim_relatedness_delta,
            "sim_relatedness_g": self.params.sim_relatedness_g,
            "N0": self.params.N0,
            "g_ne_change_start": self.params.g_ne_change_start,
        }
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
        self.slim_tree_fn = f"tmp_slim_out_single_pop_{outid}.trees"

    def _parse_slim_stdout(self):
        """parse stdout to get true ne and restart_count"""
        ne_lines = ["GEN\tNE"]
        daf_lines = [
            # header line
            "\t".join(
                ["GEN"]
                + [f"DAF_CHR{chrno}" for chrno in range(1, 1 + self.params.nchrs)]
            )
        ]
        restart_count = 0
        for line in self.slim_stdout.strip().split("\n"):
            if line.startswith("restart_count"):
                restart_count = int(line.split("\t")[1])
            elif line.startswith("True_Ne"):
                ne_lines.append(line.replace("True_Ne\t", ""))
            # elif line.startswith("in"):
            #     print(line)
            elif line.startswith("DAF"):
                # TODO: need to deal with DAF of multiple selected sites
                # DAF	1
                # _FRQ_	0.918398	0.857072	0.98368	0.772502	0.962413	0.962908	0.971316	0.997033	0.99456	0.98269	0.871414	0.9545
                # _POS_	250000	1000000	1750000	2500000	3250000	4000000	4750000	5500000	6250000	7000000	7750000	8500000
                # _FIXPOS_
                order_freq = [0.0] * self.params.nchrs
                lst = line.strip().split()
                g = lst[1]
                mut_frq = lst[lst.index("_FRQ_") + 1 : lst.index("_POS_")]
                mut_frq = [float(x) for x in mut_frq]
                mut_pos = lst[lst.index("_POS_") + 1 : lst.index("_FIXPOS_")]
                mut_pos = [int(x) for x in mut_pos]
                fix_pos = lst[lst.index("_FIXPOS_") + 1 :]
                fix_pos = [int(x) for x in fix_pos]

                for frq, pos in zip(mut_frq, mut_pos):
                    idx = int(pos) // int(self.params.seqlen)
                    order_freq[idx] = frq
                for pos in fix_pos:
                    idx = int(pos) // int(self.params.seqlen)
                    order_freq[idx] = 1.0

                frq_join_str = "\t".join([f"{f}" for f in order_freq])

                daf_lines.append(f"{g}\t{frq_join_str}")

        self.true_ne_df = pd.read_csv(io.StringIO("\n".join(ne_lines)), sep="\t")
        self.daf_df = pd.read_csv(io.StringIO("\n".join(daf_lines)), sep="\t")

        # Given slim could have restarted the selection process, dedup is needed
        self.daf_df.drop_duplicates("GEN", keep="last", inplace=True)
        self.true_ne_df.drop_duplicates("GEN", keep="last", inplace=True)

        self.restart_count = restart_count

    def simulate_genome(
        self,
        rm_slim_tree_files=False,
        slim_seed=None,
        recapitate_seed=None,
    ):
        # idx is used to indicate chrno but in this case we just assigned to zero
        # this needs to be of a integer type so it will be compatile with the default value set
        # in the
        idx = 0

        logging.info("run slim")
        self._run_slim(idx, slim_seed=slim_seed)
        logging.info("parse slim output")
        self._parse_slim_stdout()

        # load slim trees
        ts = tskit.load(self.slim_tree_fn)
        if rm_slim_tree_files:
            Path(self.slim_tree_fn).unlink()

        # suppress warning
        warnings.simplefilter("ignore", msprime.TimeUnitsMismatchWarning)

        # simplify slim trees but keep roots
        logging.info("run pyslim msprime")
        alive_inds = pyslim.individuals_alive_at(ts, 0)
        keep_inds = np.random.choice(alive_inds, self.params.nsam // 2, replace=False)
        keep_nodes = []
        for i in keep_inds:
            keep_nodes.extend(ts.individual(i).nodes)

        sts = ts.simplify(keep_nodes, keep_input_roots=True)

        # split genomewide trees to chromosomal trees
        # and then run recapitation.
        # this split is to improve the the speed of the recaptitation
        # assuming that coalescent is easier for smaller chunks
        # than the full genomes
        sts2_res = {}
        mts_res = {}
        for ichr in range(self.params.nchrs):
            chrno = ichr + 1
            L = self.params.seqlen
            intervals = [[ichr * L, (ichr + 1) * L - 1]]
            logging.info(f"\ttrim chromosome {chrno}")
            sts_chr = sts.keep_intervals(intervals).trim()

            # remove slim mutation information to get rid of this
            # following error:
            # _tskit.LibraryError: A mutation's time must be < the parent node
            # of the edge on which it occurs, or be marked as 'unknown'.
            # (TSK_ERR_MUTATION_TIME_OLDER_THAN_PARENT_NODE)
            site_list = list(range(sts_chr.num_sites))
            sts_chr = sts_chr.delete_sites(site_list)

            # recpaitate (finish coalescent)
            logging.info(f"\trecaptitate chromosome {chrno}")
            rts: tskit.TreeSequence = pyslim.recapitate(
                sts_chr,
                ancestral_Ne=self.params.N,
                recombination_rate=self.params.r,
                random_seed=recapitate_seed,
            )

            # simplify for the second time
            logging.info(f"\tsimplify chromosome {chrno}")
            cur_sample_nodes = np.nonzero(
                (rts.tables.nodes.flags == tskit.NODE_IS_SAMPLE)
                & (rts.tables.nodes.time == 0)
            )[0]
            sts2 = rts.simplify(cur_sample_nodes)

            # remove exising slim mutation (it has different ref/alt coding)
            site_list = list(range(sts2.num_sites))
            ts = sts2.delete_sites(site_list)

            logging.info(f"\tadd mutations for chromosome {chrno}")
            # mutate
            mts = msprime.sim_mutations(
                ts,
                rate=self.params.u,
                model=msprime.SLiMMutationModel(type=0),
                keep=True,
            )

            sts2_res[chrno] = sts2
            mts_res[chrno] = mts

        self.mutated_trees = mts_res
        self.simplified_trees = sts2_res


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

    # setting parameters for simulator wrapper
    simulator = SlimMsprimeSimulatorForSinglePop(params)

    # first do slim simulation and tree sequence recording, then using
    # msprime/pyslim to recapitate (finish coalescence)
    logging.info("starting simulation")
    simulator.simulate_genome(
        slim_seed=args.genome_set_id,
        recapitate_seed=args.genome_set_id * 3,
    )

    # output files
    prefix = f"{args.genome_set_id}"
    ofn_slim_restart_count = f"{prefix}.restart_count"
    ofn_true_ne = f"{prefix}.true_ne"
    ofn_daf = f"{prefix}.daf"

    # write files
    Path(ofn_slim_restart_count).write_text(f"{simulator.restart_count}")
    simulator.true_ne_df.to_csv(ofn_true_ne, sep="\t", index=None)
    simulator.daf_df.to_csv(ofn_daf, index=None, sep="\t")

    # split the genome-wide trees to trees per chromosomes

    vcf_names = []
    tree_names = []
    for ichr in range(simulator.params.nchrs):
        chrno = ichr + 1
        prefix = f"{args.genome_set_id}_{chrno}"
        ofn_vcf = f"{prefix}.vcf.gz"
        ofn_tree = f"{prefix}.trees"

        # write vcf
        logging.info(f"write vcf: chr{chrno}")
        mts = simulator.mutated_trees[chrno]
        write_peudo_homozygous_vcf(mts, chrno, ofn_vcf)
        vcf_names.append(ofn_vcf)

        # trees without neutral mutations
        simulator.simplified_trees[chrno].dump(ofn_tree)
        tree_names.append(ofn_tree)
        logging.info(f"run cut trees by interval: chr{chrno}")

    print(
        f"""
    output files:
        {ofn_slim_restart_count}
        {ofn_true_ne}
        {ofn_daf}
        {vcf_names}
        {tree_names}
    """
    )
