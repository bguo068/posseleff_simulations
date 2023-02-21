#! /usr/bin/env python3

import ibdutils.utils.ibdutils as ibdutils
import ibdutils.runner.ibdne as ibdne
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from pathlib import Path

# parse arguments
pa = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
pa.add_argument("--ibd_files", type=str, nargs=14, required=True)
pa.add_argument("--out_prefix", type=str, required=True)
pa.add_argument("--ibdne_jar", type=str, default=None)
args = pa.parse_args()

assert args.ibdne_jar is not None


# read ibd
genome_14_100 = ibdutils.Genome.get_genome("simu_14chr_100cm")
ibd = ibdutils.IBD(genome=genome_14_100, label="orig")
ibd.read_ibd(ibd_fn_lst=args.ibd_files)


# output files:
ofs_ibd_pq = f"{args.out_prefix}_ibddist_ibd.pq"
ofs_ne_orig_script = f"{args.out_prefix}_ne_orig.sh"
ofs_ne_rmpeaks_script = f"{args.out_prefix}_ne_rmpeaks.sh"


# store combined IBD for IBD distribution analysis
ibd._df.to_parquet(ofs_ibd_pq)


# remove highly relatedness samples
mat = ibd.make_ibd_matrix()
unrelated_samples = ibd.get_unrelated_samples(ibdmat=mat)
ibd.subset_ibd_by_samples(subset_samples=unrelated_samples)


# ######################################
# prepare input for IBDNe

#  remove ibd with tmrca < 1.5 (required according to IBDNe paper)
ibd.filter_ibd_by_time(min_tmrca=1.5)

# convert to heterzygous diploids
# Note: remove_hbd might not remove a lot segments as
# hdb only involves n/2 pairs of a total of n(n-1)/2 pairs (ratio: 1/(n-1)).
# When n is large, the difference is small
ibd.convert_to_heterozygous_diploids(remove_hbd=True)

# calculate coverage and remove peaks
ibd.calc_ibd_cov()
ibd.find_peaks()

ibd2 = ibd.duplicate("rmpeaks")
ibd2.remove_peaks()

# use nerunner (dry_run) to preare inputs and bash scripts
# --- for ibd before removing ibd peaks ------------------------
nerunner = ibdne.IbdNeRunner(
    ibd=ibd,
    input_folder=str(Path("ne_orig").absolute()),
    output_folder=str(Path("ne_orig").absolute()),
    mincm=2,
    minregion=10,
)
nerunner.run(nthreads=10, mem_gb=20, dry_run=True)
# move and rename the script file
Path(nerunner.ouptut_ne_fn).with_suffix(".sh").rename(ofs_ne_orig_script)
# link ibdne.jar file
if not Path(f"ne_orig/ibdne.jar").exists():
    Path(args.ibdne_jar).absolute().symlink_to(f"ne_orig/ibdne.jar")

# --- for ibd after removing ibd peaks ------------------------
nerunner = ibdne.IbdNeRunner(
    ibd=ibd,
    input_folder=str(Path("ne_rmpeaks").absolute()),
    output_folder=str(Path("ne_rmpeaks").absolute()),
    mincm=2,
    minregion=10,
)
nerunner.run(nthreads=10, mem_gb=20, dry_run=True)
# move and rename the script file
Path(nerunner.ouptut_ne_fn).with_suffix(".sh").rename(ofs_ne_rmpeaks_script)
# link ibdne.jar file
if not Path(f"ne_rmpeaks/ibdne.jar").exists():
    Path(args.ibdne_jar).absolute().symlink_to(f"ne_rmpeaks/ibdne.jar")
