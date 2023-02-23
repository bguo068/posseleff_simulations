#! /usr/bin/env python3

import ibdutils.utils.ibdutils as ibdutils
import ibdutils.runner.ibdne as ibdne
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from pathlib import Path

ibd_jar_default = str(Path(__file__).parent / "ibdne.jar")

# parse arguments
pa = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
pa.add_argument("--ibd_files", type=str, nargs=14, required=True)
pa.add_argument("--genome_set_id", type=str, required=True)
pa.add_argument("--ibdne_jar", type=str, default=ibd_jar_default)
args = pa.parse_args()

assert args.ibdne_jar is not None

idx = args.genome_set_id

# read ibd
genome_14_100 = ibdutils.Genome.get_genome("simu_14chr_100cm")
ibd = ibdutils.IBD(genome=genome_14_100, label=f"{idx}_orig")
ibd.read_ibd(ibd_fn_lst=args.ibd_files)


# output files:
ofs_ibd_pq = f"{args.genome_set_id}_ibddist_ibd.pq"


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

ibd2 = ibd.duplicate(f"{idx}_rmpeaks")
ibd2.remove_peaks()

# link ibdne.jar file
if not Path(f"ibdne.jar").exists():
    assert Path(args.ibdne_jar).exists()
    this = Path("ibdne.jar")
    target = Path(args.ibdne_jar).absolute()
    this.symlink_to(target)
    print(f"link {this} -> {target}")

# use nerunner (dry_run) to preare inputs and bash scripts
# --- for ibd before removing ibd peaks ------------------------
nerunner = ibdne.IbdNeRunner(
    ibd=ibd,
    input_folder=".",
    output_folder=".",
    mincm=2,
    minregion=10,
)
nerunner.run(nthreads=10, mem_gb=20, dry_run=True)

# --- for ibd after removing ibd peaks ------------------------
nerunner = ibdne.IbdNeRunner(
    ibd=ibd2,
    input_folder=".",
    output_folder=".",
    mincm=2,
    minregion=10,
)
nerunner.run(nthreads=10, mem_gb=20, dry_run=True)

print(
    f"""
      output files:
        {ofs_ibd_pq}
        ibdne.jar
        {idx}_orig.sh
        {idx}_orig.map
        {idx}_orig.ibd.gz
        {idx}_rmpeaks.sh
        {idx}_rmpeaks.map
        {idx}_rmpeaks.ibd.gz
      """
)
