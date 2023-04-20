#! /usr/bin/env python3

import ibdutils.utils.ibdutils as ibdutils
import ibdutils.runner.ibdne as ibdne
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from pathlib import Path

ibd_jar_default = str(Path(__file__).parent / "ibdne.jar")

# parse arguments
pa = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
pa.add_argument("--ibd_files", type=str, nargs=14, required=True)
pa.add_argument("--vcf_files", type=str, nargs=14, required=True)
pa.add_argument("--genome_set_id", type=str, required=True)
pa.add_argument("--ibdne_jar", type=str, default=ibd_jar_default)
pa.add_argument("--ibdne_mincm", type=float, default=2)
pa.add_argument("--ibdne_minregion", type=float, default=10)
pa.add_argument("--ibdne_flatmeth", type=str, default="none")
args = pa.parse_args()

assert args.ibdne_jar is not None
assert args.ibdne_flatmeth in ["none", "merge", "keep_hap_1_only"]

idx = args.genome_set_id
ibdne_mincm = args.ibdne_mincm
ibdne_minregion = args.ibdne_minregion

# output prefix string
label_str = "_".join(
    [
        str(x)
        for x in [
            args.genome_set_id,
            args.ibdne_mincm,
            args.ibdne_minregion,
            args.ibdne_flatmeth,
        ]
    ]
)

# read ibd
genome_14_100 = ibdutils.Genome.get_genome("simu_14chr_100cm")
ibd = ibdutils.IBD(genome=genome_14_100, label=f"{label_str}_orig")
ibd.read_ibd(ibd_fn_lst=args.ibd_files)
ibd.calc_ibd_cov()
ibd.find_peaks()


# output files:
of_ibddist_obj = f"{label_str}.ibddist.ibdobj.gz"


# store combined IBD for IBD distribution analysis
ibd.pickle_dump(of_ibddist_obj)


# remove highly relatedness samples
mat = ibd.make_ibd_matrix()
unrelated_samples = ibd.get_unrelated_samples(ibdmat=mat)
ibd.subset_ibd_by_samples(subset_samples=unrelated_samples)


# ######################################
# prepare input for IBDNe

#  remove ibd with tmrca < 1.5 (required according to IBDNe paper)
ibd.filter_ibd_by_time(min_tmrca=1.5)

ibd.filter_ibd_by_length(min_seg_cm=ibdne_mincm)

# calculate XiR,s
xirs_df = ibd.calc_xirs(vcf_fn_lst=args.vcf_files, min_maf=0.01)

# convert to heterzygous diploids
# Note: remove_hbd might not remove a lot segments as
# hdb only involves n/2 pairs of a total of n(n-1)/2 pairs (ratio: 1/(n-1)).
# When n is large, the difference is small
ibd.convert_to_heterozygous_diploids(remove_hbd=True)

if args.ibdne_flatmeth != "none":
    ibd.flatten_diploid_ibd(method=args.ibdne_flatmeth)

# calculate coverage and remove peaks
ibd.calc_ibd_cov()
ibd.find_peaks()

ibd.filter_peaks_by_xirs(xirs_df)


of_orig_ibdne_obj = f"{label_str}_orig.ibdne.ibdobj.gz"
ibd.pickle_dump(of_orig_ibdne_obj)

ibd2 = ibd.duplicate(f"{label_str}_rmpeaks")
ibd2.remove_peaks()
ibd2._df = ibd2.cut_and_split_ibd()
of_rmpeaks_ibdne_obj = f"{label_str}_rmpeaks.ibdne.ibdobj.gz"
ibd2.pickle_dump(of_rmpeaks_ibdne_obj)

# link ibdne.jar file
if not Path("ibdne.jar").exists():
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
    mincm=ibdne_mincm,
    minregion=ibdne_minregion,
)
nerunner.run(nthreads=10, mem_gb=20, dry_run=True)

# --- for ibd after removing ibd peaks ------------------------
nerunner = ibdne.IbdNeRunner(
    ibd=ibd2,
    input_folder=".",
    output_folder=".",
    mincm=ibdne_mincm,
    minregion=ibdne_minregion,
)
nerunner.run(nthreads=10, mem_gb=20, dry_run=True)

print(
    f"""
      output files:
        {of_ibddist_obj}
        ibdne.jar
        {idx}_orig.sh
        {idx}_orig.map
        {idx}_orig.ibd.gz
        {idx}_rmpeaks.sh
        {idx}_rmpeaks.map
        {idx}_rmpeaks.ibd.gz
        {of_orig_ibdne_obj}
        {of_rmpeaks_ibdne_obj}
      """
)
