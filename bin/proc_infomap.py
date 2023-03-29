#! /usr/bin/env python3

import ibdutils.utils.ibdutils as ibdutils
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

# parse arguments
pa = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
pa.add_argument("--ibd_files", type=str, nargs=14, required=True)
pa.add_argument("--vcf_files", type=str, nargs=14, required=True)
pa.add_argument("--genome_set_id", type=int, required=True)
args = pa.parse_args()


# read ibd
genome_14_100 = ibdutils.Genome.get_genome("simu_14chr_100cm")
ibd = ibdutils.IBD(genome=genome_14_100, label="orig")
ibd.read_ibd(ibd_fn_lst=args.ibd_files)


# output files:
of_ifm_orig_ibd_obj = f"{args.genome_set_id}_orig.ifm.ibdobj.gz"
of_ifm_rmpeaks_ibd_obj = f"{args.genome_set_id}_rmpeaks.ifm.ibdobj.gz"


# remove highly relatedness samples
mat = ibd.make_ibd_matrix()
unrelated_samples = ibd.get_unrelated_samples(ibdmat=mat)
ibd.subset_ibd_by_samples(subset_samples=unrelated_samples)


# calculate coverage and remove peaks
ibd.calc_ibd_cov()
ibd.find_peaks()

# calculate XiR,s and filter peaks
xirs_df = ibd.calc_xirs(vcf_fn_lst=args.vcf_files, min_maf=0.01)
ibd.filter_peaks_by_xirs(xirs_df)


ibd2 = ibd.duplicate("rmpeak")
ibd2.remove_peaks()

ibd.pickle_dump(of_ifm_orig_ibd_obj)
ibd2.pickle_dump(of_ifm_rmpeaks_ibd_obj)
