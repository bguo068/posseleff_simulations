#! /usr/bin/env python3

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

import ibdutils.utils.ibdutils as ibdutils

# parse arguments
pa = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
pa.add_argument("--ibd_files", type=str, nargs="+", required=True)
pa.add_argument("--vcf_files", type=str, nargs="+", required=True)
pa.add_argument("--genome_set_id", type=int, required=True)
pa.add_argument(
    "--peak_validate_meth", type=str, choices=["xirs", "ihs"], default="xirs"
)
args = pa.parse_args()

# use the number files to indicate the number of chromosome included
# in the genome_set
nchrom = len(args.ibd_files)
assert nchrom == len(args.vcf_files), "ibd /vcf file lists should be of the same length"

# genome
if nchrom == 1:
    genome = ibdutils.Genome.get_genome_simple_simu(
        r=0.01 / 15000, nchroms=1, seqlen_bp_chr=100 * 15000
    )
elif nchrom == 14:
    genome = ibdutils.Genome.get_genome("simu_14chr_100cm")
else:
    raise NotImplemented("nchrom can only be 1 or 14")

# read ibd
ibd = ibdutils.IBD(genome=genome, label="orig")
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
if args.peak_validate_meth == "xirs":
    xirs_df = ibd.calc_xirs(vcf_fn_lst=args.vcf_files, min_maf=0.01)
    ibd.filter_peaks_by_xirs(xirs_df)
elif args.peak_validate_meth == "ihs":
    ibd.calc_ihs(vcf_fn_lst=args.vcf_files, min_maf=0.01)
    ibd.filter_peaks_by_ihs(min_ihs_hits=1)
else:
    raise NotImplementedError(f"{args.peak_validate_meth} is not valid method")


ibd2 = ibd.duplicate("rmpeak")
ibd2.remove_peaks()

ibd.pickle_dump(of_ifm_orig_ibd_obj)
ibd2.pickle_dump(of_ifm_rmpeaks_ibd_obj)
