#! /usr/bin/env python3
import allel
from subprocess import run
import pandas as pd
import argparse


def get_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--tree", type=str, required=True)
    parser.add_argument("--vcf", type=str, required=True)
    parser.add_argument("--chrno", type=int, required=True)
    parser.add_argument("--r", type=float, default=0.01 / 15000)
    parser.add_argument("--seqlen", type=int, default=15000 * 100)
    parser.add_argument("--mincm", type=float, default=2.0)
    parser.add_argument("--call_hmmibd", type=int, choices=[0, 1], default=0)
    # for hmm only
    parser.add_argument("--n", type=int, default=100)
    parser.add_argument("--m", type=int, default=5)
    parser.add_argument("--genome_set_id", type=int, required=True)
    args = parser.parse_args()

    args.bp_per_cm = int(0.01 / args.r)
    args.seqlen_in_cm = args.seqlen * 100 * args.r
    print(args)

    return args


# ------------------- make genetic map ----------------------------
def write_map_file(ofs_map, args):
    with open(ofs_map, "w") as f:
        fields = [f"{args.chrno}", ".", "0", "1"]
        f.write("\t".join(fields) + "\n")
        fields = [
            f"{args.chrno}",
            ".",
            f"{args.seqlen_in_cm}",
            f"{args.seqlen}",
        ]
        f.write("\t".join(fields) + "\n")


# ------------------- prepare hmmibd input ----------------------------
def prep_hmmibd_input(vcf_fn, hmmibd_input_fn="hmm_inputs.txt"):
    calldata = allel.read_vcf(vcf_fn, fields=["samples", "CHROM", "POS", "GT"])
    samples = calldata["samples"]
    gt = calldata["calldata/GT"][:, :, 0]

    df_pos = pd.DataFrame(
        {
            "chrom": calldata["variants/CHROM"],
            "pos": calldata["variants/POS"],
        }
    )
    df_gt = pd.DataFrame(gt, columns=samples)

    df_hmm_inputs = pd.concat([df_pos, df_gt], axis=1)
    df_hmm_inputs.to_csv(hmmibd_input_fn, sep="\t", index=None)


# ------------------- run hmmibd  ----------------------------
def run_hmmibd_and_filter(hmmibd_input_fn, hmmibd_output_prefix, ofs_hmmibd, args):

    # run hmmibd
    cmd = (
        f"hmmIBD -i {hmmibd_input_fn} -o {hmmibd_output_prefix} -n {args.n} -m {args.m}"
    )
    assert run(cmd, shell=True).returncode == 0

    # read raw hmmibd results and keep only segment IBD (remove segment not IBD)
    #
    # Notes from hmmIBD github readme:
    # The file <filename>.hmm.txt contains a list of all segments,
    # A segment is one or more contiguous variant sites in the same state.
    # For each sample pair, the assigned state is either IBD state or not-IBD state.
    # The number of variants covered by the segment are also included.
    # Note that an assigned state of 0 means IBD, while 1 means not-IBD.
    # These segments represent the most probable state assignments.

    df_ibd = pd.read_csv(f"{hmmibd_output_prefix}.hmm.txt", sep="\t")[
        lambda x: (x["different"] == 0)  # keep only IBD segment; exlcude non-IBD
    ]
    df_ibd.columns = ["Id1", "Id2", "Chr", "Start", "End", "Diff", "Nsnp"]

    # update tsk sample name to nunber only names
    df_ibd["Id1"] = df_ibd["Id1"].str.replace("tsk_", "", regex=False)
    df_ibd["Id2"] = df_ibd["Id2"].str.replace("tsk_", "", regex=False)

    # add fake columns
    df_ibd["Ancestor"] = 99999
    df_ibd["Tmrca"] = 100
    df_ibd["HasMutation"] = 0

    sel_cols = ["Id1", "Id2", "Start", "End", "Ancestor", "Tmrca", "HasMutation"]
    sel_rows = ((df_ibd.End - df_ibd.Start) / args.bp_per_cm) >= args.mincm

    df_ibd.loc[sel_rows, sel_cols].to_csv(ofs_hmmibd, sep="\t", index=None)
    # f"{args.chrno}.ibd"


def run_tskibd(ofs_tskibd, args):
    sample_window = int(0.01 * args.bp_per_cm)
    assert (
        run(
            # run tskibd within subfolder to avoid folder containmination
            f"""
                mkdir tskibd_{args.chrno}
                cd tskibd_{args.chrno}
                tskibd {args.chrno} {args.bp_per_cm} {sample_window} {args.mincm} \
                    ../{args.tree}
                cd ..
                mv tskibd_{args.chrno}/{args.chrno}.ibd {ofs_tskibd}
                rm -rf tskibd_{args.chrno}
            """,
            shell=True,
            check=True,
        ).returncode
        == 0
    )


if __name__ == "__main__":
    args = get_args()

    # out file names
    hmmibd_input_fn = f"hmm_inputs_{args.chrno}.txt"
    hmmibd_output_prefix = f"hmm_out_{args.chrno}"
    ofs_map = f"{args.genome_set_id}_{args.chrno}.map"
    ofs_hmmibd = f"{args.genome_set_id}_{args.chrno}_hmmibd.ibd"
    ofs_tskibd = f"{args.genome_set_id}_{args.chrno}_tskibd.ibd"

    if args.call_hmmibd == 1:

        prep_hmmibd_input(args.vcf, hmmibd_input_fn)

        run_hmmibd_and_filter(hmmibd_input_fn, hmmibd_output_prefix, ofs_hmmibd, args)

    write_map_file(ofs_map, args)

    run_tskibd(ofs_tskibd, args)

    print(
        f"""
        output files:
            {ofs_hmmibd} ({'' if args.call_hmmibd==1 else 'Skipped'})
            {ofs_tskibd}
            {ofs_map} """
    )
