#! /usr/bin/env python3

from ibdutils.utils.ibdutils import IBD
import pandas as pd


df_list = []
for raw_gs_id in [10000, 10003, 20000, 20003]:

    label = {
        10000: "sp_neu",
        10003: "sp_s03",
        20000: "mp_s00",
        20003: "mp_s03",
    }[raw_gs_id]

    for rep in range(0, 30):
        gs_id = raw_gs_id + rep * 10

        for pm in ["iqr", "std"]:
            for xc in ["bonferroni", "bonferroni_cm", "fdr_bh"]:

                if label.startswith("sp"):
                    fn = (
                        f"resdir/{gs_id}_{label}_rep{rep}/ibdne_ibd/{gs_id}_orig"
                        ".ibdne.ibdobj.gz"
                    )
                    print(fn)
                else:
                    fn = f"resdir/{gs_id}_{label}_rep{rep}/ifm_input/{gs_id}_orig.ifm.ibdobj.gz"
                    print(fn)

                ibd = IBD.pickle_load(fn)

                # refind peaks and update xirs correction
                ibd._IBD__peaks_df_bk_with_num_xirs_hits = None
                ibd._peaks_df = None
                ibd._peaks_df_bk = None

                ibd.find_peaks(method=pm)
                ibd.calc_xirs([], multitest_correction=xc, only_adj_pvalues=True)
                ibd.filter_peaks_by_xirs(ibd._xirs_df, min_xirs_hits=1)

                xirs = ibd._xirs_df
                if ibd._peaks_df_bk is not None:
                    ibd._peaks_df = ibd._peaks_df_bk.copy()
                else:
                    ibd._peaks_df = pd.DataFrame({})

                out_peaks_df = ibd._IBD__peaks_df_bk_with_num_xirs_hits
                cols = list(ibd._peaks_df.columns)
                out_peaks_df = ibd._peaks_df.merge(
                    out_peaks_df, how="left", on=cols
                ).fillna(0)
                if out_peaks_df.shape[0] > 0:
                    out_peaks_df["GsId"] = gs_id
                    out_peaks_df["Label"] = label
                    out_peaks_df["Rep"] = rep
                    out_peaks_df["PeakFindMeth"] = pm
                    out_peaks_df["XirsCorrMeth"] = xc

                df_list.append(out_peaks_df)

                print(gs_id, label, rep, pm, xc)

pd.concat(df_list, axis=0).to_parquet("num_xirs_hits.pq")
