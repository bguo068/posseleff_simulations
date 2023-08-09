#! /usr/bin/env python3

from pathlib import Path

import pandas as pd

"""resdir/10665_sp_g080_rep66/ne_output/10665_rmpeaks.ne"""


class DataHandler:
    def __init__(self) -> None:
        self.resdir = "resdir"
        # self.nreps = 100
        self.nreps = 30
        self.orig_rm_modes = ["orig", "rmpeaks"]
        self.sp_defaults = dict(
            seqlen=100 * 15000,
            selpos=int(0.33 * 100 * 15000),
            num_origins=1,
            N=10000,
            h=0.5,
            s=0.3,
            g_sel_start=80,
            r=0.01 / 15_000,
            sim_relatedness=0,
            g_ne_change_start=200,
            N0=1000,
            u=1e-8,
            nsam=1000,  # haploid
        )

        self.sp_sets = dict(
            sp_neu=dict(s=0.0, genome_set_id=10000),
            sp_s01=dict(s=0.1, genome_set_id=10001),
            sp_s02=dict(s=0.2, genome_set_id=10002),
            sp_s03=dict(s=0.3, genome_set_id=10003),
            # sp_g040=dict(g_sel_start=40, genome_set_id=10004),
            # sp_g080=dict(g_sel_start=80, genome_set_id=10005),
            # sp_g120=dict(g_sel_start=120, genome_set_id=10006),
            # sp_o01=dict(num_origins=1, genome_set_id=10007),
            # sp_o03=dict(num_origins=3, genome_set_id=10008),
            # sp_o27=dict(num_origins=27, genome_set_id=10009),
            # sp_rel=dict(sim_relatedness=1, genome_set_id=30000),
        )

        self.mp_defaults = dict(
            seqlen=100 * 15000,
            selpos=int(100 * 15000 * 0.33),
            num_origins=1,
            N=10000,
            h=0.5,
            s=0.3,
            g_sel_start=80,
            r=0.01 / 15_000,
            sim_relatedness=0,
            mig=1e-5,
            sel_mig=0.01,
            npop=5,
            nsam=200,  # haploid
            Tsplit=500,
            u=1e-8,
        )

        self.mp_sets = dict(
            mp_s00=dict(s=0.0, genome_set_id=20000),
            mp_s01=dict(s=0.1, genome_set_id=20001),
            mp_s02=dict(s=0.2, genome_set_id=20002),
            mp_s03=dict(s=0.3, genome_set_id=20003),
            # mp_rel=dict(sim_relatedness=1, genome_set_id=30001),
        )

    def _get_ne_path_single(
        self, grp_label: str, rep_id: int, orig_rm_mode: str
    ) -> str:
        assert grp_label in self.sp_sets.keys()
        assert orig_rm_mode in self.orig_rm_modes
        label = grp_label + f"_rep{rep_id}"
        genome_set_id = self.sp_sets[grp_label]["genome_set_id"] + rep_id * 10
        """resdir/10665_sp_g080_rep66/ne_output/10665_rmpeaks.ne"""
        p = f"{self.resdir}/{genome_set_id}_{label}/ne_output/{genome_set_id}_{orig_rm_mode}.ne"
        assert Path(p).exists()
        return p

    def _get_member_path_single(
        self, grp_label: str, rep_id: int, orig_rm_mode: str
    ) -> str:
        assert grp_label in self.mp_sets.keys()
        assert orig_rm_mode in self.orig_rm_modes
        label = grp_label + f"_rep{rep_id}"
        genome_set_id = self.mp_sets[grp_label]["genome_set_id"] + rep_id * 10
        """resdir/20052_mp_s02_rep5/ifm_output/20052_rmpeaks_member.pq"""
        p = f"{self.resdir}/{genome_set_id}_{label}/ifm_output/{genome_set_id}_{orig_rm_mode}_member.pq"
        assert Path(p).exists()
        return p

    def get_ne_dataframe(self):
        df_lst = []
        for grp_label in self.sp_sets.keys():
            for orig_rm_mode in self.orig_rm_modes:
                for nrep in range(self.nreps):
                    fn = self._get_ne_path_single(
                        grp_label, nrep, orig_rm_mode)
                    df = pd.read_csv(fn, sep="\t")
                    df.columns = ["Gen", "Ne", "L95", "U95"]
                    df["Grp"] = grp_label
                    df["Rep"] = nrep
                    df["OrigRmMode"] = orig_rm_mode
                    df_lst.append(df)

        return pd.concat(df_lst, axis=0)

    def get_member_dataframe(self):
        df_lst = []
        for grp_label in self.mp_sets.keys():
            for orig_rm_mode in self.orig_rm_modes:
                for nrep in range(self.nreps):
                    fn = self._get_member_path_single(
                        grp_label, nrep, orig_rm_mode)
                    df = pd.read_parquet(fn)
                    df["Grp"] = grp_label
                    df["Rep"] = nrep
                    df["OrigRmMode"] = orig_rm_mode
                    df_lst.append(df)

        return pd.concat(df_lst, axis=0)


if __name__ == "__main__":
    dh = DataHandler()

    # ne table
    df = dh.get_ne_dataframe()
    df.to_parquet("ne.pq")

    # comm memebership table
    df = dh.get_member_dataframe()
    df.to_parquet("member.pq")
